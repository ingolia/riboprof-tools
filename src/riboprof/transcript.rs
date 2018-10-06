use std::cmp::{max, min};
use std::collections::HashMap;
use std::error::Error;
use std::fmt;
use std::hash::Hash;
use std::io;
use std::num::ParseIntError;
use std::ops::{Deref, Range};

use bio::data_structures::annot_map::AnnotMap;
use bio::io::bed;
use bio_types::annot::loc::Loc;
use bio_types::annot::pos::*;
use bio_types::annot::refids::RefIDSet;
use bio_types::annot::spliced::*;
use bio_types::strand::*;

use failure;

/// Annotation of a transcript as a `Spliced` `annot` location.
///
/// The transcript is associated with a gene (one gene may have
/// multiple transcripts) and has an optional coding sequence.
///
/// Parameterized over the data type used for identifiers (e.g.,
/// `String`, `Rc<String>`, or `Arc<String>`).
pub struct Transcript<R> {
    gene: R,
    trxname: R,
    loc: Spliced<R, ReqStrand>,
    cds: Option<Range<usize>>,
}

impl<R> Transcript<R> {
    /// Returns the spliced location of the transcript.
    pub fn loc(&self) -> &Spliced<R, ReqStrand> {
        &self.loc
    }

    /// Returns the (optional) coding sequence in transcript
    /// coordinates.
    pub fn cds_range(&self) -> &Option<Range<usize>> {
        &self.cds
    }
}

impl<R> Transcript<R>
where
    R: Deref<Target = String> + fmt::Debug,
{
    /// Returns a new `Transcript`.
    ///
    /// # Arguments
    ///
    /// * `gene` is the gene identifier for the annotation
    /// * `trxname` is the transcript identifier for the annotation
    /// * `loc` is the transcript annotation
    /// * `cds` is an optional coding sequence position, as a range of loc.
    ///
    /// # Errors
    ///
    /// An error variant is returned if a CDS is given and the `end`
    /// of the range is beyond `loc.length()` (or less than `start`).
    pub fn new(
        gene: R,
        trxname: R,
        loc: Spliced<R, ReqStrand>,
        cds: Option<Range<usize>>,
    ) -> Result<Self, TrxError> {
        if cds.as_ref().map_or(false, |rng| rng.end <= rng.start) {
            Err(TrxError::Cds(format!("Invalid CDS range {:?}", cds)))
        } else if cds.as_ref().map_or(false, |rng| rng.end > loc.length()) {
            Err(TrxError::Cds(format!(
                "CDS range {:?} extends beyond end of transcript {:?}",
                cds, loc
            )))
        } else {
            Ok(Transcript {
                gene: gene,
                trxname: trxname,
                loc: loc,
                cds: cds,
            })
        }
    }

    /// Returns the gene name for the transcript
    pub fn gene(&self) -> &str {
        &self.gene
    }

    /// Returns the name of the transcript
    pub fn trxname(&self) -> &str {
        &self.trxname
    }
}

impl<R> Transcript<R>
where
    R: Deref<Target = String> + From<String> + Eq + Clone,
{
    /// Construct a transcirpt from a 12-column BED annotation.
    ///
    /// The gene and transcript name are both taken from the BED name
    /// entry, which is required. The overall transcript annotation is
    /// determined from the BED location and strand information along
    /// with the exon "blocks" location in columns 10 through 12. The
    /// CDS is determined by the "thickStart" and "thickEnd"
    /// entries; if these are equal, then the CDS is `None.`
    ///
    /// # Arguments
    ///
    /// `record` is a BED format record containing the annotation information
    ///
    /// `refids` is a table of interned strings used for the gene and
    /// transcript name, along with the reference sequence
    /// (chromosome) name.
    ///
    /// # Errors
    ///
    /// An error variant is returned when required information is
    /// missing, unparseable, or inconsistent.
    pub fn from_bed12(record: &bed::Record, refids: &mut RefIDSet<R>) -> Result<Self, TrxError> {
        let loc = Self::loc_from_bed(record, refids)?;
        let cds = Self::cds_from_bed(record, &loc)?;
        let name = record
            .name()
            .ok_or_else(|| TrxError::bed(record, "No name"))?;

        Ok(Transcript {
            gene: refids.intern(name),
            trxname: refids.intern(name),
            loc: loc,
            cds: cds,
        })
    }

    const STRAND_COL: usize = 5;
    const THICK_START_COL: usize = 6;
    const THICK_END_COL: usize = 7;
    const BLOCK_COUNT_COL: usize = 9;
    const BLOCK_SIZES_COL: usize = 10;
    const BLOCK_STARTS_COL: usize = 11;

    fn loc_from_bed(
        record: &bed::Record,
        refids: &mut RefIDSet<R>,
    ) -> Result<Spliced<R, ReqStrand>, TrxError> {
        let block_count = record
            .aux(Self::BLOCK_COUNT_COL)
            .ok_or_else(|| TrxError::bed(record, "No splicing blocks"))?
            .parse::<usize>()
            .map_err(|err| TrxError::bed_parse(record, "Bad block count", err))?;

        let block_sizes = record
            .aux(Self::BLOCK_SIZES_COL)
            .ok_or_else(|| TrxError::bed(record, "No splicing block sizes"))?
            .split_terminator(",")
            .map(|size_str| size_str.parse::<usize>())
            .collect::<Result<Vec<usize>, ParseIntError>>()
            .map_err(|err| TrxError::bed_parse(record, "Bad block sizes", err))?;

        let block_starts = record
            .aux(Self::BLOCK_STARTS_COL)
            .ok_or_else(|| TrxError::bed(record, "No splicing block starts"))?
            .split_terminator(",")
            .map(|size_str| size_str.parse::<usize>())
            .collect::<Result<Vec<usize>, ParseIntError>>()
            .map_err(|err| TrxError::bed_parse(record, "Bad block starts", err))?;

        if block_sizes.len() != block_count || block_starts.len() != block_count {
            return Err(TrxError::bed(
                record,
                &format!(
                    "block count = {}, |sizes| = {}, |starts| = {}",
                    block_count,
                    block_sizes.len(),
                    block_starts.len()
                ),
            ));
        }

        let strand = match record
            .aux(Self::STRAND_COL)
            .ok_or_else(|| TrxError::bed(record, "No strand"))?
        {
            "+" => Ok(ReqStrand::Forward),
            "-" => Ok(ReqStrand::Reverse),
            _ => Err(TrxError::bed(record, "Bad strand")),
        }?;

        Spliced::with_lengths_starts(
            refids.intern(record.chrom()),
            record.start() as isize,
            &block_sizes,
            &block_starts,
            strand,
        ).map_err(|err| {
            TrxError::BedSplicing(format!("Splicing error on record {:?}", record), err)
        })
    }

    fn cds_from_bed(
        record: &bed::Record,
        loc: &Spliced<R, ReqStrand>,
    ) -> Result<Option<Range<usize>>, TrxError> {
        let thick_start = match record.aux(Self::THICK_START_COL) {
            Some(start_str) => start_str
                .parse::<usize>()
                .map_err(|err| TrxError::bed_parse(record, "thickStart", err))?,
            None => return Ok(None),
        };

        let thick_end = match record.aux(Self::THICK_END_COL) {
            Some(end_str) => end_str
                .parse::<usize>()
                .map_err(|err| TrxError::bed_parse(record, "thickEnd", err))?,
            None => return Ok(None),
        };

        if thick_start >= thick_end {
            return Ok(None);
        }

        // Left-most position within the location
        let left_pos =
            loc.pos_into(&Pos::new(
                loc.refid().clone(),
                thick_start as isize,
                loc.strand(),
            )).ok_or_else(|| TrxError::bed(record, "thickStart not in annot"))?
                .pos();

        // Right-most position _within_ the location
        // When CDS extends to right edge of transcript, then thickEnd is _outside_
        let right_pos = if let Some(pos) = loc.pos_into(&Pos::new(
            loc.refid().clone(),
            thick_end as isize,
            loc.strand(),
        )) {
            match loc.strand() {
                ReqStrand::Forward => pos.pos() - 1,
                ReqStrand::Reverse => pos.pos() + 1,
            }
        } else {
            loc.pos_into(&Pos::new(
                loc.refid().clone(),
                thick_end as isize - 1,
                loc.strand(),
            )).ok_or_else(|| TrxError::bed(record, "thickEnd-1 not in annot"))?
                .pos()
        };

        let start = min(left_pos, right_pos) as usize;
        let last = max(left_pos, right_pos) as usize;

        assert!(last >= start);

        Ok(Some(Range {
            start: start,
            end: last + 1,
        }))
    }
}

pub struct Transcriptome<R>
where
    R: Eq + Hash,
{
    gene_to_trxnames: HashMap<R, Vec<R>>,
    trxname_to_gene: HashMap<R, R>,
    trxname_to_transcript: HashMap<R, Transcript<R>>,
    trxname_by_location: AnnotMap<R, R>,
}

impl<R: Eq + Hash> Transcriptome<R> {
    pub fn new() -> Self {
        Transcriptome {
            gene_to_trxnames: HashMap::new(),
            trxname_to_gene: HashMap::new(),
            trxname_to_transcript: HashMap::new(),
            trxname_by_location: AnnotMap::new(),
        }
    }

    pub fn trxnames(&self) -> impl Iterator<Item = &R> {
        self.trxname_to_transcript.keys()
    }

    pub fn find_at_loc<'a: 'c, 'b: 'c, 'c, L: Loc<RefID = R>>(
        &'a self,
        loc: &'b L,
    ) -> impl Iterator<Item = &'c Transcript<R>> {
        self.trxname_by_location.find(loc).map(move |ent| {
            self.trxname_to_transcript
                .get(ent.data())
                .expect("transcript missing from map")
        })
    }
}

impl<R> Transcriptome<R>
where
    R: Deref<Target = String> + From<String> + Clone + Hash + Eq,
{
    pub fn insert(&mut self, transcript: Transcript<R>) -> Result<R, TrxError> {
        if self.trxname_to_transcript.contains_key(&transcript.trxname) {
            return Err(TrxError::TrxExists(transcript.trxname.to_string()));
        }

        let trxname = transcript.trxname.clone();

        self.trxname_to_gene
            .insert(trxname.clone(), transcript.gene.clone());
        self.gene_to_trxnames
            .entry(transcript.gene.clone())
            .or_insert(vec![])
            .push(trxname.clone());
        self.trxname_by_location
            .insert_at(trxname.clone(), &transcript.loc);

        self.trxname_to_transcript
            .insert(trxname.clone(), transcript);

        Ok(trxname)
    }

    pub fn new_from_bed<B: io::Read>(
        records: bed::Records<B>,
        refids: &mut RefIDSet<R>,
    ) -> Result<Transcriptome<R>, TrxError> {
        let mut trxome = Self::new();

        for recres in records {
            let rec = recres.map_err(|err| TrxError::BedRead(err.into()))?;
            let transcript = Transcript::from_bed12(&rec, refids)?;
            trxome.insert(transcript)?;
        }

        Ok(trxome)
    }
}

#[derive(Debug)]
pub enum TrxError {
    Bed(String),
    BedParse(String, ParseIntError),
    BedRead(failure::Error),
    BedSplicing(String, SplicingError),
    Cds(String),
    TrxExists(String),
}

impl TrxError {
    fn bed(record: &bed::Record, message: &str) -> TrxError {
        TrxError::Bed(format!("{} converting BED record {:?}", message, record))
    }

    fn bed_parse(record: &bed::Record, message: &str, parse_error: ParseIntError) -> TrxError {
        TrxError::BedParse(
            format!("{} parsing BED record {:?}", message, record),
            parse_error,
        )
    }
}

impl Error for TrxError {}

impl fmt::Display for TrxError {
    fn fmt(&self, f: &mut fmt::Formatter) -> Result<(), fmt::Error> {
        match self {
            TrxError::Bed(msg) => write!(f, "BED record to transcript: {}", msg),
            TrxError::BedParse(msg, err) => write!(
                f,
                "BED record to transcript: {}: parsing error {}",
                msg, err
            ),
            TrxError::BedRead(err) => write!(f, "Reading BED records: {}", err),
            TrxError::BedSplicing(msg, err) => write!(
                f,
                "BED record to transcript: {}: splicing error {}",
                msg, err
            ),
            TrxError::Cds(msg) => write!(f, "CDS on transcript: {}", msg),
            TrxError::TrxExists(trx) => write!(f, "Transcript already exists: {}", trx),
        }
    }
}

#[cfg(test)]
mod tests {
    extern crate csv;

    use super::*;

    use std::cell::*;
    use std::ops::*;
    use std::rc::*;

    use self::csv::Error;

    struct TestWriter {
        dest: Rc<RefCell<Vec<u8>>>,
    }

    impl io::Write for TestWriter {
        fn write(&mut self, buf: &[u8]) -> Result<usize, io::Error> {
            self.dest.borrow_mut().append(&mut buf.to_vec());
            Ok(buf.len())
        }

        fn flush(&mut self) -> Result<(), io::Error> {
            Ok(())
        }
    }

    fn record_from_str(recstr: &str) -> bed::Record {
        bed::Reader::new(recstr.as_bytes())
            .records()
            .next()
            .expect("Reading record string")
            .expect("No record read")
        //            .collect::<Result<Vec<bed::Record>, csv::Error>>()
        //            .expect("Reading record string")
    }

    fn transcript_from_str(recstr: &str) -> Transcript<Rc<String>> {
        let rec = record_from_str(recstr);
        let mut refids: RefIDSet<Rc<String>> = RefIDSet::new();
        Transcript::from_bed12(&rec, &mut refids).expect("Converting to transcript")
    }

    fn no_transcript_from_str(recstr: &str) -> bool {
        let rec = record_from_str(recstr);
        let mut refids: RefIDSet<Rc<String>> = RefIDSet::new();
        Transcript::from_bed12(&rec, &mut refids).is_err()
    }

    #[test]
    fn gene_1exon_fwd() {
        let recstr = "chr01	334	649	YAL069W	0	+	334	649	0	1	315,	0,\n";
        let trx = transcript_from_str(recstr);
        assert_eq!(trx.gene(), "YAL069W");
        assert_eq!(trx.loc().to_string(), "chr01:334-649(+)");
        assert_eq!(trx.cds_range(), &Some(0..315));
    }

    #[test]
    fn gene_1exon_rev() {
        let recstr = "chr01	1806	2169	YAL068C	0	-	1806	2169	0	1	363,	0,\n";
        let trx = transcript_from_str(recstr);
        assert_eq!(trx.gene(), "YAL068C");
        assert_eq!(trx.loc().to_string(), "chr01:1806-2169(-)");
        assert_eq!(trx.cds_range(), &Some(0..363));
    }

    #[test]
    fn gene_1exon_fwd_cds() {
        let recstr = "chr01	33364	34785	YAL061W	0	+	33447	34701	0	1	1421,	0,\n";
        let trx = transcript_from_str(recstr);
        assert_eq!(trx.gene(), "YAL061W");
        assert_eq!(trx.loc().to_string(), "chr01:33364-34785(+)");
        assert_eq!(trx.cds_range(), &Some(83..1337));
        // CDS alternatives
        let recstr = "chr01	33364	34785	YAL061W	0	+	33364	34701	0	1	1421,	0,\n";
        assert_eq!(transcript_from_str(recstr).cds_range(), &Some(0..1337));
        let recstr = "chr01	33364	34785	YAL061W	0	+	33365	34701	0	1	1421,	0,\n";
        assert_eq!(transcript_from_str(recstr).cds_range(), &Some(1..1337));
        let recstr = "chr01	33364	34785	YAL061W	0	+	33363	34701	0	1	1421,	0,\n";
        assert!(no_transcript_from_str(recstr));
        let recstr = "chr01	33364	34785	YAL061W	0	+	33447	34785	0	1	1421,	0,\n";
        assert_eq!(transcript_from_str(recstr).cds_range(), &Some(83..1421));
        let recstr = "chr01	33364	34785	YAL061W	0	+	33447	34784	0	1	1421,	0,\n";
        assert_eq!(transcript_from_str(recstr).cds_range(), &Some(83..1420));
        let recstr = "chr01	33364	34785	YAL061W	0	+	33447	34786	0	1	1421,	0,\n";
        assert!(no_transcript_from_str(recstr));
    }

    #[test]
    fn gene_1exon_rev_cds() {
        let recstr = "chr01	51775	52696	YAL049C	0	-	51854	52595	0	1	921,	0,\n";
        let trx = transcript_from_str(&recstr);
        assert_eq!(trx.gene(), "YAL049C");
        assert_eq!(trx.loc().to_string(), "chr01:51775-52696(-)");
        assert_eq!(trx.cds_range(), &Some(101..842));
        // CDS alternatives
        let recstr = "chr01	51775	52696	YAL049C	0	-	51854	52696	0	1	921,	0,\n";
        assert_eq!(transcript_from_str(recstr).cds_range(), &Some(0..842));
        let recstr = "chr01	51775	52696	YAL049C	0	-	51854	52695	0	1	921,	0,\n";
        assert_eq!(transcript_from_str(recstr).cds_range(), &Some(1..842));
        let recstr = "chr01	51775	52696	YAL049C	0	-	51854	52697	0	1	921,	0,\n";
        assert!(no_transcript_from_str(recstr));
        let recstr = "chr01	51775	52696	YAL049C	0	-	51775	52595	0	1	921,	0,\n";
        assert_eq!(transcript_from_str(recstr).cds_range(), &Some(101..921));
        let recstr = "chr01	51775	52696	YAL049C	0	-	51776	52595	0	1	921,	0,\n";
        assert_eq!(transcript_from_str(recstr).cds_range(), &Some(101..920));
        let recstr = "chr01	51775	52696	YAL049C	0	-	51774	52595	0	1	921,	0,\n";
        assert!(no_transcript_from_str(recstr));
    }

    #[test]
    fn gene_2exon_fwd() {
        let recstr = "chr01	87261	87822	YAL030W	0	+	87285	87752	0	2	126,322,	0,239,\n";
        let trx = transcript_from_str(&recstr);
        assert_eq!(trx.gene(), "YAL030W");
        assert_eq!(trx.loc().to_string(), "chr01:87261-87387;87500-87822(+)");
        assert_eq!(trx.cds_range(), &Some(24..378));
    }

    #[test]
    fn gene_2exon_rev() {
        let recstr = "chr02	2906	5009	YBL111C	0	-	2906	5009	0	2	1210,794,	0,1309,\n";
        let trx = transcript_from_str(&recstr);
        assert_eq!(trx.gene(), "YBL111C");
        assert_eq!(trx.loc().to_string(), "chr02:2906-4116;4215-5009(-)");
        assert_eq!(trx.cds_range(), &Some(0..2004));

        let recstr = "chr02	59630	60828	YBL087C	0	-	59821	60739	0	2	563,131,	0,1067,\n";
        let trx = transcript_from_str(&recstr);
        assert_eq!(trx.gene(), "YBL087C");
        assert_eq!(trx.loc().to_string(), "chr02:59630-60193;60697-60828(-)");
        assert_eq!(trx.cds_range(), &Some(89..503));
    }
}
