// Copyright 2018-2026 Nicholas Ingolia
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Annotations of spliced transcripts that contain a coding sequence
//!
//! This module provides data structures for named transcripts and
//! collections of transcripts.

use std::cmp::{max, min};
use std::collections::HashMap;
use std::fmt;
use std::hash::Hash;
use std::io;
use std::num::ParseIntError;
use std::ops::{Deref, Range};

use anyhow::{Context, Result, anyhow, bail, ensure};
use bio::data_structures::annot_map::AnnotMap;
use bio::io::bed;
use bio::io::common::Records;
use bio_types::annot::loc::Loc;
use bio_types::annot::pos::*;
use bio_types::annot::refids::RefIDSet;
use bio_types::annot::spliced::*;
use bio_types::strand::*;

/// Annotation of a transcript as a [`Spliced`] genomic location.
///
/// The transcript is associated with a gene (one gene may have
/// multiple transcripts) and has an optional coding sequence.
///
/// Parameterized over the data type used for identifiers (e.g.,
/// `String`, `Rc<String>`, or `Arc<String>`).
#[derive(Debug, Clone)]
pub struct Transcript<R> {
    gene: R,
    trxname: R,
    loc: Spliced<R, ReqStrand>,
    cds: Option<Range<usize>>,
}

impl<R> Transcript<R> {
    /// Returns the genomic location of the transcript.
    pub fn loc(&self) -> &Spliced<R, ReqStrand> {
        &self.loc
    }

    /// Returns the (optional) coding sequence in transcript
    /// coordinates.
    pub fn cds_range(&self) -> &Option<Range<usize>> {
        &self.cds
    }

    /// Does the transcript contain a coding sequence?
    pub fn is_coding(&self) -> bool {
        self.cds.is_some()
    }

    /// Does the transcript lack a coding sequence?
    pub fn is_noncoding(&self) -> bool {
        self.cds.is_none()
    }

    /// Returns a reference to the gene name.
    pub fn gene_ref(&self) -> &R {
        &self.gene
    }

    /// Returns a reference to the transcript name.
    pub fn trxname_ref(&self) -> &R {
        &self.trxname
    }
}

impl<R: Eq> Transcript<R> {
    /// Gathers transcripts into groups that share the same gene name
    ///
    /// Returns an association list of transcript groups in the form,
    /// ```
    /// [(gene1_name, [gene1_trx1, gene1_trx2, ...]), ...]
    /// ```
    pub fn group_by_gene<'a, I>(trx_iter: I) -> Vec<(&'a R, Vec<&'a Transcript<R>>)>
    where
        I: Iterator<Item = &'a Transcript<R>>,
    {
        let mut groups: Vec<(&'a R, Vec<&'a Transcript<R>>)> = Vec::new();

        for trx in trx_iter {
            let mut done = false;
            {
                if let Some(gene_trxs) = Self::find_mut(groups.as_mut_slice(), trx.gene_ref()) {
                    gene_trxs.push(trx);
                    done = true;
                }
            }

            if !done {
                groups.push((trx.gene_ref(), vec![trx]));
            }
        }

        groups
    }

    /// Find an item by name in an association list
    ///
    /// For an association list `gs = [(key1, value1), ...]`, returns a mutable reference
    /// to the `value` for the first entry where `key == g`, or `None` if no entry
    /// matches `g`.
    fn find_mut<'a, 'b, S: PartialEq, T>(gs: &'b mut [(&'a S, T)], g: &'a S) -> Option<&'b mut T> {
        for (x, xs) in gs.iter_mut() {
            if *x == g {
                return Some(xs);
            }
        }
        return None;
    }
}

impl<R> Transcript<R>
where
    R: Deref<Target = String>,
{
    /// Returns the gene name for the transcript as a `&str`.
    pub fn gene(&self) -> &str {
        self.gene.deref()
    }

    /// Returns the name of the transcript as a `&str`.
    pub fn trxname(&self) -> &str {
        self.trxname.deref()
    }

    /// Returns a triplet of `(reference, start, end)` that encompasses the
    /// genomic location of the transcript.
    ///
    /// This triplet can be used with [`rust_htslib::bam::IndexedReader::fetch()`]
    /// to retrieve all reads overlapping a transcript from an indexed BAM file.
    pub fn fetch_desc(&self) -> (&str, i64, i64) {
        (
            self.loc.refid().deref(),
            self.loc.start() as i64,
            self.loc.start() as i64 + self.loc.length() as i64,
        )
    }
}

impl<R> Transcript<R>
where
    R: Deref<Target = String> + fmt::Debug,
{
    /// Returns a new `Transcript`.
    ///
    /// ```
    /// # use std::error::Error;
    /// # extern crate bio_types;
    /// # extern crate riboprof;
    /// # use std::rc::Rc;
    /// # use bio_types::annot::loc::Loc;
    /// # use bio_types::annot::pos::Pos;
    /// # use bio_types::strand::ReqStrand;
    /// use riboprof::transcript::*;
    /// # fn main() -> Result<(), Box<Error>> {
    /// let gene = Rc::new("ENSG00000245848".to_string());
    /// let trxname = Rc::new("ENST00000498907.2".to_string());
    /// let trx = Transcript::new(gene.clone(), trxname.clone(),
    ///                           "chr19:33299934-33302565(-)".parse()?,
    ///                           Some(150..1227))?;
    /// let trx_pos = Pos::new(trxname.clone(), 507, ReqStrand::Forward);
    /// let chr_pos = trx.loc().pos_outof(&trx_pos).ok_or("Pos out of transcript bounds")?;
    /// assert_eq!(format!("{}", chr_pos), "chr19:33302057(-)");
    /// # Ok(())
    /// # }
    /// ```
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
    ) -> Result<Self> {
        if cds.as_ref().map_or(false, |rng| rng.end <= rng.start) {
            bail!("Invalid CDS range {:?}", cds);
        } else if cds.as_ref().map_or(false, |rng| rng.end > loc.length()) {
            bail!(
                "CDS range {:?} extends beyond end of transcript {:?}",
                cds,
                loc
            );
        } else {
            Ok(Transcript {
                gene: gene,
                trxname: trxname,
                loc: loc,
                cds: cds,
            })
        }
    }
}

impl<R> Transcript<R>
where
    R: Deref<Target = String> + From<String> + Eq + Clone,
{
    /// Construct a transcript from a 12-column BED annotation.
    ///
    /// The gene and transcript name are both taken from the BED name
    /// entry, which is required. The overall transcript annotation is
    /// determined from the BED location and strand information along
    /// with the exon "blocks" location in columns 10 through 12. The
    /// CDS is determined by the "thickStart" and "thickEnd"
    /// entries; if these are equal, then the CDS is `None.`
    ///
    /// ```
    /// # use std::error::Error;
    /// # extern crate bio;
    /// # extern crate bio_types;
    /// # extern crate riboprof;
    /// # use std::rc::Rc;
    /// use bio::io::bed;
    /// # use bio_types::annot::loc::Loc;
    /// # use bio_types::annot::pos::Pos;
    /// # use bio_types::annot::refids::RefIDSet;
    /// # use bio_types::strand::ReqStrand;
    /// use riboprof::transcript::*;
    /// # fn main() -> Result<(), Box<Error>> {
    /// let bed_str = "chr01	87261	87822	YAL030W	0	+	87285	87752	0	2	126,322,	0,239,\n";
    /// let bed = bed::Reader::new(bed_str.as_bytes()).records().next().ok_or("No record")??;
    /// let mut refids: RefIDSet<Rc<String>> = RefIDSet::new();
    /// let trx = Transcript::from_bed12(&bed, &mut refids)?;
    /// let trx_start = Pos::new(trx.trxname().clone(), trx.cds_range().as_ref().ok_or("No CDS")?.start as isize, ReqStrand::Forward);
    /// let chr_start = trx.loc().pos_outof(&trx_start).ok_or("Pos not within transcript")?;
    /// assert_eq!(format!("{}", chr_start), "chr01:87285(+)");
    /// # Ok(())
    /// # }
    /// ```
    /// # Arguments
    ///
    /// `record` is a BED format record containing the annotation information
    ///
    /// `refids` is a table of interned strings used for the gene and
    /// transcript name, along with the reference sequence
    /// (chromosome) name.
    ///
    /// # Errors
    /// The BED information is missing, malformed, or inconsistent.
    pub fn from_bed12(record: &bed::Record, refids: &mut RefIDSet<R>) -> Result<Self> {
        let loc =
            Self::loc_from_bed(record, refids).context("transcript location from BED record")?;
        let cds = Self::cds_from_bed(record, &loc).context("CDS location from BED record")?;
        let name = record
            .name()
            .ok_or_else(|| anyhow!("BED record with no name"))?;

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
    ) -> Result<Spliced<R, ReqStrand>> {
        let block_count = record
            .aux(Self::BLOCK_COUNT_COL)
            .ok_or_else(|| anyhow!("No splicing blocks"))?
            .parse::<usize>()
            .context("Parsing splicing block count")?;

        let block_sizes = record
            .aux(Self::BLOCK_SIZES_COL)
            .ok_or_else(|| anyhow!("No splicing block sizes"))?
            .split_terminator(",")
            .map(|size_str| size_str.parse::<usize>())
            .collect::<Result<Vec<usize>, ParseIntError>>()
            .context("Parsing splicing block lengths")?;

        let block_starts = record
            .aux(Self::BLOCK_STARTS_COL)
            .ok_or_else(|| anyhow!("No splicing block starts"))?
            .split_terminator(",")
            .map(|size_str| size_str.parse::<usize>())
            .collect::<Result<Vec<usize>, ParseIntError>>()
            .context("Parsing splicing block starts")?;

        if block_sizes.len() != block_count || block_starts.len() != block_count {
            bail!(
                "block count = {}, |sizes| = {}, |starts| = {}",
                block_count,
                block_sizes.len(),
                block_starts.len()
            );
        }

        let strand = match record
            .aux(Self::STRAND_COL)
            .ok_or_else(|| anyhow!("No strand"))?
        {
            "+" => ReqStrand::Forward,
            "-" => ReqStrand::Reverse,
            _ => bail!("Bad strand"),
        };

        Spliced::with_lengths_starts(
            refids.intern(record.chrom()),
            record.start() as isize,
            &block_sizes,
            &block_starts,
            strand,
        )
        .context("Splicing error")
    }

    fn cds_from_bed(
        record: &bed::Record,
        loc: &Spliced<R, ReqStrand>,
    ) -> Result<Option<Range<usize>>> {
        if let Some(thick_start) = record
            .aux(Self::THICK_START_COL)
            .map(|start_str| start_str.parse::<usize>())
            .transpose()
            .context("Parsing thickStart")?
        {
            if let Some(thick_end) = record
                .aux(Self::THICK_END_COL)
                .map(|end_str| end_str.parse::<usize>())
                .transpose()
                .context("Parsing thickEnd")?
            {
                if thick_start >= thick_end {
                    return Ok(None);
                }

                // Left-most position within the location
                let left_pos = loc
                    .pos_into(&Pos::new(
                        loc.refid().clone(),
                        thick_start as isize,
                        loc.strand(),
                    ))
                    .ok_or_else(|| anyhow!("thickStart not within annot"))?
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
                    ))
                    .ok_or_else(|| anyhow!("thickEnd-1 not in annot"))?
                    .pos()
                };

                let start = min(left_pos, right_pos) as usize;
                let last = max(left_pos, right_pos) as usize;

                assert!(last >= start);

                Ok(Some(Range {
                    start: start,
                    end: last + 1,
                }))
            } else {
                Ok(None)
            }
        } else {
            Ok(None)
        }
    }
}

/// Returns true when `inner` is compatible with the splicing of `outer`.
///
/// Compatible splicing means that `inner` is on the same strand as
/// `outer`, lies within it, and has a congruent splicing structure.
///
/// ```
/// # extern crate bio_types;
/// # extern crate riboprof;
/// # use std::error::Error;
/// # fn main() -> Result<(),Box<Error>> {
/// use bio_types::annot::spliced::Spliced;
/// use bio_types::strand::ReqStrand;
/// use riboprof::transcript::*;
/// let gene_loc: Spliced<String, ReqStrand> = "chr01:87261-87387;87500-87822(+)".parse()?;
/// let read1_loc = "chr01:87350-87387;87500-87513(+)".parse()?;
/// let read2_loc = "chr01:87464-87513(+)".parse()?;
/// assert!(splice_compatible(&gene_loc, &read1_loc));
/// assert!(!splice_compatible(&gene_loc, &read2_loc));
/// # Ok(())
/// # }
/// ```
pub fn splice_compatible<R: Clone + Eq>(
    outer: &Spliced<R, ReqStrand>,
    inner: &Spliced<R, ReqStrand>,
) -> bool {
    let mut inner_contigs = inner.exon_contigs().into_iter();

    let mut prev_end = if let Some(c0) = inner_contigs.next() {
        let start = if let Some(c0in) = outer.pos_into(&c0.first_pos()) {
            if c0in.strand() != ReqStrand::Forward {
                return false;
            }
            c0in.pos()
        } else {
            return false;
        };

        let end = if let Some(c0in) = outer.pos_into(&c0.last_pos()) {
            c0in.pos()
        } else {
            return false;
        };

        if 1 + end - start != c0.length() as isize {
            return false;
        }

        end
    } else {
        return false;
    };

    for c in inner_contigs {
        let start = if let Some(cin) = outer.pos_into(&c.first_pos()) {
            cin.pos()
        } else {
            return false;
        };

        if start != prev_end + 1 {
            return false;
        }

        let end = if let Some(cin) = outer.pos_into(&c.last_pos()) {
            cin.pos()
        } else {
            return false;
        };

        if 1 + end - start != c.length() as isize {
            return false;
        }

        prev_end = end;
    }

    true
}

/// A nucleotide position within a transcript
pub struct TrxPos<'a, R: 'a> {
    transcript: &'a Transcript<R>,
    pos: usize,
}

impl<'a, R: 'a> TrxPos<'a, R> {
    /// Construct a new transcript position
    pub fn new(transcript: &'a Transcript<R>, pos: usize) -> Result<Self> {
        ensure!(
            pos < transcript.loc().exon_total_length(),
            "Position {} beyond transcript length {}",
            pos,
            transcript.loc().exon_total_length()
        );

        Ok(TrxPos {
            transcript: transcript,
            pos: pos,
        })
    }

    /// Returns the reference transcript for the position
    pub fn transcript(&self) -> &'a Transcript<R> {
        self.transcript
    }

    /// Returns the offset along the transcript
    pub fn pos(&self) -> usize {
        self.pos
    }

    /// Returns the distance from the start of the transcript to the position.
    pub fn offset_from_trx_start(&self) -> isize {
        self.pos as isize
    }

    /// Returns the distance from the end of the transcript to the position.
    ///
    /// This will always be negative.
    pub fn offset_from_trx_end(&self) -> isize {
        self.transcript.loc().exon_total_length() as isize - self.pos as isize
    }

    /// Returns the distance from the start of the coding sequence to the position.
    ///
    /// This will be negative when the position is before the start of the
    /// coding sequence, 0 when it is the first base of the start codon, and positive
    /// otherwise.
    pub fn offset_from_cds_start(&self) -> Option<isize> {
        self.transcript
            .cds_range()
            .as_ref()
            .map(|cds| self.pos as isize - cds.start as isize)
    }

    /// Returns the distance from the end of the coding sequence to the position.
    ///
    /// This will be negative for all positions through the last base of
    /// the stop codon, zero for the first base after the end of the coding
    /// sequence, and positive beyond.
    pub fn offset_from_cds_end(&self) -> Option<isize> {
        self.transcript
            .cds_range()
            .as_ref()
            .map(|cds| self.pos as isize - cds.end as isize)
    }

    pub fn cds_frame(&self) -> Option<usize> {
        self.offset_from_cds_start()
            .map(|off| ((off % 3) + 3) as usize % 3)
    }
}

impl<'a, R: 'a + Eq> TrxPos<'a, R> {
    /// Constructs a transcript position given a transcript and a genomic position.
    ///
    /// If the genomic position `gpos` does not fall within the transcript
    /// and on the same strand, returns `None`.
    pub fn from_genomic_pos<'b>(
        transcript: &'a Transcript<R>,
        gpos: &'b Pos<R, ReqStrand>,
    ) -> Option<Self> {
        transcript
            .loc()
            .pos_into(gpos)
            .map_or(None, |tpos| match tpos.strand() {
                ReqStrand::Forward => Some(tpos.pos()),
                ReqStrand::Reverse => None,
            })
            .map(|pos| Self::new(transcript, pos as usize).expect("bad offset in from_genomic_pos"))
    }
}

impl<'a, R: Eq + Hash> TrxPos<'a, R> {
    /// Transcript positions for a given genomic position, across a collection of transcripts.
    pub fn transcriptome_pos<'b, 'c>(
        tome: &'b Transcriptome<R>,
        gpos: &'c Pos<R, ReqStrand>,
    ) -> impl Iterator<Item = TrxPos<'a, R>>
    where
        'b: 'a,
        'c: 'a,
    {
        tome.find_at_loc(gpos)
            .filter_map(move |trx| Self::from_genomic_pos(trx, gpos))
    }
}

/// Collection of named transcripts
///
/// Transcripts are grouped and indexed for efficient searching by
/// name, position, and insertion order.
pub struct Transcriptome<R>
where
    R: Eq + Hash,
{
    gene_to_trxnames: HashMap<R, Vec<R>>,
    trxname_to_gene: HashMap<R, R>,
    trxname_to_transcript: HashMap<R, Transcript<R>>,
    trxname_by_location: AnnotMap<R, R>,
    trxname_by_index: Vec<R>,
}

impl<R: Eq + Hash> Transcriptome<R> {
    /// Constructs a new, empty transcript collection
    pub fn new() -> Self {
        Transcriptome {
            gene_to_trxnames: HashMap::new(),
            trxname_to_gene: HashMap::new(),
            trxname_to_transcript: HashMap::new(),
            trxname_by_location: AnnotMap::new(),
            trxname_by_index: Vec::new(),
        }
    }

    /// Iterates over all transcript names, in the order of addition
    pub fn trxnames(&self) -> impl Iterator<Item = &R> {
        self.trxname_by_index.iter()
    }

    /// Creates an iterator over all transcripts overlapping a location
    ///
    /// Transcripts will be included when they overlap `loc`
    /// as per [`AnnotMap::find()`](bio::data_structures::annot_map::AnnotMap::find)
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

    /// Returns an iterator over all transcripts, in the order they were added
    pub fn transcripts(&self) -> impl Iterator<Item = &Transcript<R>> {
        self.trxname_by_index
            .iter()
            .map(|t| self.trxname_to_transcript.get(t).unwrap())
    }
}

impl<R> Transcriptome<R>
where
    R: Deref<Target = String> + From<String> + Clone + Hash + Eq,
{
    /// Inserts a transcript into the collection
    ///
    /// # Errors
    /// If the collection already contains a transcript with the same name
    pub fn insert(&mut self, transcript: Transcript<R>) -> Result<R> {
        if self.trxname_to_transcript.contains_key(&transcript.trxname) {
            bail!(
                "Transcript named {:?} already exists",
                transcript.trxname.to_string()
            );
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

        self.trxname_by_index.push(trxname.clone());

        Ok(trxname)
    }

    /// Constructs a transcript collection from BED records
    pub fn new_from_bed<B: io::Read>(
        records: Records<B, bed::Record>,
        refids: &mut RefIDSet<R>,
    ) -> Result<Transcriptome<R>> {
        let mut trxome = Self::new();

        for (lineno, recres) in records.enumerate() {
            // let rec = recres?;
            let rec = recres?;
            let transcript = Transcript::from_bed12(&rec, refids)
                .with_context(|| format!("Parsing line {}:\n{:?}", lineno + 1, rec))?;
            trxome.insert(transcript)?;
        }

        Ok(trxome)
    }
}

#[cfg(test)]
mod tests {
    extern crate csv;

    use super::*;

    use std::rc::*;

    fn record_from_str(recstr: &str) -> bed::Record {
        bed::Reader::new(recstr.as_bytes())
            .records()
            .next()
            .expect("Reading record string")
            .expect("No record read")
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

    fn transcriptome_from_str(bedstr: &str) -> Transcriptome<Rc<String>> {
        let mut refids = RefIDSet::new();
        Transcriptome::new_from_bed(bed::Reader::new(bedstr.as_bytes()).records(), &mut refids)
            .expect("Transcriptome from string")
    }

    fn transcripts_at_pos<R>(tome: &Transcriptome<R>, posstr: &str) -> Vec<String>
    where
        R: Hash + Eq + Deref<Target = String> + From<String>,
    {
        let pos: Pos<R, ReqStrand> = posstr.parse().expect("Parsing position");
        let mut trxs: Vec<String> = tome
            .find_at_loc(&pos)
            .map(|trx| trx.trxname().to_string())
            .collect();
        trxs.sort();
        trxs
    }

    #[test]
    fn transcriptome_find() {
        let beds = "\
chr01	1000	2000	AAA	0	+	1200	1800	0	1	1000,	0,
chr01	1900	2100	BBB	0	+	1950	2050	0	1	200,	0,
chr02	1500	2500	CCC	0	+	1600	2400	0	1	1000,	0,
chr02	2100	2600	DDD	0	-	2200	2500	0	1	500,	0,
chr03	500	1500	EEE	0	+	600	1200	0	2	250,450	0,550
";
        let tome = transcriptome_from_str(&beds);

        let none: Vec<String> = vec![];

        assert_eq!(transcripts_at_pos(&tome, "chr01:750(+)"), none);
        assert_eq!(transcripts_at_pos(&tome, "chr01:999(+)"), none);
        assert_eq!(transcripts_at_pos(&tome, "chr01:1000(+)"), vec!["AAA"]);
        assert_eq!(transcripts_at_pos(&tome, "chr01:1234(+)"), vec!["AAA"]);
        assert_eq!(
            transcripts_at_pos(&tome, "chr01:1950(+)"),
            vec!["AAA", "BBB"]
        );
        assert_eq!(
            transcripts_at_pos(&tome, "chr01:1999(+)"),
            vec!["AAA", "BBB"]
        );
        assert_eq!(transcripts_at_pos(&tome, "chr01:2000(+)"), vec!["BBB"]);
        assert_eq!(transcripts_at_pos(&tome, "chr01:2050(+)"), vec!["BBB"]);
        assert_eq!(transcripts_at_pos(&tome, "chr01:2099(+)"), vec!["BBB"]);
        assert_eq!(transcripts_at_pos(&tome, "chr01:2100(+)"), none);
        assert_eq!(transcripts_at_pos(&tome, "chr01:2101(+)"), none);

        assert_eq!(transcripts_at_pos(&tome, "chr02:2000(+)"), vec!["CCC"]);

        assert_eq!(transcripts_at_pos(&tome, "chr03:550(+)"), vec!["EEE"]);
        assert_eq!(transcripts_at_pos(&tome, "chr03:850(+)"), vec!["EEE"]);
        assert_eq!(transcripts_at_pos(&tome, "chr03:1450(+)"), vec!["EEE"]);
    }

    fn make_spliced(s: &str) -> Spliced<String, ReqStrand> {
        s.parse().expect("Parsing spliced")
    }

    #[test]
    fn test_splice_compatible() {
        let a = make_spliced("chr01:1000-1500;2000-2500;3000-3500(+)");
        assert!(splice_compatible(&a, &make_spliced("chr01:1100-1200(+)")));
        assert!(!splice_compatible(&a, &make_spliced("chr01:1100-1200(-)")));
        assert!(!splice_compatible(&a, &make_spliced("chr02:1100-1200(+)")));

        assert!(splice_compatible(&a, &make_spliced("chr01:1000-1200(+)")));
        assert!(!splice_compatible(&a, &make_spliced("chr01:999-1200(+)")));
        assert!(splice_compatible(&a, &make_spliced("chr01:1300-1500(+)")));
        assert!(!splice_compatible(&a, &make_spliced("chr01:1300-1501(+)")));

        assert!(splice_compatible(&a, &make_spliced("chr01:2000-2200(+)")));
        assert!(!splice_compatible(&a, &make_spliced("chr01:1999-2200(+)")));
        assert!(splice_compatible(&a, &make_spliced("chr01:2300-2500(+)")));
        assert!(!splice_compatible(&a, &make_spliced("chr01:2300-2501(+)")));

        assert!(splice_compatible(
            &a,
            &make_spliced("chr01:1200-1500;2000-2200(+)")
        ));
        assert!(!splice_compatible(
            &a,
            &make_spliced("chr01:1200-1501;2000-2200(+)")
        ));
        assert!(!splice_compatible(
            &a,
            &make_spliced("chr01:1200-1499;2000-2200(+)")
        ));
        assert!(!splice_compatible(
            &a,
            &make_spliced("chr01:1200-1500;2001-2200(+)")
        ));
        assert!(!splice_compatible(
            &a,
            &make_spliced("chr01:1200-1500;1999-2200(+)")
        ));
        assert!(!splice_compatible(
            &a,
            &make_spliced("chr01:1200-1500;1700-1800;2000-2200(+)")
        ));
        assert!(!splice_compatible(
            &a,
            &make_spliced("chr01:1200-1300;1400-1500;2000-2200(+)")
        ));

        assert!(splice_compatible(
            &a,
            &make_spliced("chr01:1200-1500;2000-2500;3000-3200(+)")
        ));
        assert!(splice_compatible(
            &a,
            &make_spliced("chr01:2200-2500;3000-3500(+)")
        ));
        assert!(splice_compatible(
            &a,
            &make_spliced("chr01:2200-2500;3000-3499(+)")
        ));
        assert!(!splice_compatible(
            &a,
            &make_spliced("chr01:2200-2500;3000-3501(+)")
        ));

        let b = make_spliced("chr01:1000-1500;2000-2500;3000-3500(-)");
        assert!(splice_compatible(&b, &make_spliced("chr01:1100-1200(-)")));
        assert!(!splice_compatible(&b, &make_spliced("chr01:1100-1200(+)")));
        assert!(!splice_compatible(&b, &make_spliced("chr02:1100-1200(-)")));

        assert!(splice_compatible(&b, &make_spliced("chr01:1000-1200(-)")));
        assert!(!splice_compatible(&b, &make_spliced("chr01:999-1200(-)")));
        assert!(splice_compatible(&b, &make_spliced("chr01:1300-1500(-)")));
        assert!(!splice_compatible(&b, &make_spliced("chr01:1300-1501(-)")));

        assert!(splice_compatible(&b, &make_spliced("chr01:2000-2200(-)")));
        assert!(!splice_compatible(&b, &make_spliced("chr01:1999-2200(-)")));
        assert!(splice_compatible(&b, &make_spliced("chr01:2300-2500(-)")));
        assert!(!splice_compatible(&b, &make_spliced("chr01:2300-2501(-)")));

        assert!(splice_compatible(
            &b,
            &make_spliced("chr01:1200-1500;2000-2200(-)")
        ));
        assert!(!splice_compatible(
            &b,
            &make_spliced("chr01:1200-1501;2000-2200(-)")
        ));
        assert!(!splice_compatible(
            &b,
            &make_spliced("chr01:1200-1499;2000-2200(-)")
        ));
        assert!(!splice_compatible(
            &b,
            &make_spliced("chr01:1200-1500;2001-2200(-)")
        ));
        assert!(!splice_compatible(
            &b,
            &make_spliced("chr01:1200-1500;1999-2200(-)")
        ));
        assert!(!splice_compatible(
            &b,
            &make_spliced("chr01:1200-1500;1700-1800;2000-2200(-)")
        ));
        assert!(!splice_compatible(
            &b,
            &make_spliced("chr01:1200-1300;1400-1500;2000-2200(-)")
        ));

        assert!(splice_compatible(
            &b,
            &make_spliced("chr01:1200-1500;2000-2500;3000-3200(-)")
        ));
        assert!(splice_compatible(
            &b,
            &make_spliced("chr01:2200-2500;3000-3500(-)")
        ));
        assert!(splice_compatible(
            &b,
            &make_spliced("chr01:2200-2500;3000-3499(-)")
        ));
        assert!(!splice_compatible(
            &b,
            &make_spliced("chr01:2200-2500;3000-3501(-)")
        ));
    }
}
