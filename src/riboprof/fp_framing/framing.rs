use std::ops::Range;
use std::rc::Rc;

use failure;

use bio_types::annot::loc::Loc;
//use bio_types::annot::pos::Pos;
use bio_types::annot::spliced::Spliced;
use bio_types::strand::*;
use rust_htslib::bam;

use bam_utils::*;
//use codon_assign::*;
use transcript::*;

pub fn record_framing(
    trxome: &Transcriptome<Rc<String>>,
    tids: &Tids<Rc<String>>,
    rec: &bam::Record,
    lengths: &Range<usize>,
    cdsbody: &(isize, isize),
    count_multi: bool,
) -> Result<BamFrameResult, failure::Error> {
    if !(is_single_hit(rec) || (count_multi && is_first_hit(rec))) {
        return Ok(BamFrameResult::MultiHit);
    }

    if let Some(fp) = bam_to_spliced(tids, &rec)? {
        let fp_len = fp.exon_total_length();

        if fp_len < lengths.start {
            return Ok(BamFrameResult::TooShort);
        } else if fp_len > lengths.end {
            return Ok(BamFrameResult::TooLong);
        }

        let ffr = footprint_framing(trxome, &fp, cdsbody);
        Ok(BamFrameResult::Fp(ffr))
    } else {
        Ok(BamFrameResult::NoHit)
    }
}

pub fn is_single_hit(rec: &bam::Record) -> bool {
    if let Some(bam::record::Aux::Integer(nh)) = rec.aux(b"NH") {
        nh == 1
    } else {
        true
    }
}

pub fn is_first_hit(rec: &bam::Record) -> bool {
    rec.aux(b"HI") == Some(bam::record::Aux::Integer(1))
}

pub enum BamFrameResult {
    NoHit,
    MultiHit,
    TooShort,
    TooLong,
    Fp(FpFrameResult),
}

impl BamFrameResult {
    pub fn aux(&self) -> Vec<u8> {
        match self {
            BamFrameResult::NoHit => b"BamNoHit".to_vec(),
            BamFrameResult::MultiHit => b"BamMultiHit".to_vec(),
            BamFrameResult::TooShort => b"BamTooShort".to_vec(),
            BamFrameResult::TooLong => b"BamTooLong".to_vec(),
            BamFrameResult::Fp(ffr) => ffr.aux(),
        }
    }
}

pub fn footprint_framing(
    trxome: &Transcriptome<Rc<String>>,
    fp: &Spliced<Rc<String>, ReqStrand>,
    cdsbody: &(isize, isize),
) -> FpFrameResult {
    let gene_sets = Transcript::group_by_gene(trxome.find_at_loc(fp));

    if gene_sets.len() > 1 {
        let is_coding: Vec<bool> = gene_sets
            .values()
            .map(|trxs| trxs.iter().any(|trx| trx.is_coding()))
            .collect();

        if is_coding.iter().all(|coding| *coding) {
            FpFrameResult::MultiCoding
        } else if is_coding.iter().any(|coding| *coding) {
            FpFrameResult::NoncodingOverlap
        } else {
            FpFrameResult::NoncodingOnly
        }
    } else if let Some((_gene, trxs)) = gene_sets.into_iter().next() {
        let coding_trxs: Vec<&Transcript<Rc<String>>> =
            trxs.into_iter().filter(|trx| trx.is_coding()).collect();

        if coding_trxs.is_empty() {
            FpFrameResult::NoncodingOnly
        } else {
            FpFrameResult::Gene(gene_framing(cdsbody, coding_trxs.as_slice(), fp))
        }
    } else {
        // gene_sets is empty
        FpFrameResult::NoGene
    }
}

pub enum FpFrameResult {
    Gene(GeneFrameResult),
    NoGene,
    NoncodingOnly,
    NoncodingOverlap,
    MultiCoding,
}

impl FpFrameResult {
    pub fn aux(&self) -> Vec<u8> {
        match self {
            FpFrameResult::Gene(gfr) => gfr.aux(),
            FpFrameResult::NoGene => "NoGene".to_string().into_bytes(),
            FpFrameResult::NoncodingOnly => "NoncodingOnly".to_string().into_bytes(),
            FpFrameResult::NoncodingOverlap => "NoncodingOverlap".to_string().into_bytes(),
            FpFrameResult::MultiCoding => "MultiCoding".to_string().into_bytes(),
        }
    }
}

/// Result from framing analysis for a footprint, relative to a gene.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GeneFrameResult {
    Good(GeneFraming),
    NoCompatible,
    Ambig,
}

impl GeneFrameResult {
    pub fn aux(&self) -> Vec<u8> {
        match self {
            GeneFrameResult::Good(gf) => gf.aux(),
            GeneFrameResult::NoCompatible => "NoCompatible".to_string().into_bytes(),
            GeneFrameResult::Ambig => "AmbigFrame".to_string().into_bytes(),
        }
    }
}

/// Framing information for a footprint, relative to a
/// gene. This includes the offsets relative to the CDS start and end,
/// as well as the reading frame position.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GeneFraming {
    gene: Rc<String>,
    vs_cds_start: Option<isize>,
    vs_cds_end: Option<isize>,
    frame: Option<usize>,
    fp_length: usize,
}

impl GeneFraming {
    pub fn aux(&self) -> Vec<u8> {
        let mut buf = self.gene.to_string();
        buf += &Self::opt_to_str(&self.vs_cds_start);
        buf += &Self::opt_to_str(&self.vs_cds_end);
        buf += &Self::opt_to_str(&self.frame.as_ref().map(|x| *x as isize));
        buf.into_bytes()
    }

    fn opt_to_str(opt: &Option<isize>) -> String {
        opt.map_or_else(|| "/*".to_string(), |x| format!("/{}", x))
    }

    pub fn vs_cds_start(&self) -> Option<isize> {
        self.vs_cds_start
    }
    pub fn vs_cds_end(&self) -> Option<isize> {
        self.vs_cds_end
    }
    pub fn frame(&self) -> Option<usize> {
        self.frame
    }
    pub fn fp_length(&self) -> usize {
        self.fp_length
    }
}

/// Computes framing information for a footprint, relative to a gene
/// defined by one or more transcripts. Framing in
pub fn gene_framing<'a>(
    cdsbody: &(isize, isize),
    trxs: &[&'a Transcript<Rc<String>>],
    fp: &Spliced<Rc<String>, ReqStrand>,
) -> GeneFrameResult {
    let gene = if trxs.len() == 0 {
        return GeneFrameResult::NoCompatible;
    } else {
        trxs[0].gene_ref().clone()
    };

    let fp_length = fp.exon_total_length();

    let termini: Vec<TrxPos<'a, Rc<String>>> = trxs
        .into_iter()
        .filter_map(move |trx| fp_into_transcript(fp, trx))
        .collect();

    if termini.is_empty() {
        GeneFrameResult::NoCompatible
    } else {
        let vs_cds_start = all_if_same(termini.iter().filter_map(TrxPos::offset_from_cds_start));
        let vs_cds_end = all_if_same(termini.iter().filter_map(TrxPos::offset_from_cds_end));
        let frames: Vec<usize> = termini
            .iter()
            .filter_map(move |trxpos| body_frame(cdsbody, trxpos))
            .collect();

        if frames.len() > 1 {
            GeneFrameResult::Ambig
        } else {
            GeneFrameResult::Good(GeneFraming {
                gene: gene,
                vs_cds_start: vs_cds_start,
                vs_cds_end: vs_cds_end.map(|x| x + 3),
                frame: all_if_same(frames.into_iter()),
                fp_length: fp_length,
            })
        }
    }
}

/// Returns `Some` if the iterator yields one or more items, all of
/// which are equal, and `None` otherwise.
fn all_if_same<T: Eq, I: Iterator<Item = T>>(mut iter: I) -> Option<T> {
    iter.next().map_or(None, |x0| {
        if iter.all(|x| x == x0) {
            Some(x0)
        } else {
            None
        }
    })
}

/// Returns the reading frame position (0, 1, or 2) of the transcript
/// position, provided it is within the "body" of the CDS. If the
/// position is outside the body, or if the transcript has no CDS at
/// all, then `None` is returned.
///
/// # Arguments
///
/// * `cdsbody` describes the body of the CDS. The first element is
/// the minimum offset relative to the start codon and the second
/// element is the maximum offset relative to the stop codon (and is
/// usually negative). For example, `(+15, -18)` would require
///
/// `cds_start + 15 <= pos <= cds_end + (-18)`
///
/// * `trxpos` is the transcript position
pub fn body_frame<'a, S>(cdsbody: &(isize, isize), trxpos: &TrxPos<'a, S>) -> Option<usize> {
    let vs_start = trxpos.offset_from_cds_start()?;
    let vs_end = trxpos.offset_from_cds_end()?;
    if vs_start >= cdsbody.0 && vs_end <= cdsbody.1 {
        trxpos.cds_frame()
    } else {
        None
    }
}

/// Returns the transcript position of the 5' end of a footprint,
/// provided it is compatible with the transcript, or `None`
/// otherwise. Compatibility is determined by `splice_compatible()`,
/// which requires the footprint to lie on the same strand and form a
/// contiguous sub-region of the overall transcript.
///
/// # Arguments
/// * `fp` is the location of the footprint
/// * `trx` is the transcript annotation
pub fn fp_into_transcript<'a>(
    fp: &Spliced<Rc<String>, ReqStrand>,
    trx: &'a Transcript<Rc<String>>,
) -> Option<TrxPos<'a, Rc<String>>> {
    if splice_compatible(&trx.loc(), fp) {
        let pos = trx
            .loc()
            .pos_into(&fp.first_pos())
            .expect("pos_into(first_pos) failed after splice_compatible() = true");
        assert!(pos.strand() == ReqStrand::Forward);
        assert!(pos.pos() >= 0);
        Some(TrxPos::new(trx, pos.pos() as usize))
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use bio::io::bed;
    use bio_types::annot::contig::*;
    use bio_types::annot::pos::*;
    use bio_types::annot::refids::RefIDSet;
    use bio_types::annot::spliced::*;

    fn fp(fp_str: &str) -> Spliced<Rc<String>, ReqStrand> {
        fp_str.parse().unwrap()
    }

    fn pos(pos_str: &str) -> Pos<Rc<String>, ReqStrand> {
        pos_str.parse().unwrap()
    }

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

    fn transcriptome_from_str(bedstr: &str) -> Transcriptome<Rc<String>> {
        let mut refids = RefIDSet::new();
        Transcriptome::new_from_bed(bed::Reader::new(bedstr.as_bytes()).records(), &mut refids)
            .expect("Transcriptome from string")
    }

    #[test]
    fn test_fp_into_trx() {
        // [87261..87387) [87387..87500) [87500..87822)
        // [0    ..126)                  [126  ..448)
        let fwd_str = "chr01	87261	87822	YAL030W	0	+	87285	87752	0	2	126,322,	0,239,\n";
        let fwd_trx = transcript_from_str(&fwd_str);

        // [2906..4116) [4116..4215) [4215..5009)
        // (2004..794]               (794 ..0]
        let rev_str = "chr02	2906	5009	YBL111C	0	-	2906	5009	0	2	1210,794,	0,1309,\n";
        let rev_trx = transcript_from_str(&rev_str);

        fn into(
            fp: &Spliced<Rc<String>, ReqStrand>,
            trx: &Transcript<Rc<String>>,
        ) -> Option<(String, usize)> {
            fp_into_transcript(fp, trx)
                .map(|trxpos| (trxpos.transcript().trxname().to_string(), trxpos.pos()))
        }

        assert_eq!(into(&fp("chr01:87260-87290(+)"), &fwd_trx), None);
        assert_eq!(
            into(&fp("chr01:87261-87290(+)"), &fwd_trx),
            Some(("YAL030W".to_string(), 0))
        );
        assert_eq!(into(&fp("chr01:87261-87290(-)"), &fwd_trx), None);
        assert_eq!(into(&fp("chr02:87261-87290(+)"), &fwd_trx), None);
        assert_eq!(
            into(&fp("chr01:87361-87387(+)"), &fwd_trx),
            Some(("YAL030W".to_string(), 100))
        );
        assert_eq!(into(&fp("chr01:87361-87388(+)"), &fwd_trx), None);
        assert_eq!(
            into(&fp("chr01:87361-87387;87500-87501(+)"), &fwd_trx),
            Some(("YAL030W".to_string(), 100))
        );
        assert_eq!(
            into(&fp("chr01:87361-87387;87499-87501(+)"), &fwd_trx),
            None
        );
        assert_eq!(
            into(&fp("chr01:87361-87388;87500-87501(+)"), &fwd_trx),
            None
        );
        assert_eq!(
            into(&fp("chr01:87500-87530(+)"), &fwd_trx),
            Some(("YAL030W".to_string(), 126))
        );
        assert_eq!(
            into(&fp("chr01:87800-87822(+)"), &fwd_trx),
            Some(("YAL030W".to_string(), 426))
        );
        assert_eq!(into(&fp("chr01:87800-87823(+)"), &fwd_trx), None);

        assert_eq!(into(&fp("chr02:4980-5010(-)"), &rev_trx), None);
        assert_eq!(
            into(&fp("chr02:4980-5009(-)"), &rev_trx),
            Some(("YBL111C".to_string(), 0))
        );
        assert_eq!(into(&fp("chr02:4980-5009(+)"), &rev_trx), None);
        assert_eq!(into(&fp("chr01:4980-5009(-)"), &rev_trx), None);
        assert_eq!(
            into(&fp("chr02:4780-4809(-)"), &rev_trx),
            Some(("YBL111C".to_string(), 200))
        );
        assert_eq!(into(&fp("chr02:4780-4790;4795-4809(-)"), &rev_trx), None);
        assert_eq!(
            into(&fp("chr02:4215-4245(-)"), &rev_trx),
            Some(("YBL111C".to_string(), 764))
        );
        assert_eq!(into(&fp("chr02:4214-4245(-)"), &rev_trx), None);
        assert_eq!(
            into(&fp("chr02:4107-4116;4215-4235(-)"), &rev_trx),
            Some(("YBL111C".to_string(), 774))
        );
        assert_eq!(into(&fp("chr02:4107-4116;4214-4235(-)"), &rev_trx), None);
        assert_eq!(into(&fp("chr02:4107-4116;4216-4235(-)"), &rev_trx), None);
        assert_eq!(into(&fp("chr02:4107-4115;4215-4235(-)"), &rev_trx), None);
        assert_eq!(into(&fp("chr02:4107-4117;4215-4235(-)"), &rev_trx), None);
        assert_eq!(
            into(&fp("chr02:4107-4116;4170-4180;4215-4235(-)"), &rev_trx),
            None
        );
        assert_eq!(
            into(&fp("chr02:3985-4116(-)"), &rev_trx),
            Some(("YBL111C".to_string(), 794))
        );
        assert_eq!(
            into(&fp("chr02:2906-2936(-)"), &rev_trx),
            Some(("YBL111C".to_string(), 1974))
        );
        assert_eq!(into(&fp("chr02:2905-2936(-)"), &rev_trx), None);
    }

    #[test]
    fn test_body_frame() {
        // CDS is [83..1337)
        let fwd_str = "chr01	33364	34785	YAL061W	0	+	33447	34701	0	1	1421,	0,\n";
        let fwd_trx = transcript_from_str(&fwd_str);

        // CDS is [101..842)
        let rev_str = "chr01	51775	52696	YAL049C	0	-	51854	52595	0	1	921,	0,\n";
        let rev_trx = transcript_from_str(&rev_str);

        assert_eq!(body_frame(&(15, -15), &TrxPos::new(&fwd_trx, 63)), None);
        assert_eq!(body_frame(&(15, -15), &TrxPos::new(&fwd_trx, 83)), None);
        assert_eq!(body_frame(&(15, -15), &TrxPos::new(&fwd_trx, 97)), None);
        assert_eq!(body_frame(&(15, -15), &TrxPos::new(&fwd_trx, 98)), Some(0));
        assert_eq!(body_frame(&(15, -15), &TrxPos::new(&fwd_trx, 99)), Some(1));
        assert_eq!(body_frame(&(15, -15), &TrxPos::new(&fwd_trx, 100)), Some(2));
        assert_eq!(body_frame(&(15, -15), &TrxPos::new(&fwd_trx, 101)), Some(0));
        assert_eq!(body_frame(&(16, -15), &TrxPos::new(&fwd_trx, 101)), Some(0));
        assert_eq!(body_frame(&(17, -15), &TrxPos::new(&fwd_trx, 101)), Some(0));
        assert_eq!(body_frame(&(15, -15), &TrxPos::new(&fwd_trx, 200)), Some(0));
        assert_eq!(
            body_frame(&(15, -15), &TrxPos::new(&fwd_trx, 1321)),
            Some(2)
        );
        assert_eq!(
            body_frame(&(15, -15), &TrxPos::new(&fwd_trx, 1322)),
            Some(0)
        );
        assert_eq!(body_frame(&(15, -15), &TrxPos::new(&fwd_trx, 1323)), None);
        assert_eq!(body_frame(&(15, -15), &TrxPos::new(&fwd_trx, 1353)), None);

        assert_eq!(body_frame(&(12, -18), &TrxPos::new(&rev_trx, 50)), None);
        assert_eq!(body_frame(&(12, -18), &TrxPos::new(&rev_trx, 101)), None);
        assert_eq!(body_frame(&(12, -18), &TrxPos::new(&rev_trx, 112)), None);
        assert_eq!(body_frame(&(12, -18), &TrxPos::new(&rev_trx, 113)), Some(0));
        assert_eq!(body_frame(&(11, -18), &TrxPos::new(&rev_trx, 113)), Some(0));
        assert_eq!(body_frame(&(12, -18), &TrxPos::new(&rev_trx, 114)), Some(1));
        assert_eq!(body_frame(&(12, -18), &TrxPos::new(&rev_trx, 115)), Some(2));
        assert_eq!(body_frame(&(12, -18), &TrxPos::new(&rev_trx, 116)), Some(0));
        assert_eq!(body_frame(&(12, -18), &TrxPos::new(&rev_trx, 215)), Some(0));
        assert_eq!(body_frame(&(12, -18), &TrxPos::new(&rev_trx, 823)), Some(2));
        assert_eq!(body_frame(&(12, -18), &TrxPos::new(&rev_trx, 824)), Some(0));
        assert_eq!(body_frame(&(12, -18), &TrxPos::new(&rev_trx, 825)), None);
    }

    #[test]
    fn gene_framing_1trx() {
        // [87261..87387) [87387..87500) [87500..87822)
        // [0    ..126)                  [126  ..448)
        // CDS is 24..378
        let fwd_str = "chr01	87261	87822	YAL030W	0	+	87285	87752	0	2	126,322,	0,239,\n";
        let fwd_trx = transcript_from_str(&fwd_str);
        assert_eq!(fwd_trx.cds_range(), &Some(24..378));

        // CDS body is (15, -15)
        fn frame(fp_str: &str, trx: &Transcript<Rc<String>>) -> String {
            let fp: Spliced<Rc<String>, ReqStrand> = fp_str.parse().expect("Error parsing fp");
            let gfr = gene_framing(&(15, -15), &vec![trx], &fp);
            String::from_utf8(gfr.aux()).expect("Bad UTF8")
        }

        assert_eq!(frame("chr01:87276-87305(+)", &fwd_trx), "YAL030W/-9/-363/*");
        assert_eq!(
            frame("chr01:87297-87325(+)", &fwd_trx),
            "YAL030W/+12/-342/*"
        );
        assert_eq!(
            frame("chr01:87299-87327(+)", &fwd_trx),
            "YAL030W/+14/-340/*"
        );
        assert_eq!(
            frame("chr01:87300-87328(+)", &fwd_trx),
            "YAL030W/+15/-339/+0"
        );
        assert_eq!(
            frame("chr01:87301-87329(+)", &fwd_trx),
            "YAL030W/+16/-338/+1"
        );
        assert_eq!(
            frame("chr01:87302-87330(+)", &fwd_trx),
            "YAL030W/+17/-337/+2"
        );
        assert_eq!(
            frame("chr01:87303-87331(+)", &fwd_trx),
            "YAL030W/+18/-336/+0"
        );

        assert_eq!(frame("chr01:87260-87287(+)", &fwd_trx), "NoCompatible");
        assert_eq!(frame("chr01:87361-87388(+)", &fwd_trx), "NoCompatible");
        assert_eq!(
            frame("chr01:87300-87315;87345-87360(+)", &fwd_trx),
            "NoCompatible"
        );
        assert_eq!(frame("chr01:87499-87526(+)", &fwd_trx), "NoCompatible");
        assert_eq!(
            frame("chr01:87375-87387;87450-87453;87500-87514(+)", &fwd_trx),
            "NoCompatible"
        );
        assert_eq!(frame("chr01:87800-87823(+)", &fwd_trx), "NoCompatible");

        assert_eq!(
            frame("chr01:87375-87387;87500-87514(+)", &fwd_trx),
            "YAL030W/+90/-264/+0"
        );
        assert_eq!(
            frame("chr01:87376-87387;87500-87516(+)", &fwd_trx),
            "YAL030W/+91/-263/+1"
        );

        assert_eq!(
            frame("chr01:87722-87750(+)", &fwd_trx),
            "YAL030W/+324/-30/+0"
        );
        assert_eq!(
            frame("chr01:87737-87765(+)", &fwd_trx),
            "YAL030W/+339/-15/+0"
        );
        assert_eq!(
            frame("chr01:87738-87765(+)", &fwd_trx),
            "YAL030W/+340/-14/*"
        );
        assert_eq!(frame("chr01:87756-87784(+)", &fwd_trx), "YAL030W/+358/+4/*");

        validate_framing(&fwd_trx, 28, (15, -15));

        // CDS is [101..842)
        let rev_str = "chr01	51775	52696	YAL049C	0	-	51854	52595	0	1	921,	0,\n";
        let rev_trx = transcript_from_str(&rev_str);

        validate_framing(&rev_trx, 28, (15, -15));
    }

    fn validate_framing(trx: &Transcript<Rc<String>>, fplen: isize, cdsbody: (isize, isize)) {
        for i in 0..(trx.loc().exon_total_length() as isize - fplen) {
            let trx_first = Pos::new(trx.trxname().clone(), i, ReqStrand::Forward);
            let chr_first = trx
                .loc()
                .pos_outof(&trx_first)
                .expect(&format!("Cannot pull fp start {} @ {} out", trx_first, i));
            let trx_last = Pos::new(trx.trxname().clone(), i + fplen - 1, ReqStrand::Forward);
            let chr_last = trx
                .loc()
                .pos_outof(&trx_last)
                .expect(&format!("Cannot pull fp last {} @ {} out", trx_last, i));
            let chr_length = (1 + chr_last.pos() - chr_first.pos()).abs() as usize;
            let chr_span = Contig::with_first_length(&chr_first, chr_length)
                .expect("Cannot make fp chr contig");
            let chr_fp = trx
                .loc()
                .contig_intersection(&chr_span)
                .expect("Cannot intersect fp chr contig");
            let trxs = vec![trx];
            let gf = match gene_framing(&cdsbody, &trxs, &chr_fp) {
                GeneFrameResult::Good(gf) => gf,
                _ => panic!("No gene framing"),
            };

            assert_eq!(&gf.gene, trx.gene_ref());
            assert_eq!(
                gf.vs_cds_start.as_ref().map(|vs_start| i - vs_start),
                trx.cds_range().as_ref().map(|cds| cds.start as isize)
            );
            assert_eq!(
                gf.vs_cds_end.as_ref().map(|vs_end| i - vs_end),
                trx.cds_range().as_ref().map(|cds| cds.end as isize)
            );
            if gf
                .vs_cds_start
                .as_ref()
                .map_or(true, |vs_start| *vs_start < cdsbody.0)
            {
                assert_eq!(gf.frame, None);
            } else if gf
                .vs_cds_end
                .as_ref()
                .map_or(true, |vs_end| *vs_end > cdsbody.1)
            {
                assert_eq!(gf.frame, None);
            } else {
                assert!(gf.frame.as_ref().map_or(false, |fr| {
                    gf.vs_cds_start
                        .as_ref()
                        .map_or(false, |vs_start| (*vs_start - *fr as isize) % 3 == 0)
                }));
            }
        }
    }
}
