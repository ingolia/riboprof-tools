use std::ops::Range;
use std::rc::Rc;

use bio_types::annot::loc::Loc;
use bio_types::annot::pos::Pos;
use bio_types::annot::spliced::Spliced;
use bio_types::strand::*;
use rust_htslib::bam;

use fp_framing::stats::*;

use codon_assign::*;
use transcript::*;

pub fn record_framing(
    trxome: &Transcriptome<Rc<String>>,
    rec: &bam::Record,
    framing_stats: &mut FramingStats,
    cdsbody: &(isize, isize),
    count_multi: bool,
) -> String {
    "N/A".to_string()
}

pub fn footprint_framing(
    trxome: &Transcriptome<Rc<String>>,
    rec: &bam::Record,
    framing_stats: &mut FramingStats,
    cdsbody: &(isize, isize),
    count_multi: bool,
) -> FpFrameResult {
    unimplemented!()
}

pub enum FpFrameResult {
    Gene(GeneFrameResult),
    NoncodingOnly,
    NoncodingOverlap,
    MultiCoding
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
}

impl GeneFraming {
    pub fn aux(&self) -> Vec<u8> {
        let mut buf = self.gene.to_string();
        buf += &Self::opt_to_str(&self.vs_cds_start);
        buf += &Self::opt_to_str(&self.vs_cds_end);
        buf += &Self::opt_to_str(&self.frame.as_ref().map(|x| *x as isize));
        buf.into_bytes()
    }

    fn opt_to_str(opt: &Option<isize>) -> String 
    {
        opt.map_or_else(|| "/*".to_string(), |x| format!("/{:+}", x))
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
                vs_cds_end: vs_cds_end,
                frame: all_if_same(frames.into_iter()),
            })
        }
    }
}

/// Returns `Some` if the iterator yields one or more items, all of
/// which are equal, and `None` otherwise.
pub fn all_if_same<T: Eq, I: Iterator<Item = T>>(mut iter: I) -> Option<T> {
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
pub fn body_frame<'a>(cdsbody: &(isize, isize), trxpos: &TrxPos<'a, Rc<String>>) -> Option<usize> {
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

        fn into<'a>(fp: &Spliced<Rc<String>, ReqStrand>, trx: &Transcript<Rc<String>>) -> Option<(String, usize)> {
            fp_into_transcript(fp, trx).map(|trxpos| (trxpos.transcript().trxname().to_string(), trxpos.pos()))
        }

        assert_eq!(into(&fp("chr01:87260-87290(+)"), &fwd_trx), None);
        assert_eq!(into(&fp("chr01:87261-87290(+)"), &fwd_trx), Some(("YAL030W".to_string(), 0)));
        assert_eq!(into(&fp("chr01:87261-87290(-)"), &fwd_trx), None);
        assert_eq!(into(&fp("chr02:87261-87290(+)"), &fwd_trx), None);
        assert_eq!(into(&fp("chr01:87361-87387(+)"), &fwd_trx), Some(("YAL030W".to_string(), 100)));
        assert_eq!(into(&fp("chr01:87361-87388(+)"), &fwd_trx), None);
        assert_eq!(into(&fp("chr01:87361-87387;87500-87501(+)"), &fwd_trx), Some(("YAL030W".to_string(), 100)));
        assert_eq!(into(&fp("chr01:87361-87387;87499-87501(+)"), &fwd_trx), None);
        assert_eq!(into(&fp("chr01:87361-87388;87500-87501(+)"), &fwd_trx), None);
        assert_eq!(into(&fp("chr01:87500-87530(+)"), &fwd_trx), Some(("YAL030W".to_string(), 126)));
        assert_eq!(into(&fp("chr01:87800-87822(+)"), &fwd_trx), Some(("YAL030W".to_string(), 426)));
        assert_eq!(into(&fp("chr01:87800-87823(+)"), &fwd_trx), None);

        assert_eq!(into(&fp("chr02:4980-5010(-)"), &rev_trx), None);
        assert_eq!(into(&fp("chr02:4980-5009(-)"), &rev_trx), Some(("YBL111C".to_string(), 0)));
        assert_eq!(into(&fp("chr02:4980-5009(+)"), &rev_trx), None);
        assert_eq!(into(&fp("chr01:4980-5009(-)"), &rev_trx), None);
        assert_eq!(into(&fp("chr02:4780-4809(-)"), &rev_trx), Some(("YBL111C".to_string(), 200)));
        assert_eq!(into(&fp("chr02:4780-4790;4795-4809(-)"), &rev_trx), None);
        assert_eq!(into(&fp("chr02:4215-4245(-)"), &rev_trx), Some(("YBL111C".to_string(), 764)));
        assert_eq!(into(&fp("chr02:4214-4245(-)"), &rev_trx), None);
        assert_eq!(into(&fp("chr02:4107-4116;4215-4235(-)"), &rev_trx), Some(("YBL111C".to_string(), 774)));
        assert_eq!(into(&fp("chr02:4107-4116;4214-4235(-)"), &rev_trx), None);
        assert_eq!(into(&fp("chr02:4107-4116;4216-4235(-)"), &rev_trx), None);
        assert_eq!(into(&fp("chr02:4107-4115;4215-4235(-)"), &rev_trx), None);
        assert_eq!(into(&fp("chr02:4107-4117;4215-4235(-)"), &rev_trx), None);
        assert_eq!(into(&fp("chr02:4107-4116;4170-4180;4215-4235(-)"), &rev_trx), None);
        assert_eq!(into(&fp("chr02:3985-4116(-)"), &rev_trx), Some(("YBL111C".to_string(), 794)));
        assert_eq!(into(&fp("chr02:2906-2936(-)"), &rev_trx), Some(("YBL111C".to_string(), 1974)));
        assert_eq!(into(&fp("chr02:2905-2936(-)"), &rev_trx), None);
    }
}
