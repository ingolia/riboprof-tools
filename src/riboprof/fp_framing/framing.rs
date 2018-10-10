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

/// Framing information for a footprint, relative to a
/// gene. This includes the offsets relative to the CDS start and end,
/// as well as the reading frame position.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GeneFraming {
    vs_cds_start: Option<isize>,
    vs_cds_end: Option<isize>,
    frame: Option<usize>,
}

/// Computes framing information for a footprint, relative to a gene
/// defined by one or more transcripts. Framing in
pub fn gene_framing<'a>(
    cdsbody: &(isize, isize),
    trxs: &[&'a Transcript<Rc<String>>],
    fp: &Spliced<Rc<String>, ReqStrand>,
) -> GeneFrameResult {
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
                vs_cds_start: vs_cds_start,
                vs_cds_end: vs_cds_end,
                frame: all_if_same(frames.into_iter()),
            })
        }
    }
}

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

}
