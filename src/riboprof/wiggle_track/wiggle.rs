use std::io::Write;

use anyhow::Result;
use bio_types::annot::pos::Pos;
use bio_types::strand::ReqStrand;
use rust_htslib::bam;

use crate::bam_utils::bam_to_spliced_tid;
use crate::codon_assign::ASites;
use crate::wiggle_track::genome_count::{GenomeCounts, WindowCounts};

pub fn fp_asite(asites: &ASites, rec: &bam::Record) -> Result<Option<Pos<u32, ReqStrand>>> {
    Ok(bam_to_spliced_tid(&rec)?.and_then(|fp| asites.a_site(fp)))
}

pub fn write_wigs<W: Write>(
    gcount: &GenomeCounts<usize>,
    outputs: (&mut W, &mut W),
    qnorm: Option<f64>,
) -> Result<()> {
    let (output_fwd, output_rev) = outputs;

    for chr in gcount.chroms() {
        write!(output_fwd, "variableStep chrom={} span=1\n", chr.name())?;
        write_wig_chrom_strand(chr.fwd(), output_fwd, qnorm)?;

        write!(output_rev, "variableStep chrom={} span=1\n", chr.name())?;
        write_wig_chrom_strand(chr.rev(), output_rev, qnorm)?;
    }

    Ok(())
}

pub fn write_wig_chrom_strand<W: Write>(
    counts: &WindowCounts<usize>,
    output: &mut W,
    qnorm: Option<f64>,
) -> Result<()> {
    for (pos, ct) in counts.sparse_pos_iter().filter(|(_, ct)| **ct > 0) {
        if let Some(scale) = qnorm {
            write!(output, "{} {:0.3}\n", pos + 1, (*ct as f64) * scale)?;
        } else {
            write!(output, "{} {}\n", pos + 1, ct)?;
        }
    }
    Ok(())
}
