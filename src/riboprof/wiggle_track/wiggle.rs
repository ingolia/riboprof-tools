use anyhow::{Result, anyhow, bail};
use bio_types::annot::loc::Loc;
use bio_types::annot::pos::Pos;
use bio_types::annot::spliced::Spliced;
use bio_types::strand::ReqStrand;
use rust_htslib::bam;

use crate::bam_utils::{Tids, bam_to_spliced_tid};
use crate::codon_assign::ASites;
use crate::transcript::{Transcript, splice_compatible};

pub fn fp_asite(asites: &ASites, rec: &bam::Record) -> Result<Option<Pos<u32, ReqStrand>>> {
    Ok(bam_to_spliced_tid(&rec)?.and_then(|fp| asites.a_site(fp)))
}
