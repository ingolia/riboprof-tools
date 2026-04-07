use std::error::Error;
use std::ops::Range;
use std::str;

use anyhow::{Result, bail};

pub mod bam_suppress_duplicates;
pub mod bam_utils;
pub mod codon_assign;
pub mod fastx_split;
pub mod fp_count;
pub mod fp_framing;
pub mod metagene;
pub mod transcript;

fn parse_pair<I>(pair_str: &str) -> Result<Range<I>>
where
    I: str::FromStr,
    I::Err: Error + Send + Sized + Sync + 'static,
{
    let strs: Vec<&str> = pair_str.split(",").collect();
    if strs.len() == 2 {
        Ok(Range {
            start: strs[0].parse()?,
            end: strs[1].parse()?,
        })
    } else {
        bail!("Expecting integer pair \"a,b\" but got \"{}\"", pair_str)
    }
}
