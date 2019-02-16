use std::error::Error;
use std::fmt;
use std::fs;
use std::path::{Path, PathBuf};

use failure;

use rust_htslib::bam;
use rust_htslib::bam::Read as BamRead;

use bam_utils::*;

pub struct CLI {
    pub bam_input: String,
    pub bam_output: String,
    pub bam_dups: Option<String>,
    pub stats: Option<String>,
    pub annotate: bool,
}

pub struct Config {
    bam_input: PathBuf,
    bam_output: PathBuf,
    bam_dups: Option<PathBuf>,
    stats: Option<PathBuf>,
    annotate: bool,
}

impl Config {
    pub fn new(cli: &CLI) -> Result<Self, failure::Error> {
        Err(failure::err_msg(format!("Unimplemented")))
    }
}

pub fn bam_suppress_duplicates(mut config: Config) -> Result<(), failure::Error> {
    Ok(())
}
