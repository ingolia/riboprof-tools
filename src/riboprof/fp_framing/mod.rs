use std::fs;
use std::io::{self, Read, Write};
use std::path::{Path, PathBuf};
use std::str;

use failure;

use bio::io::fastq;

pub struct CLI {
    pub output: String,
    pub bed: String,
    pub genes: Vec<String>,
    pub flanking: String,
    pub cdsbody: String,
    pub lengths: String,
    pub count_multi: bool,
    pub annotate: Option<String>,
}

pub struct Config {

}

impl Config {
    pub fn new(cli: &CLI) -> Result<Self, failure::Error> {
        Ok(Config{})
    }
}

pub fn fp_framing(mut config: Config) -> Result<(), failure::Error> {
    Ok( () )
}
