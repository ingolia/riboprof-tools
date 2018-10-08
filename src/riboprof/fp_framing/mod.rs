use std::error::Error;
use std::fmt;
use std::fs;
use std::io::{self, Read, Write};
use std::ops::Range;
use std::path::{Path, PathBuf};
use std::rc::Rc;
use std::str;

use failure;

use bio::io::bed;
use bio_types::annot::*;
use rust_htslib::bam;
use rust_htslib::bam::Read as BamRead;

use transcript::*;

mod stats;

use fp_framing::stats::*;

pub struct CLI {
    pub input: String,
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
    input: String,
    output: PathBuf,
    trxome: Transcriptome<Rc<String>>,
    flanking: Range<isize>,
    cdsbody: Range<isize>,
    lengths: Range<usize>,
    count_multi: bool,
    annotate: Option<PathBuf>,
}

impl Config {
    pub fn new(cli: &CLI) -> Result<Self, failure::Error> {
        let trxome = Self::read_transcriptome(&cli)?;

        Ok(Config {
            input: cli.input.to_string(),
            output: Path::new(&cli.output).to_path_buf(),
            trxome: trxome,
            flanking: Self::parse_pair(&cli.flanking)?,
            cdsbody: Self::parse_pair(&cli.cdsbody)?,
            lengths: Self::parse_pair(&cli.lengths)?,
            count_multi: cli.count_multi,
            annotate: cli.annotate.as_ref().map(|ann| Path::new(&ann).to_path_buf()),
        })
    }

    fn read_transcriptome(cli: &CLI) -> Result<Transcriptome<Rc<String>>, failure::Error> {
        // ZZZ Handle Trx->Gene mappings
        let mut refids = refids::RefIDSet::new();
        let mut trxome = Transcriptome::new();

        for recres in bed::Reader::from_file(&cli.bed)?.records() {
            let rec = recres?;
            let trx = Transcript::from_bed12(&rec, &mut refids)?;
            trxome.insert(trx)?;
        }

        Ok(trxome)
    }

    fn parse_pair<I>(pair_str: &str) -> Result<Range<I>, failure::Error>
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
            Err(failure::err_msg(format!(
                "Expecting integer pair \"a,b\" but got \"{}\"",
                pair_str
            )))
        }
    }
}

pub fn fp_framing(config: Config) -> Result<(), failure::Error> {
    let mut input = if config.input == "-" {
        bam::Reader::from_stdin()?
    } else {
        bam::Reader::from_path(Path::new(&config.input))?
    };
    
    let mut annotate = match config.annotate {
        None => None,
        Some(ref annot_file) => {
            let header = bam::Header::from_template(input.header());
            Some(bam::Writer::from_path(Path::new(&annot_file), &header)?)
        }
    };

    let mut align_stats = AlignStats::new();

    for recres in input.records() {
        let mut rec = recres?;

        let tag = record_framing(&config, &mut align_stats, &rec);

        if let Some(ref mut ann_writer) = &mut annotate {
            rec.push_aux(b"ZF", &bam::record::Aux::String(tag.as_bytes()))?; 
            ann_writer.write(&rec)?;
        }
    }

    Ok(())
}

pub fn record_framing(config: &Config, align_stats: &mut AlignStats, rec: &bam::Record) -> String {
    "N/A".to_string()
}

#[derive(Debug)]
pub enum FpFramingError {
    BadArgument(String),
}

impl Error for FpFramingError {}

impl fmt::Display for FpFramingError {
    fn fmt(&self, f: &mut fmt::Formatter) -> Result<(), fmt::Error> {
        match self {
            FpFramingError::BadArgument(msg) => write!(f, "Bad argument: {}", msg),
        }
    }
}
