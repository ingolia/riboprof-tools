use std::fs;
use std::io::Write;
use std::path::Path;
use std::rc::Rc;

use anyhow::Result;
use bio::io::bed;
use bio_types::annot::refids::RefIDSet;
use rust_htslib::bam;
use rust_htslib::bam::Read as BamRead;

use crate::bam_utils::Tids;
use crate::codon_assign::ASites;
use crate::parse_pair;
use crate::transcript::Transcriptome;

mod count;

use crate::fp_count::count::{Count, CountConfig};
// mod framing;
// mod stats;

// use crate::fp_framing::framing::*;
// use crate::fp_framing::stats::*;

pub struct CLI {
    pub input: String,
    pub output: String,
    pub bed: String,
    pub asite: String,
    pub cds_insets: String,
    // pub stats: Option<String>,
    // pub framing: Option<String>,
    pub nhits: bool,
    pub whole: bool,
    pub debug: bool,
    pub reverse: bool,
}

pub struct Config {
    input: bam::IndexedReader,
    output: Box<dyn std::io::Write>,
    trxome: Transcriptome<Rc<String>>,
    refids: RefIDSet<Rc<String>>,
    count_config: CountConfig,
    // asites: ASites,
    // cds_insets: Range<isize>,
    // nhits: bool,
    // whole: bool,
    #[allow(dead_code)]
    debug: bool,
    // strand: ReqStrand,
}

impl Config {
    pub fn new(cli: &CLI) -> Result<Self> {
        if cli.input == "-" {
            eprintln!(
                "Warning: indexed BAM file required, not interpreting {:?} as stdin",
                cli.input
            )
        }

        let input = bam::IndexedReader::from_path(Path::new(&cli.input))?;
        let mut refids = RefIDSet::new();

        let trxome =
            Transcriptome::new_from_bed(bed::Reader::from_file(&cli.bed)?.records(), &mut refids)?;
        let asites = ASites::new_from_file(&cli.asite)?;

        let output = fs::File::create(&cli.output)?;

        Ok(Config {
            input,
            output: Box::new(output),
            trxome,
            refids,
            count_config: CountConfig::new(
                asites,
                parse_pair(&cli.cds_insets)?,
                cli.nhits,
                cli.whole,
                cli.reverse,
            ),
            debug: cli.debug,
        })
    }
}

pub fn run_fp_count_from_cli(cli: &CLI) -> Result<()> {
    run_fp_count(Config::new(cli)?)
}

pub fn run_fp_count(mut config: Config) -> Result<()> {
    let tids = Tids::new(&mut config.refids, &config.input.header());

    for trx in config.trxome.transcripts() {
        let mut trx_count = Count::new();
        config.input.fetch(trx.fetch_desc())?;

        for recres in config.input.records() {
            let rec = recres?;
            trx_count.tally_record(&config.count_config, &trx, &tids, &rec)?;
        }

        write!(
            config.output,
            "{}\t{}\t{:0.1}\n",
            trx.trxname(),
            trx.loc().exon_total_length(),
            trx_count.compat()
        )?;
    }

    Ok(())
}
