extern crate failure;
#[macro_use]
extern crate clap;

extern crate bio;

use std::fs;
use std::io::Write;
use std::path::PathBuf;

use clap::{Arg, App};

#[derive(Debug)]
struct Config {
    fastx_inputs: Vec<PathBuf>,
    out_dir: PathBuf,
    min_insert: usize,
    // linker_format: LinkerFormat,
    sample_sheet: PathBuf,
    progress: Option<usize>,
}

fn main() {
    let matches = App::new("fastx-split")
        .version("0.1.0")
        .author("Nick Ingolia <ingolia@berkeley.edu>")
        .about("Split FastQ file using index and random nucleotides")
        .arg(Arg::with_name("output_dir")
             .short("o")
             .long("output-dir")
             .value_name("OUTPUT-DIR")
             .help("Output directory name")
             .takes_value(true)
             .required(true))
        .arg(Arg::with_name("min_insert")
             .short("m")
             .long("min-insert")
             .value_name("MIN-INSERT")
             .help("Minimum insert length")
             .takes_value(true)
             .default_value("0"))
        .arg(Arg::with_name("prefix")
             .short("p")
             .long("prefix")
             .value_name("PREFIX")
             .help("Prefix format string")
             .takes_value(true)
             .default_value(""))
        .arg(Arg::with_name("suffix")
             .short("x")
             .long("suffix")
             .value_name("SUFFIX")
             .help("Suffix format string")
             .takes_value(true)
             .default_value(""))
        .arg(Arg::with_name("sample_sheet")
             .short("s")
             .long("sample-sheet")
             .value_name("SAMPLESHEET.CSV")
             .help("File name of CSV-format sample sheet")
             .takes_value(true)
             .required(true))
        .arg(Arg::with_name("progress")
             .long("progress")
             .value_name("NSEQS")
             .help("Report progress every NSEQS sequences")
             .takes_value(true))
        .arg(Arg::with_name("input")
             .multiple(true)
             .required(true))
        .get_matches();
}
