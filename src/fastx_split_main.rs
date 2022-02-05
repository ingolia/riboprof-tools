extern crate failure;
#[macro_use]
extern crate clap;

extern crate riboprof;

use std::io;
use std::io::Write;
use std::process;

use clap::{App, Arg};

use riboprof::fastx_split::*;

fn main() {
    match wrapper() {
        Err(e) => {
            io::stderr().write(format!("{}\n", e).as_bytes()).unwrap();
            process::exit(1);
        }
        _ => (),
    };
}

fn wrapper() -> Result<(), failure::Error> {
    let cli = get_cli()?;
    let config = Config::new(&cli)?;
    fastx_split(config)
}

fn get_cli() -> Result<CLI, failure::Error> {
    let matches = App::new("fastx-split")
        .version("0.1.0")
        .author("Nick Ingolia <ingolia@berkeley.edu>")
        .about("Split FastQ file using index and random nucleotides")
        .arg(
            Arg::new("output_dir")
                .short('o')
                .long("output-dir")
                .value_name("OUTPUT-DIR")
                .help("Output directory name")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::new("min_insert")
                .short('m')
                .long("min-insert")
                .value_name("MIN-INSERT")
                .help("Minimum insert length")
                .takes_value(true)
                .default_value("0"),
        )
        .arg(
            Arg::new("prefix")
                .short('p')
                .long("prefix")
                .value_name("PREFIX")
                .help("Prefix format string")
                .takes_value(true)
                .default_value(""),
        )
        .arg(
            Arg::new("suffix")
                .short('x')
                .long("suffix")
                .value_name("SUFFIX")
                .help("Suffix format string")
                .takes_value(true)
                .default_value(""),
        )
        .arg(
            Arg::new("sample_sheet")
                .short('s')
                .long("sample-sheet")
                .value_name("SAMPLESHEET.CSV")
                .help("File name of CSV-format sample sheet")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::new("progress")
                .long("progress")
                .value_name("NSEQS")
                .help("Report progress every NSEQS sequences")
                .takes_value(true)
                .default_value("0"),
        )
        .arg(Arg::new("input").multiple_occurrences(true).required(true).allow_invalid_utf8(true))
        .get_matches();

    Ok(CLI {
        fastx_inputs: matches.values_of_lossy("input").unwrap(),
        output_dir: matches.value_of("output_dir").unwrap().to_string(),
        min_insert: value_t!(matches.value_of("min_insert"), usize)?,
        prefix: matches.value_of("prefix").unwrap().to_string(),
        suffix: matches.value_of("suffix").unwrap().to_string(),
        sample_sheet: matches.value_of("sample_sheet").unwrap().to_string(),
        progress: value_t!(matches.value_of("progress"), usize)?,
    })
}
