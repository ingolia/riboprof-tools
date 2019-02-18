extern crate clap;
extern crate failure;

extern crate riboprof;

use std::io;
use std::io::Write;
use std::process;

use clap::{App, Arg};

use riboprof::bam_suppress_duplicates::*;

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
    bam_suppress_duplicates(config)
}

fn get_cli() -> Result<CLI, failure::Error> {
    let matches = App::new("bam-suppress-duplicates")
        .version("0.1.0")
        .author("Nick Ingolia <ingolia@berkeley.edu>")
        .about("Suppress likely PCR duplicates based on UMIs embedded in sequence names")
        .arg(
            Arg::with_name("bam_input")
                .short("i")
                .long("input")
                .value_name("INPUT.BAM")
                .help("BAM format input file")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("bam_output")
                .short("o")
                .long("output")
                .value_name("OUTPUT.BAM")
                .help("BAM format output file")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("bam_dups")
                .short("d")
                .long("dups")
                .long("duplicates")
                .value_name("DUPLICATES.BAM")
                .help("BAM format file of duplicates")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("stats")
                .short("s")
                .long("stats")
                .long("statistics")
                .value_name("STATS.TXT")
                .help("Output file with duplicate statistics")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("annotate")
                .short("a")
                .long("annotate")
                .help("Annotate deduplicated reads"),
        )
        .get_matches();

    Ok(CLI {
        bam_input: matches.value_of("bam_input").unwrap().to_string(),
        bam_output: matches.value_of("bam_output").unwrap().to_string(),
        bam_dups: matches.value_of_lossy("bam_dups").map(|a| a.to_string()),
        stats: matches.value_of_lossy("stats").map(|a| a.to_string()),
        annotate: matches.is_present("annotate"),
    })
}
