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
        .long_about("Identifies and removes likely PCR duplicates. Duplicates are identified based on a nucleotide tag at the end of the read name stored in the BAM file, separated from the rest of the read name by a \"#\". Reads with no tag are not subject to deduplication. When multiple reads aligning to the same position share the same nucleotide tag, one is selected arbitrarily and written as the \"unique\" representative and, if specified, the rest are written to the file of duplicates. Optionally, the unique representative can be tagged with a \"ZD\" tag indicating the total number of duplicate reads (always 2 or more) at that position. Optionally, a table of duplicate suppression statistics can be written as a tab-separated file, tabulating the duplicate status of each distinct mapping site. In this statistics file, the first column is the total number of reads aligned to the site, the second is the number of unique reads, and the third is the count of distinct sites. Thus, \"1  1  234\" would indicate 234 distinct positions with a single unique read, \"2  2  17\" would indicate 17 distinct positions with two unique reads, and \"2  1  5\" would indicate 5 positions with a single duplicated read (2 reads total, 1 unique).")
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
