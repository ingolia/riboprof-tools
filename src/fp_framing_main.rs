extern crate failure;
extern crate clap;

extern crate riboprof;

use std::io;
use std::io::Write;
use std::process;

use clap::{App, Arg};

use riboprof::fp_framing::*;

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
    run_fp_framing(config)
}

fn get_cli() -> Result<CLI, failure::Error> {
    let matches = App::new("fp-framing")
        .version("0.1.0")
        .author("Nick Ingolia <ingolia@berkeley.edu>")
        .about("Calculates ribosome profiling QC information including reading frame bias and start and stop codon meta-genes")
        .arg(
            Arg::with_name("output")
                .short("o")
                .long("output")
                .value_name("OUTBASE")
                .help("Base filename for output files")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("bed")
                .short("b")
                .long("bed")
                .value_name("BED")
                .help("BED-format annotation filename")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("genes")
                .short("g")
                .long("genes")
                .value_name("GENES.TXT")
                .help("Tab-delimited table of Transcript<TAB>Gene (or just Transcript to suppress a transcript)")
                .takes_value(true)
                .multiple(true)
                .number_of_values(1),
        )
        .arg(
            Arg::with_name("flanking")
                .short("f")
                .long("flanking")
                .value_name("START,END")
                .help("Range of profiles surrounding the start and end codons")
                .takes_value(true)
                .default_value("-100,100"),
        )
        .arg(
            Arg::with_name("cdsbody")
                .short("c")
                .long("cdsbody")
                .value_name("AFTERSTART,BEFOREEND")
                .help("Offsets from the start and end of the gene for framing analysis")
                .takes_value(true)
                .default_value("34,31"),
        )
        .arg(
            Arg::with_name("lengths")
                .short("l")
                .long("lengths")
                .value_name("MINLEN,MAXLEN")
                .help("Length frange for framing analysis")
                .takes_value(true)
                .default_value("26,34"),
        )
        .arg(
            Arg::with_name("count-multi")
                .short("m")
                .long("count-multi")
                .help("Count multi-mapping reads once, at their first occurrence (i.e., HI = 0)")
        )
        .arg(
            Arg::with_name("annotate")
                .short("a")
                .long("annotate")
                .value_name("ANNOTATED.BAM")
                .help("Write output BAM file annotated wiht framing information")
                .takes_value(true)
        )
        .arg(Arg::with_name("input").value_name("INPUT.BAM").required(true))
        .get_matches();

    Ok(CLI {
        output: matches.value_of("output").unwrap().to_string(),
        bed: matches.value_of("bed").unwrap().to_string(),
        genes: matches
            .values_of_lossy("genes")
            .unwrap_or_else(|| Vec::new()),
        flanking: matches.value_of("flanking").unwrap().to_string(),
        cdsbody: matches.value_of("cdsbody").unwrap().to_string(),
        lengths: matches.value_of("lengths").unwrap().to_string(),
        count_multi: matches.is_present("count-multi"),
        annotate: matches.value_of_lossy("annotate").map(|a| a.to_string()),
        input: matches.value_of("input").unwrap().to_string(),
    })
}
