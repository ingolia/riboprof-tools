extern crate failure;
#[macro_use]
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
    fastx_split(config)
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
            Arg::with_name("annotate")
                .short("a")
                .long("annotate")
                .value_name("ANNOTATED.BAM")
                .help("Write output BAM file annotated wiht framing information")
                .takes_value(true)
        )
        
        .arg(Arg::with_name("input").value_name("INPUT.BAM").required(true))
        .get_matches();

    /// ZZZ
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
