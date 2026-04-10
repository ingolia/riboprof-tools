use std::io;
use std::io::Write;
use std::process;

use clap::{Arg, ArgAction, Command};

use anyhow::Result;
use riboprof::fp_framing::{CLI, run_fp_framing_cli};

fn main() {
    match wrapper() {
        Err(e) => {
            io::stderr().write(format!("{:#}\n", e).as_bytes()).unwrap();
            process::exit(1);
        }
        _ => (),
    };
}

fn wrapper() -> Result<()> {
    run_fp_framing_cli(get_cli()?)
}

fn get_cli() -> Result<CLI> {
    let matches = Command::new("fp-framing")
        .version("0.1.0")
        .author("Nick Ingolia <ingolia@berkeley.edu>")
        .about("Calculates ribosome profiling QC information including reading frame bias and start and stop codon meta-genes")
        .arg(
            Arg::new("output")
                .short('o')
                .long("output")
                .value_name("OUTBASE")
                .help("Base filename for output files")
                .action(ArgAction::Set)
                .required(true),
        )
        .arg(
            Arg::new("bed")
                .short('b')
                .long("bed")
                .value_name("BED")
                .help("BED-format annotation filename")
                .action(ArgAction::Set)
                .required(true),
        )
        .arg(
            Arg::new("genes")
                .short('g')
                .long("genes")
                .value_name("GENES.TXT")
                .help("Tab-delimited table of Transcript<TAB>Gene (or just Transcript to suppress a transcript)")
                .action(ArgAction::Append)
        )
        .arg(
            Arg::new("flanking")
                .short('f')
                .long("flanking")
                .value_name("START,END")
                .help("Range of profiles surrounding the start and end codons")
                .action(ArgAction::Set)
                .default_value("-100,100"),
        )
        .arg(
            Arg::new("cdsbody")
                .short('c')
                .long("cdsbody")
                .value_name("AFTERSTART,BEFOREEND")
                .help("Offsets from the start and end of the gene for framing analysis")
                .action(ArgAction::Set)
                .default_value("34,31"),
        )
        .arg(
            Arg::new("lengths")
                .short('l')
                .long("lengths")
                .value_name("MINLEN,MAXLEN")
                .help("Length frange for framing analysis")
                .action(ArgAction::Set)
                .default_value("26,34"),
        )
        .arg(
            Arg::new("count-multi")
                .short('m')
                .long("count-multi")
                .help("Count multi-mapping reads once, at their first occurrence (i.e., HI = 0)")
                .action(ArgAction::SetTrue)
        )
        .arg(
            Arg::new("annotate")
                .short('a')
                .long("annotate")
                .value_name("ANNOTATED.BAM")
                .help("Write output BAM file annotated wiht framing information")
                .action(ArgAction::Set)
        )
        .arg(Arg::new("input").value_name("INPUT.BAM").required(true))
        .get_matches();

    Ok(CLI {
        output: matches
            .get_one::<String>("output")
            .expect("output is missing")
            .to_string(),
        bed: matches
            .get_one::<String>("bed")
            .expect("bed file is missing")
            .to_string(),
        genes: matches
            .get_many::<String>("genes")
            .map_or_else(|| Vec::new(), |v| v.map(|s| s.to_string()).collect()),
        flanking: matches
            .get_one::<String>("flanking")
            .expect("flanking is missing")
            .to_string(),
        cdsbody: matches
            .get_one::<String>("cdsbody")
            .expect("cdsbody is missing")
            .to_string(),
        lengths: matches
            .get_one::<String>("lengths")
            .expect("length range is missing")
            .to_string(),
        count_multi: matches.get_flag("count-multi"),
        annotate: matches
            .get_many::<String>("annotate")
            .map(|r| r.map(|a| a.to_string()).collect()),
        input: matches
            .get_one::<String>("input")
            .expect("input is missing")
            .to_string(),
    })
}
