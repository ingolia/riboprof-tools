use std::io;
use std::io::Write;
use std::process;

use clap::{Arg, ArgAction, Command};

use anyhow::Result;
use riboprof::fp_count::{CLI, run_fp_count_from_cli};

fn main() {
    match wrapper() {
        Err(e) => {
            io::stderr().write(format!("{}\n", e).as_bytes()).unwrap();
            process::exit(1);
        }
        _ => (),
    };
}

fn wrapper() -> Result<()> {
    let cli = get_cli()?;
    run_fp_count_from_cli(&cli)
}

fn get_cli() -> Result<CLI> {
    let matches = Command::new("fp-framing")
        .version("0.1.0")
        .author("Nick Ingolia <ingolia@berkeley.edu>")
        .about("Quantifies ribosome profiling footprints on transcripts")
        .arg(
            Arg::new("output")
                .short('o')
                .long("output")
                .value_name("OUTFILE")
                .help("Output filename")
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
            Arg::new("asite")
                .short('a')
                .long("asite")
                .value_name("ASITEFILE")
                .help("A site offsets filename")
                .action(ArgAction::Set)
                .required(true),
        )
        .arg(
            Arg::new("cds_insets")
                .short('c')
                .long("cds-inset")
                .value_name("IN5',IN3'")
                .help("CDS insets for quantitation")
                .action(ArgAction::Set)
                .default_value("15,5"),
        )
        .arg(
            Arg::new("nhits")
                .short('n')
                .long("nhits")
                .help("Scale by NH of alignment")
                .action(ArgAction::SetTrue),
        )
        .arg(
            Arg::new("whole")
                .short('w')
                .long("whole")
                .help("Quantify over the whole transcript, not just the CDS")
                .action(ArgAction::SetTrue),
        )
        .arg(
            Arg::new("debug")
                .short('d')
                .long("debug")
                .help("Debug read -> feature assignment")
                .action(ArgAction::SetTrue),
        )
        .arg(
            Arg::new("reverse")
                .short('r')
                .long("reverse")
                .help("Reverse strand reads")
                .action(ArgAction::SetTrue),
        )
        .arg(Arg::new("input").value_name("INPUT.BAM").required(true))
        .get_matches();

    Ok(CLI {
        input: matches
            .get_one::<String>("input")
            .expect("input is missing")
            .to_string(),
        output: matches
            .get_one::<String>("output")
            .expect("output is missing")
            .to_string(),
        bed: matches
            .get_one::<String>("bed")
            .expect("BED file is missing")
            .to_string(),
        asite: matches
            .get_one::<String>("asite")
            .expect("A sites file is missing")
            .to_string(),
        cds_insets: matches
            .get_one::<String>("cds_insets")
            .expect("CDS insets are missing")
            .to_string(),
        // stats: matches.get_one::<String>("stats").map(|s| s.to_string()),
        // framing: matches.get_one::<String>("framing").map(|s| s.to_string()),
        nhits: matches.get_flag("nhits"),
        whole: matches.get_flag("whole"),
        debug: matches.get_flag("debug"),
        reverse: matches.get_flag("reverse"),
    })
}
