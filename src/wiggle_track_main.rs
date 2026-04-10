use std::io;
use std::io::Write;
use std::process;

use clap::{Arg, ArgAction, Command, value_parser};

use anyhow::Result;
use riboprof::wiggle_track::{CLI, run_wiggle_track_cli};

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
    run_wiggle_track_cli(&cli)
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
            Arg::new("asite")
                .short('a')
                .long("asite")
                .value_name("ASITEFILE")
                .help("A site offsets filename")
                .action(ArgAction::Set)
                .required(true),
        )
        .arg(
            Arg::new("qnorm")
                .value_parser(value_parser!(f64))
                .short('q')
                .long("qnorm")
                .value_name("QNORM")
                .help("Multiplicative scaling factor")
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("chrsizes")
                .short('x')
                .long("chrsizes")
                .help("Chromosome size output file")
                .action(ArgAction::Set),
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
        asite: matches
            .get_one::<String>("asite")
            .expect("A sites file is missing")
            .to_string(),
        qnorm: matches.get_one::<f64>("qnorm").copied(),
        chrsizes: matches.get_one::<String>("chrsizes").map(|s| s.to_string()),
    })
}
