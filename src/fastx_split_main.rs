extern crate clap;
extern crate failure;

extern crate riboprof;

use std::io;
use std::io::Write;
use std::process;

use clap::{Arg, ArgAction, Command, value_parser};

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
    let matches = Command::new("fastx-split")
        .version("0.1.0")
        .author("Nick Ingolia <ingolia@berkeley.edu>")
        .about("Split FastQ file using index and random nucleotides")
        .arg(
            Arg::new("output_dir")
                .short('o')
                .long("output-dir")
                .value_name("OUTPUT-DIR")
                .help("Output directory name")
                .action(ArgAction::Set)
                .required(true),
        )
        .arg(
            Arg::new("min_insert")
                .value_parser(value_parser!(usize))
                .short('m')
                .long("min-insert")
                .value_name("MIN-INSERT")
                .help("Minimum insert length")
                .action(ArgAction::Set)
                .default_value("0"),
        )
        .arg(
            Arg::new("prefix")
                .short('p')
                .long("prefix")
                .value_name("PREFIX")
                .help("Prefix format string")
                .action(ArgAction::Set)
                .default_value(""),
        )
        .arg(
            Arg::new("suffix")
                .short('x')
                .long("suffix")
                .value_name("SUFFIX")
                .help("Suffix format string")
                .action(ArgAction::Set)
                .default_value(""),
        )
        .arg(
            Arg::new("sample_sheet")
                .short('s')
                .long("sample-sheet")
                .value_name("SAMPLESHEET.CSV")
                .help("File name of CSV-format sample sheet")
                .action(ArgAction::Set)
                .required(true),
        )
        .arg(
            Arg::new("progress")
                .value_parser(value_parser!(usize))
                .long("progress")
                .value_name("NSEQS")
                .help("Report progress every NSEQS sequences")
                .action(ArgAction::Set)
                .default_value("0"),
        )
        .arg(Arg::new("input").action(ArgAction::Append).required(true))
        .get_matches();

    Ok(CLI {
        fastx_inputs: matches
            .get_many::<String>("input")
            .expect("inputs")
            .map(|s| s.to_string())
            .collect(),
        output_dir: matches
            .get_one::<String>("output_dir")
            .expect("output directory is missing")
            .to_string(),
        min_insert: *matches
            .get_one::<usize>("min_insert")
            .expect("minimum insert is missing"),
        prefix: matches
            .get_one::<String>("prefix")
            .expect("prefix is missing")
            .to_string(),
        suffix: matches
            .get_one::<String>("suffix")
            .expect("suffix is missing")
            .to_string(),
        sample_sheet: matches
            .get_one::<String>("sample_sheet")
            .expect("sample sheet is missing")
            .to_string(),
        progress: *matches
            .get_one::<usize>("progress")
            .expect("progress interval is missing"),
    })
}
