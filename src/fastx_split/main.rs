extern crate failure;
#[macro_use]
extern crate clap;

extern crate bio;

use std::error::Error;
use std::fs;
use std::io::{self, Read, Write};
use std::path::{Path,PathBuf};

use clap::{App, Arg, ArgMatches};

mod linkers;
mod sample_sheet;

use linkers::LinkerSpec;
use sample_sheet::*;

struct Config {
    fastx_inputs: Vec<PathBuf>,
    output_dir: PathBuf,
    min_insert: usize,
    linker_spec: LinkerSpec,
    sample_map: SampleMap<fs::File>,
    progress: Option<usize>,
}

fn main() {
    match run() {
        Err(e) => { io::stderr().write(format!("{}\n", e).as_bytes()).unwrap();
                    ::std::process::exit(1);
        },
        _ => ()
    };
}

fn run() -> Result<(), failure::Error> {
    let config = cli_matches()?;
    Ok( () )
}

fn cli_matches() -> Result<Config, failure::Error> {
    let matches = App::new("fastx-split")
        .version("0.1.0")
        .author("Nick Ingolia <ingolia@berkeley.edu>")
        .about("Split FastQ file using index and random nucleotides")
        .arg(
            Arg::with_name("output_dir")
                .short("o")
                .long("output-dir")
                .value_name("OUTPUT-DIR")
                .help("Output directory name")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("min_insert")
                .short("m")
                .long("min-insert")
                .value_name("MIN-INSERT")
                .help("Minimum insert length")
                .takes_value(true)
                .default_value("0"),
        )
        .arg(
            Arg::with_name("prefix")
                .short("p")
                .long("prefix")
                .value_name("PREFIX")
                .help("Prefix format string")
                .takes_value(true)
                .default_value(""),
        )
        .arg(
            Arg::with_name("suffix")
                .short("x")
                .long("suffix")
                .value_name("SUFFIX")
                .help("Suffix format string")
                .takes_value(true)
                .default_value(""),
        )
        .arg(
            Arg::with_name("sample_sheet")
                .short("s")
                .long("sample-sheet")
                .value_name("SAMPLESHEET.CSV")
                .help("File name of CSV-format sample sheet")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("progress")
                .long("progress")
                .value_name("NSEQS")
                .help("Report progress every NSEQS sequences")
                .takes_value(true),
        )
        .arg(Arg::with_name("input").multiple(true).required(true))
        .get_matches();

    let linker_spec = LinkerSpec::new(matches.value_of("prefix").unwrap(), matches.value_of("suffix").unwrap())?;

    let index_length = linker_spec.sample_index_length();

    let sample_sheet_txt = fs::read_to_string(matches.value_of("sample_sheet").unwrap())?;
    
    let mut sample_map = SampleMap::new(index_length);

    let output_dir = PathBuf::from(matches.value_of("output_dir").unwrap());

    fs::DirBuilder::new().recursive(true).create(output_dir.as_path())?;

    for (name, index) in parse_sample_sheet(&sample_sheet_txt)?.into_iter() {
        let mut output_path = output_dir.clone();
        output_path.push(Path::new(&name));
        output_path.set_extension("fastq");
        let output_file = fs::File::create(output_path.as_path())?;
        sample_map.insert(index.into_bytes(), true, name.to_string(), output_file)?;
    }

    Ok( Config { fastx_inputs: matches.values_of_lossy("input").unwrap().into_iter().map(PathBuf::from).collect(),
                 output_dir: output_dir,
                 min_insert: value_t!(matches.value_of("min_insert"), usize)?,
                 linker_spec: linker_spec,
                 sample_map: sample_map,
                 progress: None } )
}
