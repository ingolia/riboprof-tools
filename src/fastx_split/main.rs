extern crate failure;
#[macro_use]
extern crate clap;

extern crate bio;

use std::fs;
use std::io::{self, Read, Write};
use std::path::{Path,PathBuf};
use std::str;

use clap::{App, Arg};

use bio::io::fastq;

mod linkers;
mod sample_sheet;

use linkers::*;
use sample_sheet::*;

#[derive(Debug)]
struct Config {
    fastx_inputs: Vec<PathBuf>,
    output_dir: PathBuf,
    min_insert: usize,
    linker_spec: LinkerSpec,
    sample_map: SampleMap<Sample>,
    short_file: fastq::Writer<fs::File>,
    unknown_sample: Sample,
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
    let mut config = cli_config()?;

    for input_name in config.fastx_inputs.iter() {
        let input_reader: Box<Read> = if input_name == Path::new("-") {
            Box::new(io::stdin())
        } else {
            Box::new(fs::File::open(input_name)?)
        };

        for fqres in fastq::Reader::new(input_reader).records() {
            let fq = fqres?;

            if fq.seq().len() < config.linker_spec.linker_length() + config.min_insert {
                config.short_file.write_record(&fq)?;
            } else {
                let split = config.linker_spec.split_record(&fq)
                    .ok_or_else(|| failure::err_msg(format!("Split failed on \"{}\"", str::from_utf8(fq.seq()).unwrap_or("???"))))?;
                let sample = config.sample_map.get_mut(split.sample_index())?.unwrap_or(&mut config.unknown_sample);
                sample.handle_split_read(&fq, &split)?;
            }
        }
    }

    Ok( () )
}

fn cli_config() -> Result<Config, failure::Error> {
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

    let output_dir = PathBuf::from(matches.value_of("output_dir").unwrap());
    fs::DirBuilder::new().recursive(true).create(output_dir.as_path())?;
    
    let mut sample_map = SampleMap::new(index_length);

    let sample_sheet_txt = fs::read_to_string(matches.value_of("sample_sheet").unwrap())?;
    for (name, index) in parse_sample_sheet(&sample_sheet_txt)?.into_iter() {
        let output_file = create_fastq_writer(&output_dir, &name)?;
        let sample = Sample { name: name.to_string(), dest: output_file };
        sample_map.insert(index.into_bytes(), true, sample)?;
    }
    
    let short_file = create_fastq_writer(&output_dir, "tooshort")?;
    let unknown_sample = Sample { name: "UnknownIndex".to_string(), dest: create_fastq_writer(&output_dir, "UnknownIndex")? };

    Ok( Config { fastx_inputs: matches.values_of_lossy("input").unwrap().into_iter().map(PathBuf::from).collect(),
                 output_dir: output_dir,
                 min_insert: value_t!(matches.value_of("min_insert"), usize)?,
                 linker_spec: linker_spec,
                 sample_map: sample_map,
                 short_file: short_file,
                 unknown_sample: unknown_sample,
                 progress: None } )
}

fn create_fastq_writer(output_dir: &Path, name: &str) -> Result<fastq::Writer<fs::File>, failure::Error> {
    let mut output_path = output_dir.to_path_buf();
    output_path.push(Path::new(name));
    output_path.set_extension("fastq");
    fastq::Writer::to_file(output_path.as_path()).map_err(::std::convert::Into::into)
}

#[derive(Debug)]
struct Sample {
    name: String,
    dest: fastq::Writer<fs::File>,
}

impl Sample {
    fn handle_split_read(&mut self, fq: &fastq::Record, split: &LinkerSplit) -> Result<(), failure::Error> {
        let umi_id = format!("{}#{}", fq.id(), str::from_utf8(split.umi())?);
        let splitfq = fastq::Record::with_attrs(umi_id.as_str(), fq.desc(), split.sequence(), split.quality());
        self.dest.write_record(&splitfq)?;
        Ok( () )
    }
}
