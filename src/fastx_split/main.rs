extern crate failure;
#[macro_use]
extern crate clap;

extern crate bio;

use std::fs;
use std::io::{self, Read, Write};
use std::path::{Path, PathBuf};
use std::str;

use clap::{App, Arg};

use bio::io::fastq;

mod linkers;
mod sample;
mod sample_sheet;

use linkers::*;
use sample::*;
use sample_sheet::*;

struct Config {
    fastx_inputs: Vec<PathBuf>,
    output_dir: PathBuf,
    min_insert: usize,
    linker_spec: LinkerSpec,
    sample_map: SampleMap<Sample>,
    short_file: fastq::Writer<fs::File>,
    progress: Option<usize>,
}

fn main() {
    match run() {
        Err(e) => {
            io::stderr().write(format!("{}\n", e).as_bytes()).unwrap();
            ::std::process::exit(1);
        }
        _ => (),
    };
}

fn run() -> Result<(), failure::Error> {
    let mut config = cli_config()?;

    let mut total = 0;
    let mut tooshort = 0;

    for input_name in config.fastx_inputs.iter() {
        let input_reader: Box<Read> = if input_name == Path::new("-") {
            Box::new(io::stdin())
        } else {
            Box::new(fs::File::open(input_name)?)
        };

        for fqres in fastq::Reader::new(input_reader).records() {
            let fq = fqres?;

            total += 1;

            if fq.seq().len() < config.linker_spec.linker_length() + config.min_insert {
                config.short_file.write_record(&fq)?;
                tooshort += 1;
            } else {
                let split = config.linker_spec.split_record(&fq).ok_or_else(|| {
                    failure::err_msg(format!(
                        "Split failed on \"{}\"",
                        str::from_utf8(fq.seq()).unwrap_or("???")
                    ))
                })?;
                let mut sample = config.sample_map.get_mut(split.sample_index())?;
                sample.handle_split_read(&fq, &split)?;
            }
        }
    }

    let mut fates_path = config.output_dir.clone();
    fates_path.push("fates.txt");
    let mut fates = fs::File::create(&fates_path)?;

    for sample_rc in config.sample_map.things() {
        let sample = sample_rc.try_borrow()?;
        let mut stats_path = config.output_dir.clone();
        stats_path.push(format!("{}_stats.txt", sample.name()));
        fs::write(&stats_path, sample.stats_table())?;

        let fract = 100.0 * (sample.total() as f64) / (total as f64);
        write!(
            fates,
            "{}\t{}\t{}\t{:.2}%\n",
            sample.name(),
            str::from_utf8(sample.index())?,
            sample.total(),
            fract
        )?;
    }

    write!(
        fates,
        "short\tN/A\t{}\t{:.2}%\n",
        tooshort,
        100.0 * (tooshort as f64) / (total as f64)
    )?;

    Ok(())
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

    let linker_spec = LinkerSpec::new(
        matches.value_of("prefix").unwrap(),
        matches.value_of("suffix").unwrap(),
    )?;
    let index_length = linker_spec.sample_index_length();

    let output_dir = PathBuf::from(matches.value_of("output_dir").unwrap());
    fs::DirBuilder::new()
        .recursive(true)
        .create(output_dir.as_path())?;

    let unknown_sample = Sample::new(
        "UnknownIndex".to_string(),
        vec![b'N'; index_length],
        create_writer(&output_dir, "UnknownIndex")?,
    );

    let mut sample_map = SampleMap::new(index_length, unknown_sample);

    let sample_sheet_txt = fs::read_to_string(matches.value_of("sample_sheet").unwrap())?;
    for (name, index) in parse_sample_sheet(&sample_sheet_txt)?.into_iter() {
        let output_file = create_writer(&output_dir, &name)?;
        let sample = Sample::new(
            name.to_string(),
            index.to_string().into_bytes(),
            output_file,
        );
        sample_map.insert(index.into_bytes(), true, sample)?;
    }

    let short_file = fastq::Writer::new(create_writer(&output_dir, "tooshort")?);

    let mut mapping_file = output_dir.clone();
    mapping_file.push("mapping.txt");
    fs::write(&mapping_file, sample_map.mapping_table())?;

    Ok(Config {
        fastx_inputs: matches
            .values_of_lossy("input")
            .unwrap()
            .into_iter()
            .map(PathBuf::from)
            .collect(),
        output_dir: output_dir,
        min_insert: value_t!(matches.value_of("min_insert"), usize)?,
        linker_spec: linker_spec,
        sample_map: sample_map,
        short_file: short_file,
        progress: None,
    })
}

fn create_writer(
    output_dir: &Path,
    name: &str,
) -> Result<fs::File, failure::Error> {
    let mut output_path = output_dir.to_path_buf();
    output_path.push(Path::new(name));
    output_path.set_extension("fastq");
    fs::File::create(output_path.as_path()).map_err(::std::convert::Into::into)
}
