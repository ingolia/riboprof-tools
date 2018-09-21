use std::fs;
use std::io::{self, Read, Write};
use std::path::{Path, PathBuf};
use std::str;

use failure;

use bio::io::fastq;

mod linkers;
mod sample;
mod sample_sheet;

use fastx_split::linkers::*;
use fastx_split::sample::*;
use fastx_split::sample_sheet::*;

pub struct CLI {
    pub fastx_inputs: Vec<String>,
    pub output_dir: String,
    pub min_insert: usize,
    pub prefix: String,
    pub suffix: String,
    pub sample_sheet: String,
    pub progress: usize,
}

pub struct Config {
    fastx_inputs: Vec<PathBuf>,
    output_dir: PathBuf,
    min_insert: usize,
    linker_spec: LinkerSpec,
    sample_map: SampleMap<Sample>,
    short_file: fastq::Writer<fs::File>,
    progress: Option<usize>,
}

impl Config {
    pub fn new(cli: &CLI) -> Result<Self, failure::Error> {
        let linker_spec = LinkerSpec::new(&cli.prefix, &cli.suffix)?;
        let index_length = linker_spec.sample_index_length();

        let output_dir = Path::new(&cli.output_dir).to_path_buf();
        fs::DirBuilder::new()
            .recursive(true)
            .create(output_dir.as_path())?;

        let unknown_sample = Sample::new(
            "UnknownIndex".to_string(),
            vec![b'N'; index_length],
            Config::create_writer(&output_dir, "UnknownIndex")?,
        );

        let mut sample_map = SampleMap::new(index_length, unknown_sample);

        let sample_sheet_txt = fs::read_to_string(&cli.sample_sheet)?;
        for (name, index) in parse_sample_sheet(&sample_sheet_txt)?.into_iter() {
            let output_file = Config::create_writer(&output_dir, &name)?;
            let sample = Sample::new(
                name.to_string(),
                index.to_string().into_bytes(),
                output_file,
            );
            sample_map.insert(index.into_bytes(), true, sample)?;
        }

        let short_file = fastq::Writer::new(Config::create_writer(&output_dir, "tooshort")?);

        let mut mapping_file = output_dir.clone();
        mapping_file.push("mapping.txt");
        fs::write(&mapping_file, sample_map.mapping_table())?;

        Ok(Config {
            fastx_inputs: cli.fastx_inputs.iter().map(PathBuf::from).collect(),
            output_dir: output_dir,
            min_insert: cli.min_insert,
            linker_spec: linker_spec,
            sample_map: sample_map,
            short_file: short_file,
            progress: if cli.progress > 0 {
                Some(cli.progress)
            } else {
                None
            },
        })
    }

    fn create_writer(output_dir: &Path, name: &str) -> Result<fs::File, failure::Error> {
        let mut output_path = output_dir.to_path_buf();
        output_path.push(Path::new(name));
        output_path.set_extension("fastq");
        fs::File::create(output_path.as_path()).map_err(::std::convert::Into::into)
    }
}

pub fn split_file<P: AsRef<Path>>(
    config: &mut Config,
    input_name: P,
) -> Result<(usize, usize), failure::Error> {
    let mut total = 0;
    let mut tooshort = 0;

    let input_reader: Box<Read> = if input_name.as_ref() == Path::new("-") {
        Box::new(io::stdin())
    } else {
        Box::new(fs::File::open(&input_name)?)
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

        if config.progress.map_or(false, |nprog| total % nprog == 0) {
            print!(
                "{:7} reads from {}\n",
                total,
                input_name.as_ref().to_str().unwrap_or("???")
            );
        }
    }

    Ok((total, tooshort))
}

pub fn write_stats(config: &Config, total: usize, tooshort: usize) -> Result<(), failure::Error> {
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

pub fn fastx_split(mut config: Config) -> Result<(), failure::Error> {
    let mut total = 0;
    let mut tooshort = 0;

    for input_name in config.fastx_inputs.to_vec() {
        let (file_total, file_tooshort) = split_file(&mut config, input_name)?;
        total += file_total;
        tooshort += file_tooshort;
    }

    write_stats(&config, total, tooshort)?;

    Ok(())
}
