use std::fs;
use std::io::Write;
use std::ops::Range;
use std::path::{Path, PathBuf};
use std::rc::Rc;
use std::str;

use anyhow::{Context, Result};
use bio::io::bed;
use bio_types::annot::refids::RefIDSet;
use rust_htslib::bam;
use rust_htslib::bam::Read as BamRead;

use crate::bam_utils::*;
use crate::parse_pair;
use crate::transcript::*;

mod framing;
mod stats;

use crate::fp_framing::framing::*;
use crate::fp_framing::stats::*;

#[derive(Debug)]
pub struct CLI {
    pub input: String,
    pub output: String,
    pub bed: String,
    pub genes: Vec<String>,
    pub flanking: String,
    pub cdsbody: String,
    pub lengths: String,
    pub count_multi: bool,
    pub annotate: Option<String>,
}

pub struct Config {
    input: String,
    output: PathBuf,
    trxome: Transcriptome<Rc<String>>,
    flanking: Range<isize>,
    cdsbody: (isize, isize),
    lengths: Range<usize>,
    count_multi: bool,
    annotate: Option<PathBuf>,
}

impl Config {
    pub fn new(cli: &CLI) -> Result<Self> {
        let trxome = Self::read_transcriptome(&cli).context("Reading transcriptome")?;

        let cdsbody_range = parse_pair(&cli.cdsbody).context("Parsing CDS body argument")?;

        Ok(Config {
            input: cli.input.to_string(),
            output: Path::new(&cli.output).to_path_buf(),
            trxome: trxome,
            flanking: parse_pair(&cli.flanking).context("Parsing flanking region argument")?,
            cdsbody: (cdsbody_range.start, cdsbody_range.end),
            lengths: parse_pair(&cli.lengths).context("Parsing length range argument")?,
            count_multi: cli.count_multi,
            annotate: cli
                .annotate
                .as_ref()
                .map(|ann| Path::new(&ann).to_path_buf()),
        })
    }

    fn output_filename(&self, suffix: &str) -> PathBuf {
        let mut name_base = self.output.file_name().map_or_else(
            || "".to_string(),
            |filename| filename.to_string_lossy().to_string(),
        );
        name_base += suffix;
        let mut filepath = self.output.clone();
        filepath.set_file_name(&name_base);
        filepath
    }

    fn read_transcriptome(cli: &CLI) -> Result<Transcriptome<Rc<String>>> {
        // ZZZ Handle Trx->Gene mappings
        let mut refids = RefIDSet::new();
        let mut trxome = Transcriptome::new();

        let mut bed_reader = bed::Reader::from_file(&cli.bed)
            .with_context(|| format!("Opening BED file {:?}", cli.bed))?;

        for (lineno, recres) in bed_reader.records().enumerate() {
            let rec = recres
                .with_context(|| format!("BED file {:?}, parsing line {}", &cli.bed, lineno + 1))?;
            let trx = Transcript::from_bed12(&rec, &mut refids)
                .with_context(|| format!("BED file {:?}, line {}", &cli.bed, lineno + 1))?;
            trxome.insert(trx)?;
        }

        Ok(trxome)
    }
}

pub fn run_fp_framing_cli(cli: CLI) -> Result<()> {
    let config = Config::new(&cli)?;
    run_fp_framing(config).with_context(|| format!("{:#?}", cli))
}

pub fn run_fp_framing(config: Config) -> Result<()> {
    let mut input = if config.input == "-" {
        bam::Reader::from_stdin()?
    } else {
        bam::Reader::from_path(Path::new(&config.input))?
    };

    let tids = {
        let mut refids: RefIDSet<Rc<String>> = RefIDSet::new();
        Tids::new(&mut refids, input.header())
    };

    // Open (empty) stats output file early to detect errors before processing data.
    let stats_filename = config.output_filename("_framing_stats.txt");
    let mut stats_file = fs::File::create(&stats_filename)
        .with_context(|| format!("Creating statistics output file {:?}", &stats_filename))?;

    let mut annotate = match config.annotate {
        Option::None => None,
        Some(ref annot_file) => {
            let header = bam::Header::from_template(input.header());
            Some(
                bam::Writer::from_path(Path::new(&annot_file), &header, bam::Format::Bam)
                    .with_context(|| format!("Opening annotated BAM file {:?}", annot_file))?,
            )
        }
    };

    let mut framing_stats = FramingStats::new(&config.lengths, &config.flanking);

    for (lineno, recres) in input.records().enumerate() {
        let mut rec = recres.with_context(|| format!("Reading record {}", lineno))?;

        let res = record_framing(
            &config.trxome,
            &tids,
            &rec,
            &config.lengths,
            &config.cdsbody,
            config.count_multi,
        )?;

        framing_stats.tally_bam_frame(&res);

        if let Some(ann_writer) = &mut annotate {
            rec.push_aux(
                b"ZF",
                bam::record::Aux::String(std::str::from_utf8(&res.aux())?),
            )?;
            ann_writer
                .write(&rec)
                .context("writing annotated BAM record")?;
        }
    }

    write!(stats_file, "{}", framing_stats.align_stats().table())
        .context("writing frame_stats.txt file")?;

    fs::write(
        config.output_filename("_frame_length.txt"),
        framing_stats.frame_length_table(),
    )
    .context("writing frame_length.txt file")?;
    fs::write(
        config.output_filename("_around_start.txt"),
        framing_stats.around_start_table(),
    )
    .context("writing around_start.txt file")?;
    fs::write(
        config.output_filename("_around_end.txt"),
        framing_stats.around_end_table(),
    )
    .context("writing around_end.txt file")?;

    Ok(())
}
