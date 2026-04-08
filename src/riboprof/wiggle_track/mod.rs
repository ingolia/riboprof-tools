use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};

use anyhow::{Result, anyhow};
use rust_htslib::bam;
use rust_htslib::bam::Read as BamRead;

pub mod wiggle;
pub mod window_count;

use crate::codon_assign::ASites;

pub struct CLI {
    pub input: String,
    pub output: String,
    pub asite: String,
    pub qnorm: Option<f64>,
    pub chrsizes: Option<String>,
}

impl CLI {
    fn output_names(&self) -> Result<(PathBuf, PathBuf)> {
        let out_path = Path::new(&self.output);
        let out_stem = out_path
            .file_stem()
            .ok_or_else(|| anyhow!("Bad output name {:?}", self.output))?;

        let filename = out_stem
            .to_str()
            .ok_or_else(|| anyhow!("Bad output name {:?}", self.output))?;

        let mut out_fwd = out_path.with_file_name(format!("{}_fwd", filename));
        let mut out_rev = out_path.with_file_name(format!("{}_rev", filename));

        if let Some(out_ext) = out_path.extension() {
            out_fwd.set_extension(out_ext);
            out_rev.set_extension(out_ext);
        }

        Ok((out_fwd, out_rev))
    }
}

pub struct Config {
    input: bam::Reader,
    output_fwd: Box<dyn std::io::Write>,
    output_rev: Box<dyn std::io::Write>,
    asites: ASites,
    qnorm: Option<f64>,
    chrsizes: Option<PathBuf>,
}

impl Config {
    pub fn new(cli: &CLI) -> Result<Self> {
        let input = if cli.input == "-" {
            bam::Reader::from_stdin()?
        } else {
            bam::Reader::from_path(Path::new(&cli.input))?
        };

        let (output_fwd_name, output_rev_name) = cli.output_names()?;
        let output_fwd = Box::new(fs::File::create(&output_fwd_name)?);
        let output_rev = Box::new(fs::File::create(&output_rev_name)?);

        let asites = ASites::new_from_file(&cli.asite)?;

        Ok(Config {
            input,
            output_fwd,
            output_rev,
            asites,
            qnorm: cli.qnorm,
            chrsizes: cli.chrsizes.as_ref().map(|s| Path::new(&s).to_path_buf()),
        })
    }
}

pub fn run_wiggle_track_cli(cli: &CLI) -> Result<()> {
    run_wiggle_track(&Config::new(cli)?)
}

pub fn run_wiggle_track(config: &Config) -> Result<()> {
    if let Some(chrsizes) = &config.chrsizes {
        write_chr_sizes(chrsizes, config.input.header())?;
    }

    Ok(())
}

pub fn write_chr_sizes<P: AsRef<Path>>(filename: P, hdr: &bam::HeaderView) -> Result<()> {
    let mut chrout = fs::File::create(filename)?;
    for (tid, tname) in hdr.target_names().iter().enumerate() {
        write!(
            chrout,
            "{}\t{}\n",
            String::from_utf8_lossy(tname),
            hdr.target_len(tid as u32).ok_or_else(|| anyhow!(
                "No length for {} {}",
                tid,
                String::from_utf8_lossy(tname)
            ))?
        )?;
    }
    Ok(())
}
