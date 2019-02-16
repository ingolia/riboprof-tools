use std::error::Error;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};

use failure;

use rust_htslib::bam;
use rust_htslib::bam::Read as BamRead;

mod read_class;
mod read_group;
mod stats;

use bam_suppress_duplicates::read_class::*;
use bam_suppress_duplicates::read_group::*;
use bam_suppress_duplicates::stats::*;

pub struct CLI {
    pub bam_input: String,
    pub bam_output: String,
    pub bam_dups: Option<String>,
    pub stats: Option<String>,
    pub annotate: bool,
}

pub struct Config {
    input: bam::Reader,
    uniq_output: bam::Writer,
    dups_output: Option<bam::Writer>,
    stat_file: Option<PathBuf>,
    annotate: bool,
    stats: Stats,
}

const DEFAULT_NLIM: usize = 100; // ZZZ

impl Config {
    pub fn new(cli: &CLI) -> Result<Self, failure::Error> {
        let input = if cli.bam_input == "-" {
            bam::Reader::from_stdin()?
        } else {
            bam::Reader::from_path(Path::new(&cli.bam_input))?
        };

        let header = bam::Header::from_template(input.header());
        let uniq_out = if cli.bam_output == "-" {
            bam::Writer::from_stdout(&header)?
        } else {
            bam::Writer::from_path(Path::new(&cli.bam_output), &header)?
        };
        
        let dups_out = match cli.bam_dups {
            None => None,
            Some(ref dups_file) => {
                Some(bam::Writer::from_path(Path::new(&dups_file), &header)?)
            }
        };
        
        let stats = Stats::new(DEFAULT_NLIM);
        
        Ok(Config {
            input: input,
            uniq_output: uniq_out,
            dups_output: dups_out,
            stat_file: cli.stats.as_ref().map(|s| Path::new(&s).to_path_buf()),
            annotate: cli.annotate,
            stats: stats,
        })
    }
}

pub fn same_location(r1: &bam::Record, r2: &bam::Record) -> bool {
    (r1.tid() == r2.tid()) &&
        (r1.pos() == r2.pos()) &&
        (r1.is_reverse() == r2.is_reverse())
}

pub fn read_tag(r1: &bam::Record) -> Option<&[u8]> {
    if let Some(delim_pos) = r1.qname().iter().position(|&ch| ch == b'#') {
        Some(r1.qname().split_at(delim_pos + 1).1)
    } else {
        None
    }
}

// N.B. No read tag => never a duplicate!
pub fn same_tag(r0: &bam::Record, r1: &bam::Record) -> bool {
    if let Some(tag0) = read_tag(r0) {
        if let Some(tag1) = read_tag(r1) {
            (tag0 == tag1)
        } else {
            false
        }
    } else {
        false
    }
}

pub fn same_cigar(r0: &bam::Record, r1: &bam::Record) -> bool {
    r0.raw_cigar() == r1.raw_cigar()
}

pub fn bam_suppress_duplicates(mut config: Config) -> Result<(), failure::Error> {
    let loc_groups = ReadGroups::new(&same_location, &mut config.input)?;

    for loc_group_res in loc_groups {
        let loc_group = loc_group_res?;
        let mut cigar_classes = ReadClass::new(&same_cigar);
        cigar_classes.insert_all(loc_group.into_iter());
        for cigar_class in cigar_classes.classes() {
            let mut tag_classes = ReadClass::new(&same_tag);
            tag_classes.insert_all(cigar_class.into_iter());

            let mut n_total = 0;
            let mut n_unique = 0;
            
            for mut tag_class in tag_classes.classes() {
                if read_tag(tag_class.first().unwrap()).is_none() {
                    assert!(tag_class.len() == 1);
                    config.uniq_output.write(tag_class.first().unwrap())?;
                    config.stats.tally_untagged();
                } else {
                    let tag_class_len = tag_class.len();
                    n_total += tag_class_len;
                    n_unique += 1;

                    let (mut uniq, dups) = tag_class.split_first_mut().unwrap();

                    if config.annotate {
                        uniq.push_aux(b"ZD", &bam::record::Aux::Integer(tag_class_len as i64))?;
                    }

                    config.uniq_output.write(uniq)?;
                    if let Some(ref mut out) = config.dups_output.as_mut() {
                        for dup in dups {
                            out.write(dup)?;
                        }
                    }
                }

                config.stats.tally(n_total, n_unique);
            }
        }
    }
    
    if let Some(ref stats_file) = config.stat_file {
        let mut stats_out = fs::File::create(stats_file)?;
        stats_out.write_all(config.stats.dedup_table().as_bytes())?;
    }
    
    Ok(())
}
