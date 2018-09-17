use std::error;
use std::fmt;

use failure;

use bio::io::fastq;

#[derive(Debug, Clone, Copy, PartialOrd, Ord, PartialEq, Eq, Hash)]
enum LinkerNtSpec {
    UMI,
    SampleIndex,
}

impl LinkerNtSpec {
    pub fn new(ch: char) -> Result<Self, failure::Error> {
        match ch {
            'N' => Ok(LinkerNtSpec::UMI),
            'I' => Ok(LinkerNtSpec::SampleIndex),
            _ => Err(LinkerError::BadSpecChar(ch).into()),
        }
    }
}

impl fmt::Display for LinkerNtSpec {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            LinkerNtSpec::UMI => write!(f, "N"),
            LinkerNtSpec::SampleIndex => write!(f, "I"),
        }
    }
}

#[derive(Debug, Clone, PartialOrd, Ord, PartialEq, Eq, Hash)]
pub struct LinkerSpec {
    prefix: Vec<LinkerNtSpec>,
    suffix: Vec<LinkerNtSpec>,
    sample_index_length: usize,
    umi_length: usize,
}

impl LinkerSpec {
    pub fn new(prefix_str: &str, suffix_str: &str) -> Result<Self, failure::Error> {
        let prefix_res: Result<Vec<LinkerNtSpec>, failure::Error> =
            prefix_str.chars().map(LinkerNtSpec::new).collect();
        let suffix_res: Result<Vec<LinkerNtSpec>, failure::Error> =
            suffix_str.chars().map(LinkerNtSpec::new).collect();

        let prefix = prefix_res?;
        let suffix = suffix_res?;

        let sample_index_length = prefix
            .iter()
            .chain(suffix.iter())
            .filter(|&nt| *nt == LinkerNtSpec::SampleIndex)
            .count();
        let umi_length = prefix
            .iter()
            .chain(suffix.iter())
            .filter(|&nt| *nt == LinkerNtSpec::UMI)
            .count();

        Ok(LinkerSpec {
            prefix: prefix,
            suffix: suffix,
            sample_index_length: sample_index_length,
            umi_length: umi_length,
        })
    }

    #[allow(dead_code)]
    pub fn prefix_length(&self) -> usize {
        self.prefix.len()
    }

    #[allow(dead_code)]
    pub fn suffix_length(&self) -> usize {
        self.suffix.len()
    }

    pub fn linker_length(&self) -> usize {
        self.prefix.len() + self.suffix.len()
    }

    pub fn sample_index_length(&self) -> usize {
        self.sample_index_length
    }

    pub fn umi_length(&self) -> usize {
        self.umi_length
    }

    pub fn split_record<'a>(&self, fq: &'a fastq::Record) -> Option<LinkerSplit<'a>> {
        let sequence = fq.seq();

        if sequence.len() >= self.prefix.len() + self.suffix.len() {
            let mut umi = Vec::new();
            let mut sample_index = Vec::new();

            for i in 0..self.prefix.len() {
                match self.prefix[i] {
                    LinkerNtSpec::UMI => umi.push(sequence[i]),
                    LinkerNtSpec::SampleIndex => sample_index.push(sequence[i]),
                };
            }

            let suffix_start = sequence.len() - self.suffix.len();
            for i in 0..self.suffix.len() {
                match self.suffix[i] {
                    LinkerNtSpec::UMI => umi.push(sequence[suffix_start + i]),
                    LinkerNtSpec::SampleIndex => sample_index.push(sequence[suffix_start + i]),
                };
            }

            Some(LinkerSplit {
                umi: umi,
                sample_index: sample_index,
                sequence: &sequence[self.prefix.len()..suffix_start],
                quality: &fq.qual()[self.prefix.len()..suffix_start],
            })
        } else {
            None
        }
    }
}

impl fmt::Display for LinkerSpec {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "prefix: ")?;
        for nt in self.prefix.iter() {
            nt.fmt(f)?;
        }
        write!(f, ", suffix: ")?;
        for nt in self.suffix.iter() {
            nt.fmt(f)?;
        }
        Ok(())
    }
}

#[derive(Debug, Clone, PartialOrd, Ord, PartialEq, Eq, Hash)]
pub struct LinkerSplit<'a> {
    umi: Vec<u8>,
    sample_index: Vec<u8>,
    sequence: &'a [u8],
    quality: &'a [u8],
}

impl<'a> LinkerSplit<'a> {
    pub fn umi<'b>(&'b self) -> &'b [u8] {
        &self.umi
    }

    pub fn sample_index<'b>(&'b self) -> &'b [u8] {
        &self.sample_index
    }

    pub fn sequence(&self) -> &'a [u8] {
        self.sequence
    }

    pub fn quality(&self) -> &'a [u8] {
        self.quality
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum LinkerError {
    BadSpecChar(char),
}

impl fmt::Display for LinkerError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            LinkerError::BadSpecChar(ch) => write!(f, "Bad linker spec char \'{}\'", ch),
        }
    }
}

impl error::Error for LinkerError {}
