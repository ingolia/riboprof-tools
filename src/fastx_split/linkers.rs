use std::error;
use std::fmt;

use failure;

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
}

impl LinkerSpec {
    pub fn new(prefix_str: &str, suffix_str: &str) -> Result<Self, failure::Error> {
        let prefix: Result<Vec<LinkerNtSpec>, failure::Error> = prefix_str.chars().map(LinkerNtSpec::new).collect();
        let suffix: Result<Vec<LinkerNtSpec>, failure::Error> = suffix_str.chars().map(LinkerNtSpec::new).collect();

        Ok( LinkerSpec { prefix: prefix?, suffix: suffix? } )
    }

    pub fn prefix_length(&self) -> usize { self.prefix.len() }
    
    pub fn suffix_length(&self) -> usize { self.suffix.len() }

    pub fn linker_length(&self) -> usize { self.prefix.len() + self.suffix.len() }

    pub fn sample_index_length(&self) -> usize {
        self.prefix
            .iter()
            .chain(self.suffix.iter())
            .filter(|&nt| *nt == LinkerNtSpec::SampleIndex)
            .count()
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
        Ok( () )
    }
}

#[derive(Debug, Clone, PartialOrd, Ord, PartialEq, Eq, Hash)]
struct LinkerSplit<'a> {
    umi: Vec<u8>,
    sample_index: Vec<u8>,
    sequence: &'a [u8],
}

impl<'a> LinkerSplit<'a> {
    pub fn new(spec: &LinkerSpec, sequence: &'a [u8]) -> Option<Self> {
        if sequence.len() >= spec.prefix.len() + spec.suffix.len() {
            let mut umi = Vec::new();
            let mut sample_index = Vec::new();

            for i in 0..spec.prefix.len() {
                match spec.prefix[i] {
                    LinkerNtSpec::UMI => umi.push(sequence[i]),
                    LinkerNtSpec::SampleIndex => sample_index.push(sequence[i]),
                };
            }

            let suffix_start = sequence.len() - spec.suffix.len();
            for i in 0..spec.suffix.len() {
                match spec.suffix[i] {
                    LinkerNtSpec::UMI => umi.push(sequence[suffix_start + i]),
                    LinkerNtSpec::SampleIndex => umi.push(sequence[suffix_start + i]),
                };
            }

            Some( LinkerSplit { umi: umi, 
                                sample_index: sample_index,
                                sequence: &sequence[spec.prefix.len()..suffix_start] } )
        } else {
            None
        }
    }

    pub fn umi(&'a self) -> &'a [u8] {
        &self.umi
    }

    pub fn sample_index(&'a self) -> &'a [u8] {
        &self.sample_index
    }

    pub fn sequence(&self) -> &'a [u8] {
        self.sequence
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

impl error::Error for LinkerError { }
