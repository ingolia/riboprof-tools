use std::error;
use std::fmt;

use failure;

use bio::io::fastq;

/// Nucleotide type in the linker, either a unique molecule identifier
/// (UMI) base or a part of the sample index.
#[derive(Debug, Clone, Copy, PartialOrd, Ord, PartialEq, Eq, Hash)]
enum LinkerNtSpec {
    UMI,
    SampleIndex,
}

impl LinkerNtSpec {
    /// Create a linker nucleotide from a specification character.
    ///
    /// # Arguments
    /// * `ch` is the specification character
    ///   * `N` specifies a UMI character
    ///   * `I` specifies a sample index character
    ///
    /// # Errors
    /// An error variant is returned for any other character.
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

/// Linker sequence specification describing how bases are removed
/// from the beginning and/or the end of the sequence and converted
/// into the UMI and the sample barcode.
#[derive(Debug, Clone, PartialOrd, Ord, PartialEq, Eq, Hash)]
pub struct LinkerSpec {
    prefix: Vec<LinkerNtSpec>,
    suffix: Vec<LinkerNtSpec>,
    sample_index_length: usize,
    umi_length: usize,
}

impl LinkerSpec {
    /// Create a new linker specification from prefix and suffix
    /// specification strings.
    ///
    /// # Arguments 
    ///
    /// * `prefix_str` describes the nucleotide prefix
    /// removed from the beginning of the sequence
    /// * `suffix_str` describes the nucleotide suffix
    /// removed from the end of the sequence
    ///
    /// # Errors
    /// An error variant is returned when any of the characters in the
    /// specification strings cannot be parsed.
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

    /// Returns the length of the prefix, the number of bases that
    /// will be removed from the beginning of the raw read
    #[allow(dead_code)]
    pub fn prefix_length(&self) -> usize {
        self.prefix.len()
    }

    /// Returns the length of the suffix, the number of bases that
    /// will be removed from the end of the raw read
    #[allow(dead_code)]
    pub fn suffix_length(&self) -> usize {
        self.suffix.len()
    }

    /// Returns the total linker length (prefix + suffix),
    /// corresponding to the total number of bases that will be
    /// removed from the raw read.
    pub fn linker_length(&self) -> usize {
        self.prefix.len() + self.suffix.len()
    }

    /// Returns the length in bases of the sample index that will be
    /// constructed from the linker
    pub fn sample_index_length(&self) -> usize {
        self.sample_index_length
    }

    /// Returns the length in bases of the UMI sequence that will be
    /// constructed from the linker
    #[allow(dead_code)]
    pub fn umi_length(&self) -> usize {
        self.umi_length
    }

    /// Split a fastq record sequence according to the linker
    /// specification. If the sequence is too short to split -- if its
    /// total length is less than the total linker length -- then
    /// `None` is returned.
    /// 
    /// # Arguments
    ///
    /// * `fq` is a FastQ record
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

/// Represents the split sequence (and quality) information from a
/// FastQ record along with the sample index and UMI sequences.
#[derive(Debug, Clone, PartialOrd, Ord, PartialEq, Eq, Hash)]
pub struct LinkerSplit<'a> {
    umi: Vec<u8>,
    sample_index: Vec<u8>,
    sequence: &'a [u8],
    quality: &'a [u8],
}

impl<'a> LinkerSplit<'a> {
    /// Returns the UMI sequence
    pub fn umi<'b>(&'b self) -> &'b [u8] {
        &self.umi
    }

    /// Returns the sample index sequence
    pub fn sample_index<'b>(&'b self) -> &'b [u8] {
        &self.sample_index
    }

    /// Returns the non-linker portion of the raw input sequence
    pub fn sequence(&self) -> &'a [u8] {
        self.sequence
    }

    /// Returns the quality information for the non-linker sequence
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

#[cfg(test)]
mod tests {
    use super::*;

    const SEQ1: &[u8] = b"ACGTACGTACGTACGT";
    const SEQ2: &[u8] = b"AAAACCCCGGGGTTTT";
    const SEQ3: &[u8] = b"CGATCGATCGATCGATCG";

    fn fastq(seq: &[u8]) -> fastq::Record
    {
        let qual: Vec<u8> = (32..(32+(seq.len() as u8))).collect();
        fastq::Record::with_attrs("test_record", None, seq, &qual)
    }

    fn assert_split(raw_seq: &[u8], prefix: &str, suffix: &str, umi: &[u8], index: &[u8], sequence: &[u8], qualstart: u8) -> ()
    {
        let rec = fastq(raw_seq);
        let spec = LinkerSpec::new(prefix, suffix).unwrap();
        let split = spec.split_record(&rec).unwrap();

        assert!(split.umi() == umi);
        assert!(split.sample_index() == index);
        assert!(split.sequence() == sequence);
        assert!(split.quality()[0] == qualstart);
        assert!(split.quality().len() == split.sequence().len());
        for i in 0..(split.sequence().len()-1) {
            assert!(split.quality()[i+1] == split.quality()[i] + 1);
        }
    }

    #[test]
    fn test_0_0() {
        assert_split(SEQ1, "", "", b"", b"", SEQ1, 0+32);
        assert_split(SEQ2, "", "", b"", b"", SEQ2, 0+32);
        assert_split(SEQ3, "", "", b"", b"", SEQ3, 0+32);
        let spec = LinkerSpec::new("", "").unwrap();
        assert!(spec.prefix_length() == 0);
        assert!(spec.suffix_length() == 0);
        assert!(spec.linker_length() == 0);
        assert!(spec.sample_index_length() == 0);
        assert!(spec.umi_length() == 0);
    }

    #[test]
    fn test_iinn_iinn() {
        assert_split(SEQ1, "IINN", "IINN", b"GTGT", b"ACAC", b"ACGTACGT", 4+32);
        assert_split(SEQ2, "IINN", "IINN", b"AATT", b"AATT", b"CCCCGGGG", 4+32);
        assert_split(SEQ3, "IINN", "IINN", b"ATCG", b"CGAT", b"CGATCGATCG", 4+32);
        let spec = LinkerSpec::new("IINN", "IINN").unwrap();
        assert!(spec.prefix_length() == 4);
        assert!(spec.suffix_length() == 4);
        assert!(spec.linker_length() == 8);
        assert!(spec.sample_index_length() == 4);
        assert!(spec.umi_length() == 4);
    }

    #[test]
    fn test_nnii_0() {
        assert_split(SEQ1, "NNII", "", b"AC", b"GT", b"ACGTACGTACGT", 4+32);
        assert_split(SEQ2, "NNII", "", b"AA", b"AA", b"CCCCGGGGTTTT", 4+32);
        assert_split(SEQ3, "NNII", "", b"CG", b"AT", b"CGATCGATCGATCG", 4+32);
        let spec = LinkerSpec::new("NNII", "").unwrap();
        assert!(spec.prefix_length() == 4);
        assert!(spec.suffix_length() == 0);
        assert!(spec.linker_length() == 4);
        assert!(spec.sample_index_length() == 2);
        assert!(spec.umi_length() == 2);
    }

    #[test]
    fn test_0_ininnini() {
        assert_split(SEQ1, "", "ININNINI", b"CTAG", b"AGCT", b"ACGTACGT", 0+32);
        assert_split(SEQ2, "", "ININNINI", b"GGTT", b"GGTT", b"AAAACCCC", 0+32);
        assert_split(SEQ3, "", "ININNINI", b"TGAC", b"ACTG", b"CGATCGATCG", 0+32);
        let spec = LinkerSpec::new("", "ININNINI").unwrap();
        assert!(spec.prefix_length() == 0);
        assert!(spec.suffix_length() == 8);
        assert!(spec.linker_length() == 8);
        assert!(spec.sample_index_length() == 4);
        assert!(spec.umi_length() == 4);
    }

    #[test]
    fn test_nnn_iii() {
        assert_split(SEQ1, "NNN", "III", b"ACG", b"CGT", b"TACGTACGTA", 3 + 32);
        assert_split(SEQ2, "NNN", "III", b"AAA", b"TTT", b"ACCCCGGGGT", 3 + 32);
        assert_split(SEQ3, "NNN", "III", b"CGA", b"TCG", b"TCGATCGATCGA", 3 + 32);
        let spec = LinkerSpec::new("NNN", "III").unwrap();
        assert!(spec.prefix_length() == 3);
        assert!(spec.suffix_length() == 3);
        assert!(spec.linker_length() == 6);
        assert!(spec.sample_index_length() == 3);
        assert!(spec.umi_length() == 3);
    }

    #[test]
    fn test_i_nnnn() {
        assert_split(SEQ1, "I", "NNNN", b"ACGT", b"A", b"CGTACGTACGT", 1+32);
        assert_split(SEQ2, "I", "NNNN", b"TTTT", b"A", b"AAACCCCGGGG", 1+32);
        assert_split(SEQ3, "I", "NNNN", b"ATCG", b"C", b"GATCGATCGATCG", 1+32);
        let spec = LinkerSpec::new("I", "NNNN").unwrap();
        assert!(spec.prefix_length() == 1);
        assert!(spec.suffix_length() == 4);
        assert!(spec.linker_length() == 5);
        assert!(spec.sample_index_length() == 1);
        assert!(spec.umi_length() == 4);
    }

    const SEQ10: &[u8] = b"ACACAGTGTG";
    const SEQ11: &[u8] = b"TGCATGCATGC";
    const SEQ12: &[u8] = b"CCCTTTGGGAAA";

    #[test]
    fn test_too_long() {
        let spec = LinkerSpec::new("IIIII", "NNNNNN").unwrap();

        let rec10 = fastq(SEQ10);
        let rec11 = fastq(SEQ11);
        let rec12 = fastq(SEQ12);

        let split10 = spec.split_record(&rec10);
        let split11 = spec.split_record(&rec11);
        let split12 = spec.split_record(&rec12);

        assert!(split10 == None);

        assert!(split11.as_ref().map(LinkerSplit::umi) == Some(b"GCATGC"));
        assert!(split11.as_ref().map(LinkerSplit::sample_index) == Some(b"TGCAT"));
        assert!(split11.as_ref().map(LinkerSplit::sequence) == Some(b""));

        assert!(split12.as_ref().map(LinkerSplit::umi) == Some(b"GGGAAA"));
        assert!(split12.as_ref().map(LinkerSplit::sample_index) == Some(b"CCCTT"));
        assert!(split12.as_ref().map(LinkerSplit::sequence) == Some(b"T"));
    }
}
