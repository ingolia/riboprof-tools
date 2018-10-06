use std::error::Error;
use std::fmt;
use std::fs;
use std::hash::Hash;
use std::num::ParseIntError;
use std::path::Path;
use std::str::FromStr;

use failure;
use regex::Regex;

use bio_types::annot::loc::*;
use bio_types::annot::pos::*;
use bio_types::strand::*;

use transcript::*;

pub struct TrxPos<'a, R: 'a> {
    transcript: &'a Transcript<R>,
    pos: usize,
}

impl<'a, R: 'a> TrxPos<'a, R> {
    pub fn new(transcript: &'a Transcript<R>, pos: usize) -> Self {
        TrxPos {
            transcript: transcript,
            pos: pos,
        }
    }

    pub fn transcript(&self) -> &'a Transcript<R> {
        self.transcript
    }
    pub fn pos(&self) -> usize {
        self.pos
    }

    pub fn offset_from_trx_start(&self) -> isize {
        self.pos as isize
    }
    pub fn offset_from_trx_end(&self) -> isize {
        self.transcript.loc().length() as isize - self.pos as isize
    }

    pub fn offset_from_cds_start(&self) -> Option<isize> {
        self.transcript
            .cds_range()
            .as_ref()
            .map(|cds| cds.start as isize - self.pos as isize)
    }

    pub fn offset_from_cds_end(&self) -> Option<isize> {
        self.transcript
            .cds_range()
            .as_ref()
            .map(|cds| cds.end as isize - self.pos as isize)
    }
}

impl<'a, R: 'a + Eq> TrxPos<'a, R> {
    pub fn from_genomic_pos<'b>(
        transcript: &'a Transcript<R>,
        gpos: &'b Pos<R, ReqStrand>,
    ) -> Option<Self> {
        transcript
            .loc()
            .pos_into(gpos)
            .map_or(None, |tpos| match tpos.strand() {
                ReqStrand::Forward => Some(tpos.pos()),
                ReqStrand::Reverse => None,
            })
            .map(|pos| Self::new(transcript, pos as usize))
    }
}

impl<'a, R: Eq + Hash> TrxPos<'a, R> {
    pub fn transcriptome_pos<'b, 'c>(tome: &'b Transcriptome<R>, gpos: &'c Pos<R, ReqStrand>) -> impl Iterator<Item = TrxPos<'a, R>>
        where 'b: 'a, 'c: 'a
    {
        tome.find_at_loc(gpos).filter_map(move |trx| Self::from_genomic_pos(trx, gpos))
    }
}

/// Mapping of A site positions within a footprint, based on fragment
/// length.
#[derive(Debug, Clone)]
pub struct ASites {
    a_site_offsets: Vec<Option<usize>>,
}

impl ASites {
    /// Construct an `ASites` mapping based on a table of A site
    /// offsets. The file is a tab-delimited table of length / offset
    /// pairs, e.g.,
    ///
    /// `27      14`
    /// `28      15`
    ///
    /// # Arguments
    ///
    /// `path` specifies the path for the file to read
    ///
    /// # Errors
    ///
    /// An error variant is returned when an `io::Error` arises
    /// reading the file or when an `ASiteParseError` arises in
    /// parsing it.
    pub fn new_from_file<P: AsRef<Path>>(path: P) -> Result<Self, failure::Error> {
        Self::from_str(&fs::read_to_string(path)?).map_err(|e| e.into())
    }

    /// Returns the A site offset within a footprint fragment as a
    /// function of its length.
    ///
    /// # Arguments
    ///
    /// `len` is the fragment length
    ///
    /// ```
    /// # use riboprof::codon_assign::*;
    /// # use riboprof::codon_assign::ASiteParseError;
    /// # fn try_main() -> Result<(), Box<ASiteParseError>> {
    /// let asites = "27\t14\n28\t15\n".parse::<ASites>()?;
    /// assert_eq!(asites.offset(26), None);
    /// assert_eq!(asites.offset(27), Some(14));
    /// assert_eq!(asites.offset(28), Some(15));
    /// assert_eq!(asites.offset(29), None);
    /// # Ok(())
    /// # }
    /// # fn main() { try_main().unwrap(); }
    /// ```
    pub fn offset(&self, len: usize) -> Option<usize> {
        self.a_site_offsets
            .get(len)
            .unwrap_or(&None)
            .as_ref()
            .map(|o| *o)
    }

    /// Returns the A site position from a footprint location.
    ///
    /// # Arguments
    ///
    /// `fp` is the location of a footprint fragment
    /// ```
    /// # extern crate bio_types;
    /// # extern crate riboprof;
    /// # use std::error::Error;
    /// # use riboprof::codon_assign::*;
    /// # use riboprof::codon_assign::ASiteParseError;
    /// # use bio_types::annot::pos::*;
    /// # use bio_types::annot::contig::*;
    /// # use bio_types::strand::*;
    /// # fn try_main() -> Result<(), Box<Error>> {
    /// let asites = "27\t14\n28\t15\n".parse::<ASites>()?;
    /// let fp1 = "chr2:300000-300027(+)".parse::<Contig<String,ReqStrand>>()?;
    /// assert_eq!(asites.a_site(fp1), Some("chr2:300014(+)".parse::<Pos<String,ReqStrand>>()?));
    /// let fp2 = "chr3:200000-200023(-)".parse::<Contig<String,ReqStrand>>()?;
    /// assert_eq!(asites.a_site(fp2), None);
    /// # Ok(())
    /// # }
    /// # fn main() { try_main().unwrap(); }
    /// ```
    pub fn a_site<L>(&self, fp: L) -> Option<Pos<L::RefID, ReqStrand>>
    where
        L: Loc,
        L::Strand: Into<ReqStrand> + Copy,
        L::RefID: Clone,
    {
        match self.offset(fp.length()) {
            Some(offset) => fp.pos_outof(&Pos::new((), offset as isize, ReqStrand::Forward)),
            None => None,
        }
    }
}

impl FromStr for ASites {
    type Err = ASiteParseError;

    fn from_str(table: &str) -> Result<Self, Self::Err> {
        let mut offsets = Vec::new();
        let re = Regex::new("^(\\d+)\t(\\d+)$").unwrap();

        for line in table.lines().map(str::trim_right) {
            let cap = re
                .captures(line)
                .ok_or_else(|| ASiteParseError::BadLine(line.to_string()))?;
            let len = cap[1]
                .parse::<usize>()
                .map_err(|e| ASiteParseError::BadLength(e, cap[1].to_owned()))?;
            let off = cap[2]
                .parse::<usize>()
                .map_err(|e| ASiteParseError::BadOffset(e, cap[1].to_owned()))?;

            while offsets.len() <= len {
                offsets.push(None);
            }
            offsets[len] = Some(off);
        }

        Ok(ASites {
            a_site_offsets: offsets,
        })
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ASiteParseError {
    BadLine(String),
    BadLength(ParseIntError, String),
    BadOffset(ParseIntError, String),
}

impl Error for ASiteParseError {}

impl fmt::Display for ASiteParseError {
    fn fmt(&self, f: &mut fmt::Formatter) -> Result<(), fmt::Error> {
        match self {
            ASiteParseError::BadLine(line) => write!(f, "Bad A sites line: \"{}\"", &line),
            ASiteParseError::BadLength(err, line) => {
                write!(f, "Error parsing length \"{}\": {}", line, err)
            }
            ASiteParseError::BadOffset(err, line) => {
                write!(f, "Error parsing offset \"{}\": {}", line, err)
            }
        }
    }
}
