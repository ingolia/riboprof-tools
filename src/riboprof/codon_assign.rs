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

#[cfg(test)]
mod tests {
    extern crate csv;

    use super::*;

    use std::cell::*;
    use std::ops::*;
    use std::rc::*;

    use self::csv::Error;

    use bio::io::bed;
    use bio_types::annot::refids::*;

    fn transcriptome_from_str(bedstr: &str) -> Transcriptome<Rc<String>> {
        let mut refids = RefIDSet::new();
        Transcriptome::new_from_bed(bed::Reader::new(bedstr.as_bytes()).records(), &mut refids)
            .expect("Transcriptome from string")
    }

    fn trxpos_at_pos<R>(tome: &Transcriptome<R>, posstr: &str) -> Vec<(String, usize)>
    where
        R: Hash + Eq + Deref<Target = String> + From<String>,
    {
        let pos: Pos<R, ReqStrand> = posstr.parse().expect("Parsing position");
        let mut trxposns: Vec<(String, usize)> = TrxPos::transcriptome_pos(tome, &pos)
            .map(|trxpos| {
                (
                    trxpos.transcript().trxname().deref().to_string(),
                    trxpos.pos(),
                )
            })
            .collect();
        trxposns.sort();
        trxposns
    }

    #[test]
    fn transcriptome_positions() {
        let beds = "\
chr01	1000	2000	AAA	0	+	1200	1800	0	1	1000,	0,
chr01	1900	2100	BBB	0	+	1950	2050	0	1	200,	0,
chr02	1500	2500	CCC	0	+	1600	2400	0	1	1000,	0,
chr02	2100	2600	DDD	0	-	2200	2500	0	1	500,	0,
chr03	500	1500	EEE	0	+	600	1200	0	2	250,450	0,550
";
        let tome = transcriptome_from_str(&beds);

        let none: Vec<(String, usize)> = Vec::new();
        let aaa = "AAA".to_string();
        assert_eq!(
            trxpos_at_pos(&tome, "chr01:1000(+)"),
            vec![("AAA".to_string(), 0)]
        );
        assert_eq!(
            trxpos_at_pos(&tome, "chr01:1234(+)"),
            vec![("AAA".to_string(), 234)]
        );
        assert_eq!(trxpos_at_pos(&tome, "chr01:999(+)"), none);
        assert_eq!(trxpos_at_pos(&tome, "chr01:1000(-)"), none);
        assert_eq!(trxpos_at_pos(&tome, "chr01:1234(-)"), none);

        assert_eq!(
            trxpos_at_pos(&tome, "chr01:1899(+)"),
            vec![("AAA".to_string(), 899)]
        );
        assert_eq!(
            trxpos_at_pos(&tome, "chr01:1900(+)"),
            vec![("AAA".to_string(), 900), ("BBB".to_string(), 0)]
        );
        assert_eq!(
            trxpos_at_pos(&tome, "chr01:1901(+)"),
            vec![("AAA".to_string(), 901), ("BBB".to_string(), 1)]
        );
        assert_eq!(
            trxpos_at_pos(&tome, "chr01:1999(+)"),
            vec![("AAA".to_string(), 999), ("BBB".to_string(), 99)]
        );
        assert_eq!(
            trxpos_at_pos(&tome, "chr01:2000(+)"),
            vec![("BBB".to_string(), 100)]
        );
        assert_eq!(
            trxpos_at_pos(&tome, "chr01:2001(+)"),
            vec![("BBB".to_string(), 101)]
        );
        assert_eq!(
            trxpos_at_pos(&tome, "chr01:2099(+)"),
            vec![("BBB".to_string(), 199)]
        );
        assert_eq!(trxpos_at_pos(&tome, "chr01:2100(+)"), none);

        assert_eq!(
            trxpos_at_pos(&tome, "chr02:2000(+)"),
            vec![("CCC".to_string(), 500)]
        );
        assert_eq!(
            trxpos_at_pos(&tome, "chr02:2200(+)"),
            vec![("CCC".to_string(), 700)]
        );
        assert_eq!(
            trxpos_at_pos(&tome, "chr02:2400(+)"),
            vec![("CCC".to_string(), 900)]
        );

        assert_eq!(trxpos_at_pos(&tome, "chr02:2000(-)"), none);
        assert_eq!(trxpos_at_pos(&tome, "chr02:2099(-)"), none);
        assert_eq!(
            trxpos_at_pos(&tome, "chr02:2100(-)"),
            vec![("DDD".to_string(), 499)]
        );
        assert_eq!(
            trxpos_at_pos(&tome, "chr02:2200(-)"),
            vec![("DDD".to_string(), 399)]
        );
        assert_eq!(
            trxpos_at_pos(&tome, "chr02:2400(-)"),
            vec![("DDD".to_string(), 199)]
        );
        assert_eq!(
            trxpos_at_pos(&tome, "chr02:2599(-)"),
            vec![("DDD".to_string(), 0)]
        );
        assert_eq!(trxpos_at_pos(&tome, "chr02:2600(-)"), none);

        assert_eq!(trxpos_at_pos(&tome, "chr03:499(+)"), none);
        assert_eq!(
            trxpos_at_pos(&tome, "chr03:500(+)"),
            vec![("EEE".to_string(), 0)]
        );
        assert_eq!(
            trxpos_at_pos(&tome, "chr03:555(+)"),
            vec![("EEE".to_string(), 55)]
        );
        assert_eq!(
            trxpos_at_pos(&tome, "chr03:749(+)"),
            vec![("EEE".to_string(), 249)]
        );
        assert_eq!(trxpos_at_pos(&tome, "chr03:750(+)"), none);
        assert_eq!(trxpos_at_pos(&tome, "chr03:800(+)"), none);
        assert_eq!(trxpos_at_pos(&tome, "chr03:1000(+)"), none);
        assert_eq!(trxpos_at_pos(&tome, "chr03:1049(+)"), none);
        assert_eq!(
            trxpos_at_pos(&tome, "chr03:1050(+)"),
            vec![("EEE".to_string(), 250)]
        );
        assert_eq!(
            trxpos_at_pos(&tome, "chr03:1250(+)"),
            vec![("EEE".to_string(), 450)]
        );
        assert_eq!(
            trxpos_at_pos(&tome, "chr03:1450(+)"),
            vec![("EEE".to_string(), 650)]
        );
        assert_eq!(
            trxpos_at_pos(&tome, "chr03:1499(+)"),
            vec![("EEE".to_string(), 699)]
        );
        assert_eq!(trxpos_at_pos(&tome, "chr03:1500(+)"), none);
        assert_eq!(trxpos_at_pos(&tome, "chr03:2000(+)"), none);
    }
}
