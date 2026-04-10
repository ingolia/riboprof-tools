use std::ops::Range;
use std::rc::Rc;

use anyhow::{Result, anyhow};
use bio_types::annot::loc::Loc;
use bio_types::annot::spliced::Spliced;
use bio_types::strand::ReqStrand;
use rust_htslib::bam;

use crate::bam_utils::{Tids, bam_to_spliced};
use crate::codon_assign::ASites;
use crate::transcript::{Transcript, splice_compatible};

#[derive(PartialEq, Eq, PartialOrd, Ord, Debug, Hash)]
pub enum Compat {
    IncompatStrand,
    IncompatSplice,
    IncompatASite,
    IncompatNoCDS,
    IncompatUpstream,
    IncompatDownstream,
    Compat,
}

pub struct CountConfig {
    asites: ASites,
    cds_insets: Range<isize>,
    nhits: bool,
    whole: bool,
    strand: ReqStrand,
}

impl CountConfig {
    /// Creates a read counting configuration
    ///
    /// # Arguments
    /// * `asites` is the length-based A-site assignment table
    /// * `cds_insets` are the codons excluded from the start and end of the CDS
    /// * `nhits` controls whether multimapping reads are scaled based on NH
    /// * `whole` controls whether reads are mapped to the whole transcript,
    ///   versus just the CDS
    /// * `reverse` indicates that footprints align to the reverse strand
    pub fn new(
        asites: ASites,
        cds_insets: Range<isize>,
        nhits: bool,
        whole: bool,
        reverse: bool,
    ) -> Self {
        CountConfig {
            asites,
            cds_insets,
            nhits,
            whole,
            strand: if reverse {
                ReqStrand::Reverse
            } else {
                ReqStrand::Forward
            },
        }
    }

    /// Returns the count weight for a BAM alignment record
    ///
    /// If `nhits` is set on `CountConfig`, the count weight is,
    /// `1.0 / (NH)`, where `NH` is the value in the NH aux field.
    /// An error is returned if NH is missing or non-positive.
    ///
    /// If `nhits` is not set, the count weight is always 1.0
    /// # Arguments
    /// * `rec` is the record to be counted
    pub fn record_weight(&self, rec: &bam::Record) -> Result<f64> {
        if self.nhits {
            let aux_nh = rec.aux(b"NH")?;
            if let bam::record::Aux::I32(nh) = aux_nh {
                if nh < 1 {
                    Err(anyhow!("Bad NH = {} < 1", nh))
                } else {
                    Ok((nh as f64).recip())
                }
            } else {
                Err(anyhow!("Bad NH record {:?}", aux_nh))
            }
        } else {
            Ok(1.0)
        }
    }

    /// Compatibility of the BAM alignment with transcript
    ///
    /// The BAM alignment is convered into a genomic location and
    /// used to determine the compatibility of the footprint with
    /// the transcript as per `fp_compat()`.
    /// # Arguments
    /// * `trx` is the transcript annotation
    /// * `tids` is the BAM header target ID-to-name mapping
    /// * `rec` is the BAM alignment record of the footprint
    pub fn bam_compat(
        &self,
        trx: &Transcript<Rc<String>>,
        tids: &Tids<Rc<String>>,
        rec: &bam::Record,
    ) -> Result<Compat> {
        if let Some(fp) = bam_to_spliced(tids, &rec)? {
            Ok(self.fp_compat(trx, fp))
        } else {
            Ok(Compat::IncompatSplice)
        }
    }

    /// Compatibility of a footprint with the transcript
    ///
    /// * If the footprint location is on the opposite strand from the transcript,
    ///   then `Compat::IncompatStrand` is returned
    /// * If the splicing structure of the footprint location is not compatible
    ///   with the transcript, then `Compat::IncompatSplice` is returned.
    /// * If the A site cannot be determined from the footprint length,
    ///   then `Comapt::IncompatASite` is returned.
    /// Otherwise, the A site is used to determine compatibility as per `asite_compat()`.
    /// # Arguments
    /// * `trx` is the transcript
    /// * `fp` is the genomic location of the footprint alignment
    pub fn fp_compat(
        &self,
        trx: &Transcript<Rc<String>>,
        fp: Spliced<Rc<String>, ReqStrand>,
    ) -> Compat {
        let strand = fp.strand().on_strand(self.strand);
        let fp_str = fp.into_stranded(strand);

        if splice_compatible(&trx.loc(), &fp_str) {
            if let Some(asite) = self.asites.a_site(fp_str) {
                let pos = trx
                    .loc()
                    .pos_into(&asite)
                    .expect("pos_into(asite) failed after splice_compatible(...) = true")
                    .pos();
                self.asite_compat(&trx.cds_range(), pos)
            } else {
                Compat::IncompatASite
            }
        } else {
            if trx.loc().strand() != fp_str.strand() {
                Compat::IncompatStrand
            } else {
                Compat::IncompatSplice
            }
        }
    }

    /// Compatibility of A site with transcript
    ///
    /// If `whole` is set on `CountConfig`, all A sites are considered `Compat::Compat`.
    ///
    /// If `whole` is not set and `cds_range` is `None`, all A sites are
    /// considered `Compat::IncompatNoCDS`.
    ///
    /// Otherwise:
    /// * The nucleotide position of the A site must fall downstream of the first
    ///   nucleotide in the CDS by at least as many _codons_ as the CDS insets start.
    ///   If not, the A site is considered `Compat::IncompatUpstream`.
    /// * It must fall upstream of the last nucleotide of the CDS by at least as many
    ///   _codons_ as the CDS insets end. If not, it is `Compat::IncompatDownstream`.
    /// If both criteria are satisfied, it is `Compat::Compat`.
    ///
    /// # Arguments
    /// * `cds_range` is the CDS annotated on the transcript
    /// * `asite` is the position of the A site on the transcript
    pub fn asite_compat(&self, cds_range: &Option<Range<usize>>, asite: isize) -> Compat {
        if self.whole {
            Compat::Compat
        } else {
            if let Some(cds) = cds_range {
                if (asite - cds.start as isize) < (self.cds_insets.start * 3) {
                    Compat::IncompatUpstream
                } else if (cds.end as isize - asite) < (self.cds_insets.end * 3) {
                    Compat::IncompatDownstream
                } else {
                    Compat::Compat
                }
            } else {
                Compat::IncompatNoCDS
            }
        }
    }
}

#[derive(Clone, Debug)]
pub struct Count {
    compat: f64,
    strand: f64,
    upstream: f64,
    downstream: f64,
    asite: f64,
    incompat: f64,
}

impl Count {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn compat(&self) -> f64 {
        self.compat
    }

    pub fn strand(&self) -> f64 {
        self.strand
    }

    pub fn upstream(&self) -> f64 {
        self.upstream
    }

    pub fn downstream(&self) -> f64 {
        self.downstream
    }

    pub fn asite(&self) -> f64 {
        self.asite
    }

    pub fn incompat(&self) -> f64 {
        self.incompat
    }

    pub fn total(&self) -> f64 {
        self.compat + self.strand + self.upstream + self.downstream + self.asite
    }

    pub fn report(&self) -> String {
        fn report_one(name: &str, num: f64, ttl: f64) -> String {
            format!("{:14} {:12.1}  {:0.3}\n", name, num, num / ttl)
        }
        let mut out = String::new();
        out.push_str(&report_one("compatible:", self.compat(), self.total()));
        // &format!("compatible:   {:12.1}", self.compat()));
        out.push_str(&report_one("wrong strand:", self.strand(), self.total()));
        out.push_str(&report_one("upstream:", self.upstream(), self.total()));
        out.push_str(&report_one("downstream:", self.downstream(), self.total()));
        out.push_str(&report_one("no a site:", self.asite(), self.total()));
        out.push_str(&report_one("incompatible:", self.incompat(), self.total()));
        out.push_str(&report_one("TOTAL", self.total(), self.total()));
        out
    }

    pub fn tally_record(
        &mut self,
        config: &CountConfig,
        trx: &Transcript<Rc<String>>,
        tids: &Tids<Rc<String>>,
        rec: &bam::Record,
    ) -> Result<()> {
        let compat = config.bam_compat(trx, tids, rec)?;
        let weight = config.record_weight(&rec)?;

        match compat {
            Compat::Compat => self.compat += weight,
            Compat::IncompatStrand => self.strand += weight,
            Compat::IncompatUpstream => self.upstream += weight,
            Compat::IncompatDownstream => self.downstream += weight,
            Compat::IncompatASite => self.asite += weight,
            Compat::IncompatNoCDS => self.incompat += weight,
            Compat::IncompatSplice => self.incompat += weight,
        }

        Ok(())
    }
}

impl std::default::Default for Count {
    fn default() -> Self {
        Count {
            compat: 0.0,
            strand: 0.0,
            upstream: 0.0,
            downstream: 0.0,
            asite: 0.0,
            incompat: 0.0,
        }
    }
}

impl std::ops::AddAssign<&Count> for Count {
    fn add_assign(&mut self, rhs: &Count) {
        self.compat += rhs.compat;
        self.strand += rhs.strand;
        self.upstream += rhs.upstream;
        self.downstream += rhs.downstream;
        self.asite += rhs.asite;
        self.incompat += rhs.incompat;
    }
}
