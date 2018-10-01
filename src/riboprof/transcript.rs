use std::cmp::{min,max};
use std::collections::HashMap;
use std::error::Error;
use std::fmt; 
use std::hash::Hash;
use std::num::ParseIntError;
use std::ops::{Deref, Range};

use bio_types::annot::loc::Loc;
use bio_types::annot::pos::*;
use bio_types::annot::refids::RefIDSet;
use bio_types::annot::spliced::*;
use bio_types::strand::*;
use bio::data_structures::annot_map::AnnotMap;
use bio::io::bed;

pub struct Transcript<R> {
    gene: R,
    trx: R,
    loc: Spliced<R, ReqStrand>,
    cds: Option<Range<usize>>,
}

impl <R> Transcript<R> {
    pub fn loc(&self) -> &Spliced<R, ReqStrand> { &self.loc }
    pub fn cds_range(&self) -> &Option<Range<usize>> { &self.cds }
}

impl <R> Transcript<R>
    where R: Deref<Target = String>
{
    pub fn gene(&self) -> &str { &self.gene }
    pub fn trx(&self) -> &str { &self.trx }

}

impl <R> Transcript<R>
    where R: Deref<Target = String> + From<String> + Eq + Clone
{
    pub fn from_bed12(record: &bed::Record, refids: &mut RefIDSet<R>) -> Result<Self, TrxError> {
        let loc = Self::loc_from_bed(record, refids)?;
        let cds = Self::cds_from_bed(record, &loc)?;
        let name = record.name().ok_or_else(|| TrxError::bed(record, "No name"))?;

        Ok( Transcript { gene: refids.intern(name), trx: refids.intern(name), loc: loc, cds: cds } )
    }

    const STRAND_COL: usize = 5;
    const THICK_START_COL: usize = 6;
    const THICK_END_COL: usize = 7;
    const BLOCK_COUNT_COL: usize = 9;
    const BLOCK_SIZES_COL: usize = 10;
    const BLOCK_STARTS_COL: usize = 11;

    fn loc_from_bed(record: &bed::Record, refids: &mut RefIDSet<R>) -> Result<Spliced<R, ReqStrand>, TrxError> {
        let block_count = record.aux(Self::BLOCK_COUNT_COL)
            .ok_or_else(|| TrxError::bed(record, "No splicing blocks"))?
            .parse::<usize>()
            .map_err(|err| TrxError::bed_parse(record, "Bad block count", err))?;

        let block_sizes = record.aux(Self::BLOCK_SIZES_COL)
            .ok_or_else(|| TrxError::bed(record, "No splicing block sizes"))?
            .split_terminator(",")
            .map(|size_str| size_str.parse::<usize>())
            .collect::<Result<Vec<usize>, ParseIntError>>()
            .map_err(|err| TrxError::bed_parse(record, "Bad block sizes", err))?;

        let block_starts = record.aux(Self::BLOCK_STARTS_COL)
            .ok_or_else(|| TrxError::bed(record, "No splicing block starts"))?
            .split_terminator(",")
            .map(|size_str| size_str.parse::<usize>())
            .collect::<Result<Vec<usize>, ParseIntError>>()
            .map_err(|err| TrxError::bed_parse(record, "Bad block starts", err))?;

        if block_sizes.len() != block_count || block_starts.len() != block_count {
            return Err(TrxError::bed(record, &format!("block count = {}, |sizes| = {}, |starts| = {}", block_count, block_sizes.len(), block_starts.len())));
        }

        let strand = match record.aux(Self::STRAND_COL)
            .ok_or_else(|| TrxError::bed(record, "No strand"))? {
                "+" => Ok(ReqStrand::Forward),
                "-" => Ok(ReqStrand::Reverse),
                _ => Err(TrxError::bed(record, "Bad strand")),
            }?;
        
        Spliced::with_lengths_starts(refids.intern(record.chrom()), 
                                     record.start() as isize,
                                     &block_sizes,
                                     &block_starts,
                                     strand)
            .map_err(|err| TrxError::BedSplicing(format!("Splicing error on record {:?}", record), err))
    }
    
    fn cds_from_bed(record: &bed::Record, loc: &Spliced<R, ReqStrand>) -> Result<Option<Range<usize>>, TrxError> {
        let thick_start = match record.aux(Self::THICK_START_COL) {
            Some(start_str) => start_str.parse::<usize>()
                .map_err(|err| TrxError::bed_parse(record, "thickStart", err))?,
            None => return Ok(None),
        };

        let thick_end = match record.aux(Self::THICK_END_COL) {
            Some(end_str) => end_str.parse::<usize>()
                .map_err(|err| TrxError::bed_parse(record, "thickEnd", err))?,
            None => return Ok(None),
        };

        if thick_start >= thick_end {
            return Ok(None);
        }

        // thick_end position is outside of loc when CDS extends to
        // the end of the transcript. The last position is guaranteed
        // to be exonic (not intronic), so work with thick_end - 1.

        let left_pos = loc.pos_into(&Pos::new(loc.refid().clone(), thick_start as isize, loc.strand()))
            .ok_or_else(|| TrxError::bed(record, "thickStart not in annot"))?;

        let right_pos = loc.pos_into(&Pos::new(loc.refid().clone(), thick_end as isize - 1, loc.strand()))
            .ok_or_else(|| TrxError::bed(record, "thickEnd-1 not in annot"))?;

        let start = min(left_pos.pos(), right_pos.pos()) as usize;
        let last = max(left_pos.pos(), right_pos.pos()) as usize;

        assert!(last >= start);

        Ok(Some(Range{start: start, end: last+1}))
    }
}

pub struct Transcriptome<R>
    where R: Eq + Hash
{
    gene_to_trxs: HashMap<R,Vec<R>>,
    trx_to_gene: HashMap<R,R>,
    trx_to_transcript: HashMap<R,Transcript<R>>,
    trx_by_location: AnnotMap<R,R>
}

impl <R: Eq + Hash> Transcriptome<R> {
    pub fn new() -> Self {
        Transcriptome { gene_to_trxs: HashMap::new(),
                        trx_to_gene: HashMap::new(),
                        trx_to_transcript: HashMap::new(),
                        trx_by_location: AnnotMap::new() }
    }
}

impl <R> Transcriptome<R>
    where R: Deref<Target = String> + From<String> + Clone + Hash + Eq
{
    pub fn insert(&mut self, transcript: Transcript<R>) -> Result<R, TrxError> {
        if self.trx_to_transcript.contains_key(&transcript.trx) {
            return Err(TrxError::TrxExists(transcript.trx.to_string()));
        }

        let trx = transcript.trx.clone();

        self.trx_to_gene.insert(trx.clone(), transcript.gene.clone());
        self.gene_to_trxs.entry(transcript.gene.clone()).or_insert(vec![]).push(trx.clone());
        self.trx_by_location.insert_at(trx.clone(), &transcript.loc);

        self.trx_to_transcript.insert(trx.clone(), transcript);

        Ok(trx)
    }
}

#[derive(Debug)]
pub enum TrxError {
    Bed(String),
    BedParse(String, ParseIntError),
    BedSplicing(String, SplicingError),
    TrxExists(String),
}

impl TrxError {
    fn bed(record: &bed::Record, message: &str) -> TrxError {
        TrxError::Bed(format!("{} converting BED record {:?}", message, record))
    }

    fn bed_parse(record: &bed::Record, message: &str, parse_error: ParseIntError) -> TrxError {
        TrxError::BedParse(format!("{} parsing BED record {:?}", message, record), parse_error)
    }
}

impl Error for TrxError {

}

impl fmt::Display for TrxError {
    fn fmt(&self, f: &mut fmt::Formatter) -> Result<(), fmt::Error>
    {
        match self {
            TrxError::Bed(msg) => write!(f, "BED record to transcript: {}", msg),
            TrxError::BedParse(msg, err) => write!(f, "BED record to transcript: {}: parsing error {}", msg, err),
            TrxError::BedSplicing(msg, err) => write!(f, "BED record to transcript: {}: splicing error {}", msg, err),
            TrxError::TrxExists(trx) => write!(f, "Transcript already exists: {}", trx),
        }
    }
}
