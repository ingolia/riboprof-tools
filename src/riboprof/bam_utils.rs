use std::ops::Deref;

use failure;

use bio_types::annot::refids::RefIDSet;
use bio_types::annot::spliced::Spliced;
use bio_types::strand::ReqStrand;
use rust_htslib::bam::{HeaderView,record::Cigar,record::CigarStringView};
use rust_htslib::bam;

pub struct Tids<R> {
    tids: Vec<R>
}

impl <R: Deref<Target = String> + From<String> + Clone> Tids<R> {
    pub fn new(refids: &mut RefIDSet<R>, header: &HeaderView) -> Self {
        let mut tids = Vec::with_capacity(header.target_count() as usize);

        for (tid, target_name) in header.target_names().into_iter().enumerate() {
            let target_string = String::from_utf8_lossy(target_name);
            let target_rc = refids.intern(&target_string);
            assert!(tids.len() == tid);
            tids.push(target_rc);
        }

        Tids { tids: tids }
    }
}

impl <R> Tids<R> {
    pub fn get(&self, tid: u32) -> Option<&R> {
        self.tids.get(tid as usize)
    }
}

pub fn bam_to_spliced<R>(tids: &Tids<R>, record: &bam::Record) -> Result<Option<Spliced<R, ReqStrand>>, failure::Error>
    where R: Clone
{
    if record.tid() < 0 {
        return Ok(None);
    }

    let (lengths, starts) = cigar_to_lengths_starts(&record.cigar());

    let refid = tids.get(record.tid() as u32).ok_or_else(|| failure::err_msg(format!("BAM target ID {} out of range", record.tid())))?;

    let strand = if record.is_reverse() { ReqStrand::Reverse } else { ReqStrand::Forward };

    let spliced = Spliced::with_lengths_starts(refid.clone(), record.pos() as isize, lengths.as_slice(), starts.as_slice(), strand)?;

    Ok(Some(spliced))
}

pub fn cigar_to_lengths_starts(cigar_string: &CigarStringView) -> (Vec<usize>, Vec<usize>) {
    let mut starts = Vec::new();
    let mut lengths = Vec::new();

    let mut curr_start = 0;
    let mut curr_end = 0;

    for cigar in cigar_string.iter() {
        match cigar {
            Cigar::Match(len) => curr_end += len,
            Cigar::Ins(_) => (),
            Cigar::Del(len) => curr_end += len,
            Cigar::RefSkip(len) => {
                if curr_end > curr_start {
                    starts.push(curr_start as usize);
                    lengths.push((curr_end - curr_start) as usize);
                }
                curr_start = curr_end + len;
                curr_end = curr_start;
            },
            Cigar::SoftClip(_) => (),
            Cigar::HardClip(_) => (),
            Cigar::Pad(_) => (),
            Cigar::Equal(len) => curr_end += len,
            Cigar::Diff(len) => curr_end += len,
        };
    }

    if curr_end > curr_start {
        starts.push(curr_start as usize);
        lengths.push((curr_end - curr_start) as usize);
    }

    (lengths, starts)
}
