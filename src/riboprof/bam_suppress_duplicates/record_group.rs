use std::cmp::Ordering;

use failure;

use rust_htslib::bam;
use rust_htslib::prelude::*;

/// Groups of records from a sorted BAM file. Record groups must be
/// sorted in ascending order based on the grouping key.
pub struct RecordGroups<'a> {
    bam_reader: &'a mut bam::Reader,
    next_record: Option<bam::Record>,
    group_order: &'a Fn(&bam::Record, &bam::Record) -> Ordering,
}

impl<'a> RecordGroups<'a> {
    /// Create a grouping iterator that uses a provided equivalence /
    /// ordering function to collect individual records into
    /// groups. Records are grouped when they are `Ordering::Equal`
    /// and must be non-decreasing, i.e., the next record after a
    /// group of equivalent records must be `Ordering::Greater`,
    /// relative to that group.
    ///
    /// # Arguments
    ///
    /// * `group_order` provides an ordering for record groups.
    /// 
    /// * `bam_reader` iterates over individual records.
    ///
    /// # Errors
    /// 
    /// An error variant is returned when an error arises reading the
    /// first record from the nested `bam_reader` iterator.
    pub fn new(
        group_order: &'a Fn(&bam::Record, &bam::Record) -> Ordering,
        bam_reader: &'a mut bam::Reader,
    ) -> Result<Self, failure::Error> {
        let mut rg = RecordGroups {
            bam_reader: bam_reader,
            next_record: None,
            group_order: group_order,
        };
        rg.next_record = rg.read_next_record()?;
        Ok(rg)
    }

    /// Creates a record grouping iterator for BAM records sorted by
    /// location according to the ordering provided by `samtools
    /// sort`, as implemented by `RecordGroup::cmp_location`.
    ///
    /// # Arguments
    /// 
    /// * `bam_reader` iterates over individual records.
    ///
    /// # Errors
    ///
    /// An error variant is returned when an error arises reading the
    /// first record from the nested `bam_reader` iterator.
    pub fn new_by_location(bam_reader: &'a mut bam::Reader) -> Result<Self, failure::Error> {
        Self::new(&Self::cmp_location, bam_reader)
    }

    fn read_next_record(&mut self) -> Result<Option<bam::Record>, failure::Error> {
        let mut rec = bam::Record::new();
        match self.bam_reader.read(&mut rec) {
            Ok(()) => Ok(Some(rec)),
            Err(bam::ReadError::NoMoreRecord) => Ok(None),
            Err(e) => Err(e.into()),
        }
    }

    fn read_group(&mut self, curr: bam::Record) -> Result<Vec<bam::Record>, failure::Error> {
        let curr_ref = curr.clone();
        let mut group = Vec::new();
        group.push(curr);

        loop {
            let next = self.read_next_record()?;
            if let Some(rec) = next {
                match (self.group_order)(&curr_ref, &rec) {
                    Ordering::Less => {
                        self.next_record = Some(rec);
                        break;
                    }
                    Ordering::Equal => group.push(rec),
                    Ordering::Greater => {
                        return Err(format_err!(
                            "Records out of order: {:?} > {:?}",
                            curr_ref,
                            rec
                        ));
                    }
                }
            } else {
                self.next_record = None;
                break;
            }
        }

        Ok(group)
    }

    /// Compare BAM records according to location. This ordering
    /// should be consistent with the ordering from `samtools
    /// sort`. Reads are ordered first by reference target sequence
    /// ID, then by position, and then by strand.
    pub fn cmp_location(r1: &bam::Record, r2: &bam::Record) -> Ordering {
        match r1.tid().cmp(&r2.tid()) {
            Ordering::Equal => match r1.pos().cmp(&r2.pos()) {
                Ordering::Equal => r1.is_reverse().cmp(&r2.is_reverse()),
                ordering => ordering,
            },
            ordering => ordering,
        }
    }
}

impl<'a> Iterator for RecordGroups<'a> {
    type Item = Result<Vec<bam::Record>, failure::Error>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(curr) = self.next_record.take() {
            Some(self.read_group(curr))
        } else {
            None
        }
    }
}
