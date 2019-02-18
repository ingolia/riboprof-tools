use std::cmp::Ordering;

use failure;

use rust_htslib::bam;
use rust_htslib::prelude::*;

pub struct ReadGroups<'a> {
    bam_reader: &'a mut bam::Reader,
    next_record: Option<bam::Record>,
    group_order: &'a Fn(&bam::Record, &bam::Record) -> Ordering,
}

impl<'a> ReadGroups<'a> {
    pub fn new(
        group_order: &'a Fn(&bam::Record, &bam::Record) -> Ordering,
        bam_reader: &'a mut bam::Reader,
    ) -> Result<Self, failure::Error> {
        let mut rg = ReadGroups {
            bam_reader: bam_reader,
            next_record: None,
            group_order: group_order,
        };
        rg.next_record = rg.read_next_record()?;
        Ok(rg)
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
}

impl<'a> Iterator for ReadGroups<'a> {
    type Item = Result<Vec<bam::Record>, failure::Error>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(curr) = self.next_record.take() {
            Some(self.read_group(curr))
        } else {
            None
        }
    }
}
