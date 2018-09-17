use std::collections::hash_map::Entry;
use std::collections::HashMap;
use std::fmt;
use std::fs;
use std::str;

use failure;

use bio::io::fastq;

use linkers::*;

#[derive(Debug)]
pub struct Sample {
    name: String,
    dest: fastq::Writer<fs::File>,
    total: usize,
    umi_count: HashMap<Vec<u8>, usize>,
}

impl Sample {
    pub fn new(name: String, dest: fastq::Writer<fs::File>) -> Self {
        Sample {
            name: name,
            dest: dest,
            total: 0,
            umi_count: HashMap::new(),
        }
    }

    pub fn handle_split_read(
        &mut self,
        fq: &fastq::Record,
        split: &LinkerSplit,
    ) -> Result<(), failure::Error> {
        let umi_id = format!("{}#{}", fq.id(), str::from_utf8(split.umi())?);
        let splitfq = fastq::Record::with_attrs(
            umi_id.as_str(),
            fq.desc(),
            split.sequence(),
            split.quality(),
        );

        self.total += 1;
        *self.umi_count.entry(split.umi().to_vec()).or_insert(0) += 1;

        self.dest.write_record(&splitfq)?;
        Ok(())
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn total(&self) -> usize {
        self.total
    }

    pub fn stats_table(&self) -> String {
        let umi_length = self.umi_count.keys().next().map_or(0, |umi| umi.len());
        let mut table = String::new();

        for umi in Self::all_umis(umi_length) {
            table.push_str(&format!(
                "{}\t{}\n",
                str::from_utf8(&umi).unwrap_or("???"),
                self.umi_count.get(&umi).unwrap_or(&0)
            ));
        }

        table
    }

    fn all_umis(len: usize) -> Vec<Vec<u8>> {
        let mut umis = vec![b"".to_vec()];
        for _ in 0..len {
            umis = umis.iter().flat_map(Sample::extend_umi).collect();
        }
        umis
    }

    fn extend_umi(umi: &Vec<u8>) -> Vec<Vec<u8>> {
        [b'A', b'C', b'G', b'T', b'N']
            .into_iter()
            .map(|nt| {
                let mut ext = umi.to_vec();
                ext.push(*nt);
                ext
            })
            .collect()
    }
}

impl fmt::Display for Sample {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.name)
    }
}
