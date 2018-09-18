use std::collections::HashMap;
use std::fmt;
use std::fs;
use std::io;
use std::str;

use failure;

use bio::io::fastq;

use linkers::*;

/// Collected information about one particular sample
pub struct Sample {
    name: String,
    index: Vec<u8>,
    dest: fastq::Writer<Box<io::Write>>,
    total: usize,
    umi_count: HashMap<Vec<u8>, usize>,
}

impl Sample {
    /// Creates new sample information
    ///
    /// # Arguments
    ///
    /// * `name` is the display name for the sample
    ///
    /// * `index` is the sample index sequence
    ///
    /// * `dest` is the output writer for processed fastq records for this sample
    pub fn new<W: io::Write + 'static>(name: String, index: Vec<u8>, dest: W) -> Self {
        Sample {
            name: name,
            index: index,
            dest: fastq::Writer::new(Box::new(dest)),
            total: 0,
            umi_count: HashMap::new(),
        }
    }

    /// Handle a fastq record after linker trimming. This function
    /// will write a new fastq record to the sample output writer,
    /// using the trimmed sequence and quality. The UMI will be
    /// appended to the record `id`, after a `#` character. This
    /// function does not check the sample index in the `LinkerSplit`
    /// result.
    ///
    /// The `Sample` also collects statistics on the total number of
    /// reads, and the number of reads per UMI.
    ///
    /// # Arguments
    ///
    /// * `fq` is an input fastq record
    ///
    /// * `split` contains the results of linker trimming and processing
    ///
    /// # Errors
    ///
    /// An error variant is returned when problems arise in writing
    /// the processed fastq record to the output file.
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

    /// Returns the name of the sample
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Returns the index of the sample
    pub fn index(&self) -> &[u8] {
        &self.index
    }

    /// Returns the total number of reads handled for the sample
    pub fn total(&self) -> usize {
        self.total
    }

    /// Returns a table of the number of reads per UMI
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
