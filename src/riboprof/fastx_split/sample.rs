use std::collections::HashMap;
use std::fmt;
use std::io;
use std::str;

use failure;

use bio::io::fastq;

use fastx_split::linkers::*;

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

#[cfg(test)]
mod tests {
    use super::*;

    use std::cell::*;
    use std::ops::*;
    use std::rc::*;

    use fastx_split::linkers::*;

    struct TestWriter {
        dest: Rc<RefCell<Vec<u8>>>,
    }

    impl io::Write for TestWriter {
        fn write(&mut self, buf: &[u8]) -> Result<usize, io::Error> {
            self.dest.borrow_mut().append(&mut buf.to_vec());
            Ok(buf.len())
        }

        fn flush(&mut self) -> Result<(), io::Error> {
            Ok(())
        }
    }

    #[test]
    fn sample_output() {
        let outbuf = Rc::new(RefCell::new(Vec::new()));

        {
            let writer = TestWriter {
                dest: outbuf.clone(),
            };
            let mut sample = Sample::new("One".to_string(), b"ACGT".to_vec(), writer);

            let linker_spec = LinkerSpec::new("NN", "NNIIII").unwrap();

            let rec1 =
                fastq::Record::with_attrs("test_record", None, b"ACGTACGTACGTACGT", &vec![40; 16]);
            let spl1 = linker_spec.split_record(&rec1).unwrap();
            sample.handle_split_read(&rec1, &spl1).unwrap();
            assert!(sample.total() == 1);

            let rec2 =
                fastq::Record::with_attrs("another", None, b"TGTGCGAGCTAGTCACTC", &vec![37; 18]);
            let spl2 = linker_spec.split_record(&rec2).unwrap();
            sample.handle_split_read(&rec2, &spl2).unwrap();
            assert!(sample.total() == 2);
        }

        let mut exp = b"@test_record#ACGT\nGTACGTAC\n+\n((((((((\n".to_vec();
        let mut exp2 = b"@another#TGTC\nTGCGAGCTAG\n+\n%%%%%%%%%%\n".to_vec();
        exp.append(&mut exp2);

        assert!(outbuf.borrow().as_slice() == exp.as_slice());
    }

    #[test]
    fn sample_umi_counts() {
        let linker_spec = LinkerSpec::new("", "NN").unwrap();

        let mut sample = Sample::new("Two".to_string(), Vec::new(), io::sink());

        for nt1 in b"ACCGGGTTTT" {
            for nt2 in b"AAAACCCGGT" {
                let mut seq = b"TGGTGCCGCAAC".to_vec();
                seq.push(*nt1);
                seq.push(*nt2);
                let rec = fastq::Record::with_attrs("test", None, &seq, &vec![40; seq.len()]);
                let spl = linker_spec.split_record(&rec).unwrap();
                sample.handle_split_read(&rec, &spl).unwrap();
            }
        }

        let mut exp = "AA\t4\nAC\t3\nAG\t2\nAT\t1\nAN\t0\n".to_string();
        exp.push_str("CA\t8\nCC\t6\nCG\t4\nCT\t2\nCN\t0\n");
        exp.push_str("GA\t12\nGC\t9\nGG\t6\nGT\t3\nGN\t0\n");
        exp.push_str("TA\t16\nTC\t12\nTG\t8\nTT\t4\nTN\t0\n");
        exp.push_str("NA\t0\nNC\t0\nNG\t0\nNT\t0\nNN\t0\n");

        assert!(sample.stats_table() == exp);
    }
}
