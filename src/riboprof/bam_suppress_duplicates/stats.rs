use std::fmt::Write;

pub struct Stats {
    nlim: usize,
    counts: Vec<u64>,

    total_reads_count: u64,
    unique_reads_count: u64,
    total_sites_count: u64,
    dupl_sites_count: u64,
    untagged_count: u64,
}

impl Stats {
    pub fn new(nlim: usize) -> Self {
        Stats {
            nlim: nlim,
            counts: vec![0; nlim * nlim],
            total_reads_count: 0,
            unique_reads_count: 0,
            total_sites_count: 0,
            dupl_sites_count: 0,
            untagged_count: 0,
        }
    }

    fn index(&self, ntotal: usize, nunique: usize) -> usize {
        (if ntotal >= self.nlim {
            (self.nlim - 1)
        } else {
            ntotal
        }) * self.nlim + (if nunique >= self.nlim {
            (self.nlim - 1)
        } else {
            nunique
        })
    }

    pub fn untagged_reads(&self) -> u64 {
        self.untagged_count
    }
    pub fn total_reads(&self) -> u64 {
        self.total_reads_count
    }
    pub fn unique_reads(&self) -> u64 {
        self.unique_reads_count
    }
    pub fn dupl_reads(&self) -> u64 {
        self.total_reads_count - self.unique_reads_count
    }
    pub fn total_sites(&self) -> u64 {
        self.total_sites_count
    }
    pub fn dupl_sites(&self) -> u64 {
        self.dupl_sites_count
    }

    pub fn tally(&mut self, ntotal: usize, nunique: usize) {
        let idx = self.index(ntotal, nunique);
        *self.counts.get_mut(idx).unwrap() += 1;
        self.total_reads_count += ntotal as u64;
        self.unique_reads_count += nunique as u64;
        self.total_sites_count += 1;
        self.dupl_sites_count += if ntotal > nunique { 1 } else { 0 };
    }

    pub fn tally_untagged(&mut self) {
        self.untagged_count += 1;
    }

    pub fn dedup_table(&self) -> String {
        let mut table = "ttl\tuniq\tcount\n".to_string();

        if self.untagged_count > 0 {
            write!(table, "0\t0\t{}\n", self.untagged_count).unwrap();
        }

        for ttl in 0..(self.nlim - 1) {
            for uniq in 0..(self.nlim - 1) {
                let ct = self.counts[self.index(ttl, uniq)];
                if ct > 0 {
                    write!(table, "{}\t{}\t{}\n", ttl, uniq, ct).unwrap();
                }
            }
        }

        table
    }
}
