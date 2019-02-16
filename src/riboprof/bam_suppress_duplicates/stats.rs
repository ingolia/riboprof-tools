use std::fmt::Write;

pub struct Stats {
    nlim: usize,
    counts: Vec<u64>, 
    untagged_count: u64
}

impl Stats {
    pub fn new(nlim: usize) -> Self {
        Stats {
            nlim: nlim,
            counts: vec![0; nlim * nlim],
            untagged_count: 0,
        }
    }

    fn index(&self, ntotal: usize, nunique: usize) -> usize {
        (if ntotal >= self.nlim { (self.nlim - 1) } else { ntotal }) * self.nlim
            + (if nunique >= self.nlim { (self.nlim - 1) } else { nunique })
    }
    
    pub fn tally(&mut self, ntotal: usize, nunique: usize) {
        let idx = self.index(ntotal, nunique);
        *self.counts.get_mut(idx).unwrap() += 1;
    }

    pub fn tally_untagged(&mut self) {
        self.untagged_count += 1;
    }

    pub fn dedup_table(&self) -> String {
        let mut table = "ttl\tuniq\tcount\n".to_string();

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
   
