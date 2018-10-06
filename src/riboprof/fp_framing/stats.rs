pub struct AnnotStats {
    no_gene: usize,
    noncoding: usize,
    noncoding_overlap: usize,
    multi_coding: usize,
    incompatible: usize,
    ambig: usize,
    good: usize
}

impl AnnotStats {
    pub fn new() -> Self {
        AnnotStats { no_gene: 0, noncoding: 0, noncoding_overlap: 0, multi_coding: 0, incompatible: 0, ambig: 0, good: 0 }
    }

    pub fn no_gene(&self) -> usize { self.no_gene }
    pub fn noncoding(&self) -> usize { self.noncoding }
    pub fn noncoding_overlap(&self) -> usize { self.noncoding_overlap }
    pub fn multi_coding(&self) -> usize { self.multi_coding }
    pub fn incompatible(&self) -> usize { self.incompatible }
    pub fn ambig(&self) -> usize { self.ambig }
    pub fn good(&self) -> usize { self.good }

    pub fn tally_no_gene(&mut self) { self.no_gene += 1 }
    pub fn tally_noncoding(&mut self) { self.noncoding += 1 }
    pub fn tally_noncoding_overlap(&mut self) { self.noncoding_overlap += 1 }
    pub fn tally_multi_coding(&mut self) { self.multi_coding += 1 }
    pub fn tally_incompatible(&mut self) { self.incompatible += 1 }
    pub fn tally_ambig(&mut self) { self.ambig += 1 }
    pub fn tally_good(&mut self) { self.good += 1 }

    pub fn total(&self) -> usize {
        self.no_gene + self.noncoding + self.noncoding_overlap + self.multi_coding + self.incompatible + self.ambig + self.good
    }
}

pub struct AlignStats {
    unmapped: usize,
    short: usize,
    long: usize,
    multi_hit: usize,
    annot_stats: AnnotStats
}

impl AlignStats {
    pub fn new() -> Self {
        AlignStats { unmapped: 0, short: 0, long: 0, multi_hit: 0, annot_stats: AnnotStats::new() }
    }

    pub fn unmapped(&self) -> usize { self.unmapped }
    pub fn short(&self) -> usize { self.short }
    pub fn long(&self) -> usize { self.long }
    pub fn multi_hit(&self) -> usize { self.multi_hit }

    pub fn annot_stats(&self) -> &AnnotStats { &self.annot_stats }

    pub fn tally_unmapped(&mut self)  { self.unmapped += 1 }
    pub fn tally_short(&mut self)  { self.short += 1 }
    pub fn tally_long(&mut self)  { self.long += 1 }
    pub fn tally_multi_hit(&mut self)  { self.multi_hit += 1 }

    pub fn annot_stats_mut(&mut self) -> &mut AnnotStats { &mut self.annot_stats }

    pub fn total(&self) -> usize {
        self.unmapped + self.short + self.long + self.multi_hit + self.annot_stats.total()
    }
}
