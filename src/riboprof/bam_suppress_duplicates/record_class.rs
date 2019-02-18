use rust_htslib::bam;

/// Group BAM records according to an equivalence function. All
/// records are read and classified, and this function does not
/// require the reads to be sorted in any way. The number of
/// equivalence comparisons, and the minimum running time, scales
/// according to _N*M_, where _N_ is the number of records and _M_ is
/// the number of groups.
pub struct RecordClass<'a> {
    classes: Vec<Vec<bam::Record>>,
    same_class: &'a Fn(&bam::Record, &bam::Record) -> bool,
}

impl<'a> RecordClass<'a> {
    /// Create a new BAM record classifier.
    ///
    /// # Arguments
    ///
    /// * `same_class` specifies an equivalence function over BAM
    /// records for grouping.
    pub fn new(same_class: &'a Fn(&bam::Record, &bam::Record) -> bool) -> Self {
        RecordClass {
            classes: Vec::new(),
            same_class: same_class,
        }
    }

    /// Inserts a BAM record.
    ///
    /// # Arguments
    ///
    /// * `r` is the record to be added
    pub fn insert(&mut self, r: bam::Record) {
        for c in self.classes.iter_mut() {
            if (self.same_class)(c.first().unwrap(), &r) {
                c.push(r);
                return;
            }
        }

        self.classes.push(vec![r]);
    }

    /// Consumes an iterator over BAM records and classifies all of
    /// them.
    ///
    /// # Arguments
    ///
    /// * `iter` yields BAM records to be classified
    pub fn insert_all<I>(&mut self, iter: I)
    where
        I: Iterator<Item = bam::Record>,
    {
        for r in iter {
            self.insert(r);
        }
    }

    /// Record classification. Returns a `Vec` of record classes, each
    /// of which is a non-empty `Vec`.
    pub fn classes(self) -> Vec<Vec<bam::Record>> {
        self.classes
    }
}
