use rust_htslib::bam;

pub struct ReadClass<'a> {
    classes: Vec<Vec<bam::Record>>,
    same_class: &'a Fn(&bam::Record, &bam::Record) -> bool,
}

impl<'a> ReadClass<'a> {
    pub fn new(same_class: &'a Fn(&bam::Record, &bam::Record) -> bool) -> Self {
        ReadClass {
            classes: Vec::new(),
            same_class: same_class,
        }
    }

    pub fn insert(&mut self, r: bam::Record) {
        for c in self.classes.iter_mut() {
            if (self.same_class)(c.first().unwrap(), &r) {
                c.push(r);
                return;
            }
        }

        self.classes.push(vec![r]);
    }

    pub fn insert_all<I>(&mut self, iter: I)
    where
        I: Iterator<Item = bam::Record>,
    {
        for r in iter {
            self.insert(r);
        }
    }

    pub fn classes(self) -> Vec<Vec<bam::Record>> {
        self.classes
    }
}
