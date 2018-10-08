use std::ops::Range;
use std::rc::Rc;

use rust_htslib::bam;

use fp_framing::stats::*;

use transcript::*;

pub fn record_framing(trxome: &Transcriptome<Rc<String>>, 
                      rec: &bam::Record,
                      framing_stats: &mut FramingStats,
                      cdsbody: &Range<isize>,
                      count_multi: bool)
                      -> String
{
    "N/A".to_string()
}


