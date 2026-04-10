use anyhow::{Result, anyhow, bail};
use bio_types::annot::loc::Loc;
use bio_types::annot::pos::Pos;
use bio_types::strand::ReqStrand;
use rust_htslib::bam;

pub struct GenomeCounts<N> {
    chr_counts: Vec<ChrCounts<N>>,
}

impl<N> GenomeCounts<N> {
    pub fn new(hdr: &bam::HeaderView) -> Result<Self> {
        let mut chr_counts = Vec::with_capacity(hdr.target_count() as usize);
        for tid in 0..hdr.target_count() {
            chr_counts.push(ChrCounts::new(
                String::from_utf8(hdr.tid2name(tid as u32).to_vec())?,
                hdr.target_len(tid as u32)
                    .ok_or_else(|| anyhow!("No target len for tid {}", tid))?,
            ));
        }
        Ok(GenomeCounts { chr_counts })
    }

    pub fn chroms(&self) -> &[ChrCounts<N>] {
        &self.chr_counts
    }
}

impl<N: Default> GenomeCounts<N> {
    pub fn update<F>(&mut self, pos: Pos<u32, ReqStrand>, upd: F) -> Result<()>
    where
        F: FnOnce(&mut N),
    {
        self.chr_counts
            .get_mut(*pos.refid() as usize)
            .ok_or_else(|| anyhow!("Bad target id {}", pos.refid()))?
            .update(pos.pos() as usize, pos.strand(), upd)
    }
}

pub struct ChrCounts<N> {
    name: String,
    counts_fwd: WindowCounts<N>,
    counts_rev: WindowCounts<N>,
}

impl<N> ChrCounts<N> {
    pub fn new(name: String, length: u64) -> Self {
        ChrCounts {
            name,
            counts_fwd: WindowCounts::new(length as usize),
            counts_rev: WindowCounts::new(length as usize),
        }
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn fwd(&self) -> &WindowCounts<N> {
        &self.counts_fwd
    }

    pub fn rev(&self) -> &WindowCounts<N> {
        &self.counts_rev
    }
}

impl<N: Default> ChrCounts<N> {
    pub fn update<F>(&mut self, pos: usize, strand: ReqStrand, upd: F) -> Result<()>
    where
        F: FnOnce(&mut N),
    {
        match strand {
            ReqStrand::Forward => &mut self.counts_fwd,
            ReqStrand::Reverse => &mut self.counts_rev,
        }
        .update(pos, upd)
    }
}

pub struct WindowCounts<N> {
    counts: Vec<Option<Vec<N>>>,
    length: usize,
}

impl<N> WindowCounts<N> {
    const WINDOW_BITS: usize = 12;

    pub fn new(length: usize) -> Self {
        let n_windows = 1 + (length >> Self::WINDOW_BITS);
        let mut counts = Vec::with_capacity(n_windows);
        for _ in 0..n_windows {
            counts.push(None)
        }
        WindowCounts { counts, length }
    }

    fn window_len(&self) -> usize {
        1 << Self::WINDOW_BITS
    }

    pub fn sparse_pos_iter(&self) -> impl Iterator<Item = (usize, &N)> {
        self.counts
            .iter()
            .enumerate()
            .filter_map(|(widx, opt_window)| {
                opt_window.as_ref().map(|window| {
                    let wbase = widx * (self.window_len());
                    window
                        .iter()
                        .enumerate()
                        .map(move |(wpos, obj)| (wpos + wbase, obj))
                })
            })
            .flatten()
    }
}

impl<N: Default> WindowCounts<N> {
    pub fn update<F>(&mut self, pos: usize, upd: F) -> Result<()>
    where
        F: FnOnce(&mut N),
    {
        if pos >= self.length {
            bail!("pos {} > length {}", pos, self.length);
        }

        let widx = pos >> Self::WINDOW_BITS;
        let woff = pos - (widx << Self::WINDOW_BITS);

        if self.counts[widx].is_none() {
            let mut v = Vec::with_capacity(self.window_len());
            v.resize_with(self.window_len(), Default::default);
            self.counts[widx] = Some(v);
        }

        if let Some(window) = self.counts[widx].as_mut() {
            upd(&mut window[woff]);
        }

        Ok(())
    }
}

pub fn tally(ctptr: &mut usize) -> () {
    *ctptr += 1;
}
