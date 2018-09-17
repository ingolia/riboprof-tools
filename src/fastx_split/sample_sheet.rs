use std::cell::*;
use std::collections::hash_map::Entry;
use std::collections::HashMap;
use std::error;
use std::fmt;
use std::rc::Rc;
use std::str;

use failure;

#[derive(Debug, Clone)]
pub struct SampleMap<T> {
    index_length: usize,
    index_map: HashMap<Vec<u8>, SampleEntry<T>>,
    unknown: SampleEntry<T>,
    entries: Vec<SampleEntry<T>>,
}

#[derive(Debug, Clone)]
struct SampleEntry<T> {
    true_index: Vec<u8>,
    thing: Rc<RefCell<T>>,
}

impl<T> SampleEntry<T> {
    pub fn new<I: AsRef<[u8]>>(true_index: I, thing: &Rc<RefCell<T>>) -> Self {
        SampleEntry {
            true_index: true_index.as_ref().to_vec(),
            thing: thing.clone(),
        }
    }
}

impl<T: fmt::Display> fmt::Display for SampleEntry<T> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.thing.borrow())
    }
}

impl<T> SampleMap<T> {
    pub fn new(index_length: usize, unknown: T) -> Self {
        let unknown_rcrc = Rc::new(RefCell::new(unknown));
        let unknown_index = vec![b'N'; index_length];

        SampleMap {
            index_length: index_length,
            index_map: HashMap::new(),
            unknown: SampleEntry::new(&unknown_index, &unknown_rcrc),
            entries: vec![SampleEntry::new(&unknown_index, &unknown_rcrc)],
        }
    }

    pub fn insert(
        &mut self,
        index: Vec<u8>,
        allow_mismatch: bool,
        thing: T,
    ) -> Result<Rc<RefCell<T>>, failure::Error> {
        if index.len() != self.index_length {
            return Err(SampleError::IndexBadLength(self.index_length, index).into());
        }

        let rcrc = Rc::new(RefCell::new(thing));

        self.insert_index(index.clone(), SampleEntry::new(index.as_slice(), &rcrc))?;

        if allow_mismatch {
            for mm in 0..index.len() {
                for nt in [b'A', b'C', b'G', b'T'].iter() {
                    if index[mm] != *nt {
                        let mut index_mut = index.clone();
                        index_mut[mm] = *nt;
                        self.insert_index(index_mut, SampleEntry::new(index.as_slice(), &rcrc))?;
                    }
                }
            }
        }

        self.entries.push(SampleEntry::new(index.as_slice(), &rcrc));

        Ok(rcrc)
    }

    fn insert_index(
        &mut self,
        index: Vec<u8>,
        entry: SampleEntry<T>,
    ) -> Result<(), failure::Error> {
        match self.index_map.entry(index) {
            Entry::Occupied(occ) => Err(SampleError::IndexClash(occ.key().to_vec())),
            Entry::Vacant(vac) => Ok(vac.insert(entry)),
        }?;
        Ok(())
    }

    #[allow(dead_code)]
    pub fn get(&self, index: &[u8]) -> Result<Ref<T>, failure::Error> {
        if index.len() != self.index_length {
            return Err(SampleError::IndexBadLength(self.index_length, index.to_vec()).into());
        }

        let entry = self.index_map.get(index).unwrap_or(&self.unknown);
        let thing = entry.thing.try_borrow()?;
        Ok(thing)
    }

    pub fn get_mut(&mut self, index: &[u8]) -> Result<RefMut<T>, failure::Error> {
        if index.len() != self.index_length {
            return Err(SampleError::IndexBadLength(self.index_length, index.to_vec()).into());
        }

        let entry = self.index_map.get(index).unwrap_or(&self.unknown);
        let thing = entry.thing.try_borrow_mut()?;
        Ok(thing)
    }

    pub fn things(&self) -> Vec<Rc<RefCell<T>>> {
        let mut things = Vec::new();
        for entry in self.entries.iter() {
            things.push(entry.thing.clone());
        }
        things
    }
}

impl<T: fmt::Display> SampleMap<T> {
    pub fn mapping_table(&self) -> String {
        let mut table = String::new();
        for (index, entry) in self.index_map.iter() {
            table.push_str(&format!(
                "{}\t{}\t{}\n",
                str::from_utf8(index).unwrap(),
                entry,
                str::from_utf8(entry.true_index.as_slice()).unwrap()
            ));
        }
        table
    }
}

pub fn parse_sample_sheet(sheet: &str) -> Result<Vec<(String, String)>, failure::Error> {
    sheet.lines().map(parse_sample_line).collect()
}

fn parse_sample_line(line: &str) -> Result<(String, String), failure::Error> {
    let mut field_iter = line.split(',');
    let name = field_iter
        .next()
        .ok_or_else(|| SampleError::BadSheetLine(line.to_string()))?;
    let idx = field_iter
        .next()
        .ok_or_else(|| SampleError::BadSheetLine(line.to_string()))?;
    Ok((name.to_string(), idx.to_string()))
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum SampleError {
    BadSheetLine(String),
    IndexBadLength(usize, Vec<u8>),
    IndexClash(Vec<u8>),
}

impl fmt::Display for SampleError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            SampleError::BadSheetLine(line) => write!(f, "Bad sample sheet line: \"{}\"", line),
            SampleError::IndexBadLength(ilen, idx) => write!(
                f,
                "Index length wrong: index \"{}\" but length {}",
                str::from_utf8(idx).unwrap_or("???"),
                ilen
            ),
            SampleError::IndexClash(idx) => write!(
                f,
                "Index clash: index {}",
                str::from_utf8(idx).unwrap_or("???")
            ),
        }
    }
}

impl error::Error for SampleError {}
