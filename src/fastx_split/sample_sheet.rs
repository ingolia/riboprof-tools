use std::collections::hash_map::Entry;
use std::collections::HashMap;
use std::error;
use std::fmt;
use std::str;

use failure;

#[derive(Debug, Clone)]
pub struct SampleMap<T> {
    index_length: usize,
    index_map: HashMap<Vec<u8>, usize>,
    sample_index: Vec<Vec<u8>>,
    sample_thing: Vec<T>,
}

impl<T> SampleMap<T> {
    pub fn new(index_length: usize) -> Self {
        SampleMap {
            index_length: index_length,
            index_map: HashMap::new(),
            sample_index: Vec::new(),
            sample_thing: Vec::new(),
        }
    }

    pub fn insert(
        &mut self,
        index: Vec<u8>,
        allow_mismatch: bool,
        thing: T,
    ) -> Result<(), failure::Error> {
        if index.len() != self.index_length {
            return Err(SampleError::IndexBadLength(self.index_length, index).into());
        }

        let i = self.sample_index.len();
        self.sample_index.push(index.to_vec());
        self.sample_thing.push(thing);

        self.insert_index(index.clone(), i)?;

        if allow_mismatch {
            for mm in 0..index.len() {
                for nt in [b'A', b'C', b'G', b'T'].iter() {
                    if index[mm] != *nt {
                        let mut index_mut = index.clone();
                        index_mut[mm] = *nt;
                        self.insert_index(index_mut, i)?;
                    }
                }
            }
        }

        Ok(())
    }

    fn insert_index(&mut self, index: Vec<u8>, entry: usize) -> Result<(), failure::Error> {
        match self.index_map.entry(index) {
            Entry::Occupied(occ) => Err(SampleError::IndexClash(occ.key().to_vec())),
            Entry::Vacant(vac) => Ok(vac.insert(entry)),
        }?;
        Ok(())
    }

    #[allow(dead_code)]
    pub fn get(&self, index: &[u8]) -> Result<Option<&T>, failure::Error> {
        if index.len() != self.index_length {
            return Err(SampleError::IndexBadLength(self.index_length, index.to_vec()).into());
        }

        let entry = self.index_map.get(index);

        Ok(entry.map(|entry| &self.sample_thing[*entry]))
    }

    pub fn get_mut(&mut self, index: &[u8]) -> Result<Option<&mut T>, failure::Error> {
        if index.len() != self.index_length {
            return Err(SampleError::IndexBadLength(self.index_length, index.to_vec()).into());
        }

        let entry = self.index_map.get(index).map(|entry| *entry);

        Ok(entry.map(move |entry| &mut self.sample_thing[entry]))
    }
}

impl<T: fmt::Display> SampleMap<T> {
    pub fn mapping_table(&self) -> String {
        let mut table = String::new();
        for (index, entry) in self.index_map.iter() {
            table.push_str(&format!(
                "{}\t{}\t{}\n",
                str::from_utf8(index).unwrap(),
                self.sample_thing[*entry],
                str::from_utf8(&self.sample_index[*entry]).unwrap()
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
