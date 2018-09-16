use std::collections::HashMap;
use std::collections::hash_map::Entry;
use std::error;
use std::fmt;
use std::str;

use failure;

pub struct SampleMap<T> {
    index_length: usize,
    index_map: HashMap<Vec<u8>, usize>,
    sample_index: Vec<Vec<u8>>,
    sample_name: Vec<String>,
    sample_thing: Vec<T>,
}

impl <T> SampleMap<T> {
    pub fn new(index_length: usize) -> Self {
        SampleMap { index_length: index_length,
                    index_map: HashMap::new(),
                    sample_index: Vec::new(),
                    sample_name: Vec::new(),
                    sample_thing: Vec::new()
        }
    }

    pub fn insert(&mut self, index: Vec<u8>, allow_mismatch: bool, name: String, thing: T) -> Result<(), failure::Error> {
        if index.len() != self.index_length {
            return Err(SampleError::IndexBadLength(self.index_length, index).into());
        }

        let i = self.sample_index.len();
        self.sample_index.push(index.to_vec());
        self.sample_name.push(name);
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

        Ok( () )
    }

    fn insert_index(&mut self, index: Vec<u8>, entry: usize) -> Result<(), failure::Error> {
        match self.index_map.entry(index) {
            Entry::Occupied(occ) => {
                Err(SampleError::IndexClash(occ.key().to_vec(), 
                                            self.sample_name.get(*occ.get()).map_or("???", String::as_str).to_string(),
                                            self.sample_name.get(entry).map_or("???", String::as_str).to_string()))
            }

            Entry::Vacant(vac) => Ok( vac.insert(entry) ),
        }?;
        Ok( () )
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
    IndexClash(Vec<u8>, String, String),
}

impl fmt::Display for SampleError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            SampleError::BadSheetLine(line) => write!(f, "Bad sample sheet line: \"{}\"", line),
            SampleError::IndexBadLength(ilen, idx) => write!(f, "Index length wrong: index {} but length {}", 
                                                             str::from_utf8(idx).unwrap_or("???"), ilen),
            SampleError::IndexClash(idx, sn1, sn2) => write!(f, "Index clash: index {} for samples \"{}\" and \"{}\"",
                                                             str::from_utf8(idx).unwrap_or("???"), sn1, sn2),
            _ => Ok( () ),
        }
    }
}

impl error::Error for SampleError { }
