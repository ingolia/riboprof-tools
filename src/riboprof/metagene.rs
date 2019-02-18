use std::default::Default;
//use std::error::Error;
//use std::fmt;
use std::iter;
use std::slice;
use std::vec;

//use failure;

#[derive(Clone, Debug)]
pub struct LenProfile<T> {
    short: T,
    len_vec: Vec<T>,
    long: T,
    minlen: usize,
}

impl<T: Clone> LenProfile<T> {
    pub fn new(minlen: usize, maxlen: usize, initial: T) -> Self {
        let nlen = if minlen > maxlen {
            0
        } else {
            1 + maxlen - minlen
        };

        LenProfile {
            short: initial.clone(),
            len_vec: vec![initial.clone(); nlen],
            long: initial.clone(),
            minlen: minlen,
        }
    }
}

impl<T: Default> LenProfile<T> {
    pub fn new_with_default(minlen: usize, maxlen: usize) -> Self {
        let nlen = if minlen > maxlen {
            0
        } else {
            1 + maxlen - minlen
        };

        let mut len_vec = Vec::new();
        for _i in 0..nlen {
            len_vec.push(Default::default());
        }

        LenProfile {
            short: Default::default(),
            len_vec: len_vec,
            long: Default::default(),
            minlen: minlen,
        }
    }
}

impl<T> LenProfile<T> {
    pub fn get(&self, len: usize) -> &T {
        if len < self.minlen {
            &self.short
        } else {
            self.len_vec.get(len - self.minlen).unwrap_or(&self.long)
        }
    }

    pub fn get_mut(&mut self, len: usize) -> &mut T {
        if len < self.minlen {
            &mut self.short
        } else {
            self.len_vec
                .get_mut(len - self.minlen)
                .unwrap_or(&mut self.long)
        }
    }

    pub fn named_iter(&self) -> impl Iterator<Item = (String, &T)> {
        let minlen = self.minlen;

        iter::once((format!("<{}", self.minlen), &self.short))
            .chain(
                self.len_vec
                    .iter()
                    .enumerate()
                    .map(move |(i, x)| (format!("{}", i + minlen), x)),
            )
            .chain(iter::once((
                format!("≥{}", self.minlen + self.len_vec.len()),
                &self.long,
            )))
    }

    pub fn iter(&self) -> impl Iterator<Item = &T> {
        self.into_iter()
    }
}

impl<'a, T> IntoIterator for &'a LenProfile<T> {
    type Item = &'a T;
    type IntoIter =
        iter::Chain<iter::Chain<iter::Once<&'a T>, slice::Iter<'a, T>>, iter::Once<&'a T>>;

    fn into_iter(self) -> Self::IntoIter {
        iter::once(&self.short)
            .chain(self.len_vec.iter())
            .chain(iter::once(&self.long))
    }
}

impl<T> IntoIterator for LenProfile<T> {
    type Item = T;
    type IntoIter = iter::Chain<iter::Chain<iter::Once<T>, vec::IntoIter<T>>, iter::Once<T>>;

    fn into_iter(self) -> Self::IntoIter {
        iter::once(self.short)
            .chain(self.len_vec.into_iter())
            .chain(iter::once(self.long))
    }
}

#[derive(Clone, Debug)]
pub struct Frame<T> {
    frames: Vec<T>,
}

impl<T: Clone> Frame<T> {
    pub fn new(initial: T) -> Self {
        Frame {
            frames: vec![initial.clone(), initial.clone(), initial],
        }
    }
}

impl<T: Default> Frame<T> {
    pub fn new_with_default() -> Self {
        Frame {
            frames: vec![Default::default(), Default::default(), Default::default()],
        }
    }
}

impl<T> Frame<T> {
    pub fn to_frame<I>(i: I) -> usize
    where
        I: Into<isize>,
    {
        (((i.into() % 3) + 3) % 3) as usize
    }

    pub fn get<I: Into<isize>>(&self, pos: I) -> &T {
        &self.frames[Self::to_frame(pos)]
    }

    pub fn get_mut<I: Into<isize>>(&mut self, pos: I) -> &mut T {
        &mut self.frames[Self::to_frame(pos)]
    }

    pub fn frame_iter(&self) -> impl Iterator<Item = (isize, &T)> {
        self.frames.iter().enumerate().map(|(i, x)| (i as isize, x))
    }

    pub fn iter(&self) -> impl Iterator<Item = &T> {
        self.into_iter()
    }
}

impl<'a, T> IntoIterator for &'a Frame<T> {
    type Item = &'a T;
    type IntoIter = slice::Iter<'a, T>;

    fn into_iter(self) -> Self::IntoIter {
        self.frames.iter()
    }
}

impl<T> IntoIterator for Frame<T> {
    type Item = T;
    type IntoIter = vec::IntoIter<T>;

    fn into_iter(self) -> Self::IntoIter {
        self.frames.into_iter()
    }
}

#[derive(Clone, Debug)]
pub struct Metagene<T> {
    pos_vec: Vec<T>,
    start: isize,
}

impl<T: Clone> Metagene<T> {
    pub fn new(start: isize, len: usize, initial: T) -> Self {
        Metagene {
            pos_vec: vec![initial; len],
            start: start,
        }
    }
}

impl<T: Default> Metagene<T> {
    pub fn new_with_default(start: isize, len: usize) -> Self {
        let mut pos_vec = Vec::new();
        for _ in 0..len {
            pos_vec.push(Default::default());
        }
        Metagene {
            pos_vec: pos_vec,
            start: start,
        }
    }
}

impl<T> Metagene<T> {
    pub fn get(&self, pos: isize) -> Option<&T> {
        if pos < self.start {
            None
        } else {
            self.pos_vec.get((pos - self.start) as usize)
        }
    }

    pub fn get_mut(&mut self, pos: isize) -> Option<&mut T> {
        if pos < self.start {
            None
        } else {
            self.pos_vec.get_mut((pos - self.start) as usize)
        }
    }

    pub fn pos_iter(&self) -> impl Iterator<Item = (isize, &T)> {
        self.pos_vec
            .iter()
            .enumerate()
            .map(move |(i, x)| (i as isize + self.start, x))
    }

    pub fn iter(&self) -> impl Iterator<Item = &T> {
        self.into_iter()
    }
}

impl<'a, T> IntoIterator for &'a Metagene<T> {
    type Item = &'a T;
    type IntoIter = slice::Iter<'a, T>;

    fn into_iter(self) -> Self::IntoIter {
        self.pos_vec.iter()
    }
}

impl<T> IntoIterator for Metagene<T> {
    type Item = T;
    type IntoIter = vec::IntoIter<T>;

    fn into_iter(self) -> Self::IntoIter {
        self.pos_vec.into_iter()
    }
}
