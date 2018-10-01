use std::default::Default;
//use std::error::Error;
//use std::fmt;
//use std::iter;

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
}

// impl <T> Index<isize> for Metagene<T> {
//     type Output = Option<T>;

//     fn index(&self, pos: isize) -> &Option<T> {

//     }
// }