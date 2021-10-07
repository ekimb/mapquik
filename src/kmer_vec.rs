use std::hash::{Hash, Hasher};
use std::vec::Vec;
use std::cmp::Ordering;

#[derive(Clone, Debug)]
pub struct KmerVec {
    pub hashes: Vec<u64>,
    pub start: usize,
    pub end: usize,
    pub offset: usize
}

pub fn get(a: &KmerVec) -> Vec<u64> {
    a.hashes.to_vec()
}

impl KmerVec {
    pub fn suffix(&self) -> KmerVec {
        let mut res = KmerVec {hashes: self.hashes.clone(), start: self.start, end: self.end, offset: self.offset};
        res.hashes.remove(0);
        res
    }
    
    pub fn prefix(&self) -> KmerVec {
        let mut res = KmerVec {hashes: self.hashes.clone(), start: self.start, end: self.end, offset: self.offset};
        res.hashes.pop();
        res
    }
    
    pub fn reverse(&self) -> KmerVec {
        let mut res = KmerVec {hashes: self.hashes.clone(), start: self.start, end: self.end, offset: self.offset};
        res.hashes.reverse();
        res
    }

    pub fn normalize(&self) -> (KmerVec,bool) {
        let rev = self.reverse();
        if *self < rev {(self.clone(), false)}
        else {(rev, true)}
    }

    pub fn make_from(ar: &[u64], s: usize, e: usize, o: usize) -> KmerVec {
        KmerVec{hashes: Vec::from(ar), start: s, end: e, offset: o}
    }

    pub fn print_as_string(&self) -> String {
        format!("{:?}", &self.hashes)
    }

    pub fn minimizers(&self) -> &Vec<u64> {
        &self.hashes
    }
}

impl PartialEq for KmerVec {
    fn eq(&self, other: &KmerVec) -> bool {
        self.hashes == other.hashes
    }
}

impl Eq for KmerVec {
}

impl Hash for KmerVec {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.hashes.hash(state);
    }
}

impl Default for KmerVec {
    fn default() -> Self{KmerVec{hashes: vec![], start: 0, end: 0, offset: 0}}
}

impl Ord for KmerVec {
    fn cmp(&self, other: &Self) -> Ordering {
        self.hashes.cmp(&other.hashes)
    }
}


impl PartialOrd for KmerVec{
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
