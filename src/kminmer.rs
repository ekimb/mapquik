// kminmer.rs
// Contains the "Kminmer" struct, a k-mer in minimizer-space.

use crate::{utils::pretty_minvec, Params};
use std::cmp::Ordering;
use std::collections::{HashSet, hash_map::DefaultHasher};
use std::hash::{Hash, Hasher};
use std::vec::Vec;
use array_tool::vec::{Intersect, Union};

#[derive(Clone, Debug)]
pub struct Kminmer {
    mers: Vec<u64>, // Raw Vec of minimizer hashes
    pub start: usize, // Start location
    pub end: usize, // End location
    pub offset: usize, // Offset (index in the k-min-mer array)
    pub rev: bool, // Strand direction
}

impl Kminmer {

    // Create a new Kminmer object.
    pub fn new(mers: &[u64], start: usize, end: usize, offset: usize) -> Self {
        let mut obj = Kminmer {
            mers: Vec::from(mers),
            start,
            end,
            offset,
            rev: false,
        };
        obj.normalize();
        obj     
    }

    // Obtain the canonical Kminmer for this object.
    pub fn normalize(&mut self) {
        let mut rev_mers = self.mers.clone();
        rev_mers.reverse();
        if rev_mers < self.mers {
            self.mers = rev_mers.to_vec();
            self.rev = true;
        }
    }

    // Pretty-print a Kminmer.
    pub fn print(&self) -> String {
        pretty_minvec(&self.mers)
    }

    // Obtain a raw Vec of minimizer hashes.
    pub fn mers(&self) -> Vec<u64> {
        self.mers.to_vec()
    }

    // Hash the Vec of minimizer hashes to a u64 (this is used throughout the reference processing).
    pub fn get_hash(&self) -> u64 {
        let mut hash = DefaultHasher::new();
        self.mers.hash(&mut hash);
        hash.finish()
    }
}

// Various impls for Kminmer.
impl PartialEq for Kminmer {
    fn eq(&self, other: &Kminmer) -> bool {
        self.mers == other.mers
    }
}

impl Eq for Kminmer {
}

impl Hash for Kminmer {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.mers.hash(state);
    }
}

impl Default for Kminmer {
    fn default() -> Self{Kminmer{mers: vec![], start: 0, end: 0, offset: 0, rev: false}}
}

impl Ord for Kminmer {
    fn cmp(&self, other: &Self) -> Ordering {
        self.mers.cmp(&other.mers)
    }
}

impl PartialOrd for Kminmer{
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
