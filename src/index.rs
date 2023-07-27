// index.rs
// Contains the "Index" and "Entry" structs, which describe how reference k-min-mers are stored. 

use crate::{KH};
use rust_seq2kminmers::{KminmerType as Kminmer, KminmerHash};
use dashmap::{DashMap, ReadOnlyView};
use std::hash::BuildHasherDefault;
use core::hash::Hasher;

// from https://github.com/Manishearth/trashmap/blob/master/src/lib.rs
#[derive(Default)]
pub struct KnownHasher {
    hash: Option<KH>,
}

impl Hasher for KnownHasher {
    #[inline]
    fn write(&mut self, _: &[u8]) {
        panic!("KnownHasher must be called with known u64 hash values")
    }

#[inline]
    fn write_u32(&mut self, i: u32) {
        debug_assert!(self.hash.is_none());
        self.hash = Some(i as KH);
    }


    #[inline]
    fn write_u64(&mut self, i: u64) {
        debug_assert!(self.hash.is_none());
        self.hash = Some(i as KH);
    }

    #[inline]
    fn finish(&self) -> u64 {
        self.hash.expect("Nothing was hashed")
    }
}

// An Entry object holds information for a reference k-min-mer without storing the minimizer hashes themselves.
#[derive(Clone, Debug, PartialEq)]
pub struct Entry {
    pub id: usize, // Reference ID
    pub start: usize, // Start location
    pub end: usize, // End location
    pub offset: usize, // K-min-mer offset (index in the k-min-mer array)
    pub rc: bool, // Strand direction
}
impl Entry {

    // Create a new Entry.
    pub fn new(id: usize, start: usize, end: usize, offset: usize, rc: bool) -> Self {
        Entry {id, start, end, offset, rc}
    }

    pub fn new_with_mer(id: usize, mer: &KminmerHash) -> Self {
        Entry {id, start: mer.start, end: mer.end, offset:mer.offset, rc: mer.rev}
    }

    // An empty Entry.
    pub fn empty() -> Self {
        Entry {id: 0, start: 0, end: 0, offset: 0, rc: false}
    }

    // Check if this Entry is Empty.
    pub fn is_empty(&self) -> bool {
        self.end == 0
    }
}

// An Index object is a mapping of k-min-mer hashes (see kminmer.rs) to a single Entry (multiple Entries are not allowed).
pub struct Index {
    //pub index: Arc<DashMap<H, Entry, BuildHasherDefault<FxHasher64>>>
    pub index: DashMap<KH, Entry, BuildHasherDefault<KnownHasher>>,
}
impl Index {

    // Create a new Index.
    pub fn new() -> Self {
        //let hasher = BuildHasherDefault::<FxHasher64>::default();
        let hasher = BuildHasherDefault::<KnownHasher>::default();
        let map = DashMap::with_capacity_and_hasher(39821990/* number of kminmers in CHM13V2 with default params*/,
                                                                                         hasher);
        Index {
            index: map,
        }
    }

    pub fn get_count(&self) -> usize {
        self.index.iter().fold(0, |acc, x| (if !x.value().is_empty() {return acc + 1;} else {return acc;}))
    }

    // Add an Entry to the Index. If an Entry for the hash h already exists, insert None to prevent duplicates.
    pub fn add(&self, h: KH, id: usize, start: usize, end: usize, offset: usize, rc: bool) {
        let e = self.index.insert(h, Entry::new(id, start, end, offset, rc));
        if e.is_some() {self.index.insert(h, Entry::empty());}
    }

    pub fn add_with_mer(&self, id: usize, mer: &Kminmer) {
        let h = mer.hash;
        let e = self.index.insert(h, Entry::new_with_mer(id, mer));
        if e.is_some() {self.index.insert(h, Entry::empty());}
    }
}
            

pub struct ReadOnlyIndex {
    pub read_only_index : ReadOnlyView<KH, Entry, BuildHasherDefault<KnownHasher>>
}
impl ReadOnlyIndex {
    pub fn new(index : DashMap<KH, Entry, BuildHasherDefault<KnownHasher>>) -> Self {
        ReadOnlyIndex {
            read_only_index: index.into_read_only()
        }
    }
    // Return the Entry associated with the k-min-mer hash h, or None if none.
    pub fn get(&self, h: &KH) -> Option<&Entry> {
        let e = self.read_only_index.get(h);
        if let Some(r) = e {
            if !r.is_empty() {
                return Some(r);
            }
        }
        None
    }
    
}
