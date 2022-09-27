// index.rs
// Contains the "Index" and "Entry" structs, which describe how reference k-min-mers are stored. 

use crate::Kminmer;
use dashmap::DashMap;
use std::sync::Arc;

// An Entry object holds information for a reference k-min-mer without storing the minimizer hashes themselves.
#[derive(Clone, Debug, PartialEq)]
pub struct Entry {
    pub id: String, // Reference ID
    pub start: usize, // Start location
    pub end: usize, // End location
    pub offset: usize, // K-min-mer offset (index in the k-min-mer array)
    pub rc: bool, // Strand direction
}
impl Entry {

    // Create a new Entry.
    pub fn new(id: &str, start: usize, end: usize, offset: usize, rc: bool) -> Self {
        Entry {id: id.to_string(), start, end, offset, rc}
    }

    // Create a new Entry given a tuple.
    pub fn from_tuple(t: (String, usize, usize, usize, bool)) -> Self {
        Entry {id: t.0.to_string(), start: t.1, end: t.2, offset: t.3, rc: t.4}
    }

    // Output a Raw tuple.
    pub fn expand(&self) -> (String, usize, usize, usize, bool) {
        (self.id.to_string(), self.start, self.end, self.offset, self.rc)
    }

    // An empty Entry.
    pub fn empty() -> Self {
        Entry {id: String::new(), start: 0, end: 0, offset: 0, rc: false}
    }

    // Check if this Entry is Empty.
    pub fn is_empty(&self) -> bool {
        self.id.is_empty()
    }
}

// An Index object is a mapping of k-min-mer hashes (see kminmer.rs) to a single Entry (multiple Entries are not allowed).
pub struct Index {
    pub index: Arc<DashMap<u64, Entry>>
}
impl Index {

    // Create a new Index.
    pub fn new() -> Self {
        Index {index: Arc::new(DashMap::new())}
    }

    // Return the Entry associated with the k-min-mer hash h, or None if none.
    pub fn get(&self, h: &u64) -> Option<Entry> {
        let e = self.index.get(h);
        if e.is_some() {
            let r = e.unwrap().clone();
            if !r.is_empty() {return Some(r);}
        }
        None
    }

    // Add an Entry to the Index. If an Entry for the hash h already exists, insert None to prevent duplicates.
    pub fn add(&self, h: u64, id: &str, start: usize, end: usize, offset: usize, rc: bool) {
        let e = self.index.insert(h, Entry::new(id, start, end, offset, rc));
        if e.is_some() {self.index.insert(h, Entry::empty());}
    }

    // Get entry associated with a Kminmer q by getting its hash, or an empty Entry if None.
    pub fn get_entry(&self, q: &Kminmer) -> (bool, Entry) {
        let e = self.get(&q.get_hash());
        let r = match e.is_some() {
            true => (true, e.unwrap()),
            false => (false, Entry::empty())
        };
        r
    }
}
