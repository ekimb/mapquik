use petgraph::graph::DiGraph;
use petgraph::visit::EdgeRef;
use std::path::PathBuf;
use std::fs::File;
use std::error::Error;
use std::io::Write;
use super::Kmer;
use crate::kmer_vec::get;
use std::collections::HashMap;
use strsim::levenshtein;
use petgraph::graph::NodeIndex;
use crate::utils::revcomp;
use crate::HashEntry;
//     struct QueryMatch {kminmer: Kmer, target: Vec<HashEntry>, query: HashEntry, unique: bool};
// struct HashEntry {origin: String, seqlen: u32, shift: (u16, u16), seq_rev: bool} 

pub fn determine_orientation(query: &HashEntry, target: &HashEntry) -> (&'static str) {
    if query.seq_rev == target.seq_rev {return "+";}
    else {return "-";}
}   