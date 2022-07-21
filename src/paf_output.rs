use bio::alignment::pairwise::banded::*;
use bio::alignment::AlignmentOperation::*;
use crate::HashEntry;
//     struct QueryMatch {kminmer: Kmer, target: Vec<HashEntry>, query: HashEntry, unique: bool};
// struct HashEntry {origin: String, seqlen: u32, shift: (u16, u16), seq_rev: bool} 

pub fn determine_orientation(query: &HashEntry, target: &HashEntry) -> &'static str {
    if query.seq_rev == target.seq_rev {return "+";}
    else {return "-";}
}   

pub fn align_bases(query: &HashEntry, target: &HashEntry) -> (String, usize, usize) {
    let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
    // gap open score: -5, gap extension score: -1
    let l = 12;  // kmer match length
    let w = 6;  // Window size for creating the band
    let mut aligner = Aligner::new(-5, -1, score, l, w);
    let alignment = aligner.semiglobal(query.seq.as_bytes(), target.seq.as_bytes());
    let cigar = alignment.cigar(false);
    //println!("{}", alignment.pretty(query.seq.as_bytes(), target.seq.as_bytes()));
    let residue_match = alignment.operations.iter().filter(|&n| *n == Match).count();
    let mut block_len = query.seq.len();
    if target.seq.len() > query.seq.len() {block_len = target.seq.len();}
    (cigar, block_len, residue_match)
}