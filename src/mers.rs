// mers.rs
// Contains the "Match", "Offset", and "AlignCand" types, along with driver functions for obtaining reference and query k-min-mers, Hits, Chains, and final coordinates.

use crate::{Chain, Entry, Hit, File, KminmerType, Kminmer, Index, Params, kminmer_mapq, Stats};
use std::borrow::Cow;
use std::cmp;
use std::collections::{hash_map::DefaultHasher, HashMap, HashSet, VecDeque};
use std::hash::{Hash, Hasher};
use std::io::Write;
use dashmap::{DashMap, DashSet};
use rust_seq2kminmers::{KminmersIterator, FH, HashMode};

// A final Match: (Query ID, ref ID, query length, reference length, query start position, query end position, reference start position, reference end position, score, strand direction, MAPQ score)
pub type Match = (String, String, usize, usize, usize, usize, usize, usize, usize, bool, usize);

// An interval that needs to be aligned.
pub type Offset = (usize, usize);

// A tuple of (Query interval, query ID, reference interval, strand direction) for alignment.
pub type AlignCand = (Offset, String, Offset, bool);

// Extract k-min-mers from reference. We don't store k-min-mer objects or hashes in a Vec, but rather immediately insert into the Index.
pub fn ref_extract(seq_id: &str, inp_seq_raw: &[u8], params: &Params, mers_index: &Index) -> usize {
    let l = params.l;
    let k = params.k;
    if inp_seq_raw.len() < l+k-1 {
        return 0;
    }
    let density = params.density as FH;
    let mode = if params.use_hpc {HashMode::Hpc} else {HashMode::Regular};
    let iter = KminmersIterator::new(inp_seq_raw, l, k, density, mode).unwrap();
    let mut count = 0;
    for kminmer in iter {
        //println!("{:?}", kminmer);
        // Add a reference k-min-mer to the Index.
        mers_index.add(kminmer.get_hash(), seq_id, kminmer.start, kminmer.end, kminmer.offset, kminmer.rev);
        count += 1;
        //eprintln!("{}\r", count);
    }
    count
}

// Extract k-min-mers from the query. We need to store Kminmer objects for the query in order to compute Hits.
pub fn extract<'a>(seq_id: &str, inp_seq_raw: &'a [u8], params: &Params) -> Option<KminmersIterator<'a>> {
    let l = params.l;
    let k = params.k;
    if inp_seq_raw.len() < l+k-1 {
        return None;
    }
    let density = params.density as FH;
    let mode = if params.use_hpc {HashMode::Hpc} else {HashMode::Regular};
    return Some(KminmersIterator::new(inp_seq_raw, l, k, density, mode).unwrap());
}

// Generates raw Vecs of Hits by matching query k-min-mers to Entries from the Index.
pub fn chain_hits(query_id: &str, query_it_raw: &mut Option<KminmersIterator>, index: &Index, params: &Params, q_len: usize) -> HashMap<String, Vec<Hit>> {
    let mut hits_per_ref = HashMap::<String, Vec<Hit>>::new();
    let l = params.l;
    let k = params.k;
    if query_it_raw.is_none() {return hits_per_ref;}
    let mut stats = Stats::new(query_id);
    let mut query_it = query_it_raw.as_mut().unwrap().peekable();
    while let Some(q) = query_it.next() {
        let re = index.get(&q.get_hash());
        if let Some(r) = re {
            let mut h = Hit::new(query_id, &q, &r, params);
            h.extend(&mut query_it, index, &r, params, q_len);
            stats.add(&r);
            hits_per_ref.entry(r.id).or_insert(Vec::new()).push(h);
        }
    }
    stats.finalize();
    hits_per_ref
}


// Extract raw Vecs of Hits, construct a Chain, and obtain a final Match (and populate alignment DashMaps with intervals if necessary).
pub fn find_hits(q_id: &str, q_len: usize, q_str: &[u8], ref_lens: &DashMap<String, usize>, mers_index: &Index, params: &Params, aln_coords: &DashMap<String, Vec<AlignCand>>) -> Option<Match> {
    let mut kminmers = extract(q_id, q_str, params);
    let hits_per_ref = chain_hits(q_id, &mut kminmers, mers_index, params, q_len);
    let mut final_matches = Vec::<Match>::new();    
    for e in hits_per_ref.iter() {
        let (r_id, hits_raw) = e;
        let r_len = *ref_lens.get(r_id).unwrap();
        let mut c = Chain::new(&hits_raw);
        let mut v = c.get_match(&r_id, r_len, q_id, q_len, params);
        if let Some(m) = v {final_matches.push(m);}
    }
    let matches_len = final_matches.len();
    return match matches_len {
        0 => None,
        1 => Some(final_matches[0].clone()),
        _ => determine_best_match(&mut final_matches, matches_len),
    };
       /* let (v, c) = &final_matches[0];
        if params.a {
            let (q_coords, r_coords) = c.get_remaining_seqs(&v);
            for i in 0..q_coords.len() {
                let q_coord_tup = q_coords[i];
                let r_coord_tup = r_coords[i];
                aln_coords.get_mut(&v.1).unwrap().push(((r_coord_tup.0, r_coord_tup.1), q_id.to_string(), (q_coord_tup.0, q_coord_tup.1), v.9));
            }
        }*/
}

pub fn determine_best_match(matches: &mut Vec<Match>, matches_len: usize) -> Option<Match> {
    matches.sort_by(|a, b| a.8.cmp(&b.8));
    let mut max_score = &matches[matches_len - 1].8;
    return match max_score == &matches[matches_len - 2].8 {
        true => None,
        false => Some(matches[matches_len - 1].clone()),
    };
}


// Populate a PAF file with final Matches from the queries.
pub fn output_paf(all_matches: &DashSet<(String, Option<Match>)>, paf_file: &mut File, unmap_file: &mut File, params: &Params) {
    for e in all_matches.iter() {
        let id = &e.0;
        let v_opt = &e.1;
        if v_opt.is_none() {
            //write!(unmap_file, "{}\n", id).expect("Error writing line.");
            continue;
        }
        let v = v_opt.as_ref().unwrap();
        let query_len = v.2;
        let ref_len = v.3;
        let q_start = v.4;
        let q_end = v.5;
        let mut r_start = v.6;
        let mut r_end = v.7;
        let score = v.8;
        let rc : &str = match v.9 {true => "-", false => "+"};
        let mapq = v.10;
       /* if mapq <= 1 {
            write!(unmap_file, "{}\n", id).expect("Error writing line.");
        }*/
        let paf_line = format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n", v.0, query_len, q_start, q_end, rc, v.1, ref_len, r_start, r_end, score, ref_len, mapq);
        write!(paf_file, "{}", paf_line).expect("Error writing line.");
    }
}

