// mers.rs
// Contains the "Match", "Offset", and "AlignCand" types, along with driver functions for obtaining reference and query k-min-mers, Matches, Chains, and final coordinates.

use crate::{Entry, r#match::Match, File, Index, ReadOnlyIndex, Params, Stats, PseudoChainCoords, PseudoChainCoordsTuple, chain::Chain};
use std::borrow::Cow;
use std::cmp;
use std::collections::{hash_map::DefaultHasher, HashMap, HashSet, VecDeque};
use std::hash::{Hash, Hasher};
use std::io::Write;
use dashmap::{DashMap, DashSet};
use rust_seq2kminmers::{KminmersIterator, FH, HashMode, Kminmer};

// An interval that needs to be aligned.
pub type Offset = (usize, usize);

// A tuple of (Query interval, query ID, reference interval, strand direction) for alignment.
pub type AlignCand = (Offset, String, Offset, bool);

// Extract k-min-mers from reference. We don't store k-min-mer objects or hashes in a Vec, but rather immediately insert into the Index.
pub fn ref_extract(ref_idx: usize, inp_seq_raw: &[u8], params: &Params, mers_index: &Index) -> usize {
    let l = params.l;
    let k = params.k;
    if inp_seq_raw.len() < l+k-1 {
        return 0;
    }
    let density = params.density as FH;
    let mode = if params.use_simd {
        if params.use_hpc {HashMode::HpcSimd} else {HashMode::Simd}
    } else { 
        if params.use_hpc {HashMode::Hpc}     else {HashMode::Regular}
    };
    let iter = KminmersIterator::new(inp_seq_raw, l, k, density, mode).unwrap();
    let mut count = 0;
    for kminmer in iter {
        //println!("{:?}", kminmer);
        // Add a reference k-min-mer to the Index.
        mers_index.add(kminmer.get_hash(), ref_idx, kminmer.start, kminmer.end, kminmer.offset, kminmer.rev);
        count += 1;
        //eprintln!("{}\r", count);
    }
    count
}

// Extract k-min-mers from the query. We need to store Kminmer objects for the query in order to compute Matches.
pub fn extract<'a>(seq_id: &str, inp_seq_raw: &'a [u8], params: &Params) -> Option<KminmersIterator<'a>> {
    let l = params.l;
    let k = params.k;
    if inp_seq_raw.len() < l+k-1 {
        return None;
    }
    let density = params.density as FH;
    let mode = if params.use_simd {
        if params.use_hpc {HashMode::HpcSimd} else {HashMode::Simd}
    } else { 
        if params.use_hpc {HashMode::Hpc}     else {HashMode::Regular}
    };
    return Some(KminmersIterator::new(inp_seq_raw, l, k, density, mode).unwrap());
}

// Generates raw Vecs of Matches by matching query k-min-mers to Entries from the Index.
pub fn chain_matches(query_id: &str, query_it_raw: &mut Option<KminmersIterator>, index: &ReadOnlyIndex, params: &Params, q_len: usize) -> HashMap<usize, Vec<Match>> {
    let mut matches_per_ref = HashMap::<usize, Vec<Match>>::new();
    let l = params.l;
    let k = params.k;
    if query_it_raw.is_none() {return matches_per_ref;}
    let mut stats = Stats::new(query_id);
    let mut query_it = query_it_raw.as_mut().unwrap().peekable();
    while let Some(q) = query_it.next() {
        let re = index.get(&q.get_hash());
        if let Some(r) = re {
            let mut h = Match::new(&q, &r);
            h.extend(&mut query_it, index, &r);
            stats.add(&r);
            matches_per_ref.entry(r.id).or_insert(Vec::new()).push(h);
        }
    }
    stats.finalize();
    matches_per_ref
}


// Extract raw Vecs of Matches, construct a Chain, and obtain a final Match (and populate alignment DashMaps with intervals if necessary).
pub fn find_matches(q_id: &str, q_len: usize, q_str: &[u8], ref_map: &DashMap<usize, (String, usize)>, mers_index: &ReadOnlyIndex, params: &Params, aln_coords: &DashMap<String, Vec<AlignCand>>) -> Option<String> {
    let mut kminmers = extract(q_id, q_str, params);
    let matches_per_ref = chain_matches(q_id, &mut kminmers, mers_index, params, q_len);
    let mut all_pseudocoords = Vec::<(PseudoChainCoordsTuple, usize)>::new();    
    for e in matches_per_ref.iter() {
        let (r_id, matches_raw) = e;
        let mut c = Chain::new(matches_raw);
        let mut tp = c.get_match(&params);
        let mut len_score = c.len();
        if let Some(t) = tp {all_pseudocoords.push(((*r_id, t), len_score));}
    }
    let coords_count = all_pseudocoords.len();
    let paf_line = match coords_count {
        0 => None,
        1 => Some(find_coords(q_id, q_len, ref_map,  &all_pseudocoords[0].0, all_pseudocoords[0].1)),
        _ => determine_best_match(q_id, q_len, ref_map, &all_pseudocoords, coords_count),
    };
    if params.a {
        // TODO that code needs to be updated.. seems we don't remember the full chain anymore, and
        // get_remaining_seqs needs a very particular format
        //let (v, c) = m;
        /*let (q_coords, r_coords) = c.get_remaining_seqs(&v);
        for i in 0..q_coords.len() {
            let q_coord_tup = q_coords[i];
            let r_coord_tup = r_coords[i];
            aln_coords.get_mut(&v.1).unwrap().push(((r_coord_tup.0, r_coord_tup.1), q_id.to_string(), (q_coord_tup.0, q_coord_tup.1), v.9));
        }*/
    }
    paf_line
}

pub fn determine_best_match(q_id: &str, q_len: usize, ref_map: &DashMap<usize, (String, usize)>, all_pseudocoords: &Vec<(PseudoChainCoordsTuple, usize)>, coords_count: usize) -> Option<String> {
    let (max_i, next_max_i, max_count, next_max_count) = find_largest_two_chains(all_pseudocoords, coords_count);
    if max_count == next_max_count {return None;}
    else {return Some(find_coords(q_id, q_len, ref_map, &all_pseudocoords[max_i].0, all_pseudocoords[max_i].1));}
}

pub fn find_largest_two_chains(all_pseudocoords: &Vec<(PseudoChainCoordsTuple, usize)>, coords_count: usize) -> (usize, usize, usize, usize) {
    let mut max = 0;
    let mut max_count = 0;
    let mut second_max = 0;
    let mut second_max_count = 0;
    for i in 0..coords_count {
        let (r_id, coord) = all_pseudocoords[i].0;
        let count = coord.5;
        if count > max_count {
            second_max = max;
            second_max_count = max_count;
            max = i;
            max_count = count;
        } else if count > second_max_count {
            second_max = i;
            second_max_count = count;
        }
    }
    (max, second_max, max_count, second_max_count)
}

pub fn find_coords(q_id: &str, q_len: usize, ref_map: &DashMap<usize, (String, usize)>, t: &PseudoChainCoordsTuple, len_score: usize) -> String {
    let (r_idx, coords) = *t;
    let rtup = ref_map.get(&r_idx).unwrap();
    let r_id = &rtup.0;
    let r_len = rtup.1;
    let (rc, q_start, q_end, r_start, r_end, score, mapq) = coords;
    let mut final_r_start = r_start;
    let mut final_r_end = r_end;
    let mut final_q_start = q_start;
    let mut final_q_end = q_end;
    let mut exc_s = 0;
    let mut exc_e = 0;

    if !rc {
        if r_start >= q_start {
            final_r_start = r_start - q_start;
            exc_s = q_start;
        }
        else {
            final_r_start = 0;
            exc_s = r_start;
        }
        if r_end + (q_len - q_end - 1) <= r_len - 1 {
            final_r_end = r_end + (q_len - q_end - 1);
            exc_e = q_len - q_end - 1;
        }
        else {
            final_r_end = r_len - 1;
            exc_e = r_len - r_end - 1;
        }
    }
    else {
        if r_end + q_start <= r_len - 1 {
            final_r_end = r_end + q_start;
            exc_s = q_start;
        }
        else {
            final_r_end = r_len - 1;
            exc_s = r_len - r_end - 1;
        }
        if r_start >= (q_len - q_end - 1) {
            final_r_start = r_start - (q_len - q_end - 1);
            exc_e = q_len - q_end - 1;
        }
        else {
            final_r_start = 0;
            exc_e = r_start;
        }
    }
    final_q_start = q_start - exc_s;
    final_q_end = q_end + exc_e;
    let rc_s : &str = match rc {true => "-", false => "+"};
    let paf_line = format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", q_id, q_len, final_q_start, final_q_end, rc_s, r_id, r_len, final_r_start, final_r_end, score, r_len, mapq);
    paf_line
}
