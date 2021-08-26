use std::collections::VecDeque;
use crate::DashMap;
use crate::File;


use std::io::Write;

use std::collections::HashMap;
use std::collections::HashSet;
use super::Params;
pub type Match = (String, String, usize, usize, usize, usize, usize, usize, usize, bool);
pub type Mer = (u64, usize, usize, usize, bool);
#[derive(Clone, Debug)]
pub struct Hit {
    query_id: String,
    ref_id: String,
    query_s: usize,
    query_e: usize,
    ref_s: usize,
    ref_e: usize,
    hit_count: usize,
    is_rc: bool,
    query_span: usize,
    ref_span: usize,
    query_offset: usize,
    ref_offset: usize,

}

const seq_nt4_table: [u8; 256] =
   [0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4];
// copy from http://www.cse.yorku.ca/~oz/hash.html:

pub fn hash(mut key: u64, mask: u64) -> u64 {
        key = (!key + (key << 21)) & mask;
        key = key ^ key >> 24;
        key = ((key + (key << 3)) + (key << 8)) & mask;
        key = key ^ key >> 14;
        key = ((key + (key << 2)) + (key << 4)) & mask;
        key = key ^ key >> 28;
        key = (key + (key << 31)) & mask;
        return key;
}

pub fn update_window(q: &mut VecDeque<u64>, q_pos: &mut VecDeque<usize>, q_min_val: u64, q_min_pos: i32, new_strobe_hashval: u64, i: usize, new_minimizer: bool) -> (u64, i32, bool) {
    q.pop_front();
    let popped_index = q_pos.pop_front();
    q.push_back(new_strobe_hashval);
    q_pos.push_back(i);
    let mut min_val = q_min_val;
    let mut min_pos = q_min_pos;
    let mut new_minim = new_minimizer;
    if min_pos == popped_index.unwrap() as i32 {
        min_val = u64::max_value();
        min_pos = i as i32;
        for j in (0..q.len()).rev() {
            if q[j] < min_val {
                min_val = q[j];
                min_pos = q_pos[j] as i32;
                new_minim = true;
            }
        }
    }
    else if new_strobe_hashval < min_val { // the new value added to queue is the new minimum
        min_val = new_strobe_hashval;
        min_pos = i as i32;
        new_minim = true;
    }
    (min_val, min_pos, new_minim)
}

pub fn extract_mers(seq: &[u8], params: &Params) -> (Vec<u64>, Vec<usize>) {
    let l = params.l;
    let hash_bound = ((params.density as f64) * 4_usize.pow(l as u32) as f64) as u64; 
    let s = params.s;
    let wmin = params.wmin;
    let wmax = params.wmax;
    let smask : u64 = ((1 as u64) << 2*s) - 1;
    let lmask : u64 = ((1 as u64) << 2*l) - 1;
    let t = 3;
    let mut hash_count = 0;
    let mut seq_hashes = Vec::new();
    let mut pos_to_seq_coord = Vec::new();
    let mut qs = VecDeque::<u64>::new();
    let mut qs_pos = VecDeque::<usize>::new();
    let seq_len = seq.len();
    let mut qs_size = 0;
    let mut qs_min_val = u64::max_value();
    let mut qs_min_pos : i32 = -1;
    let mut xl : [u64; 2] = [0; 2];
    let mut xs : [u64; 2] = [0; 2];
    let mut lp = 0;
    let lshift : u64 = (l as u64 - 1) * 2;
    let sshift : u64 = (s as u64 - 1) * 2;
    for i in 0..seq.len() {
        let c = seq_nt4_table[seq[i] as usize];
        if c < 4 {
            xl[0] = (xl[0] << 2 | c as u64) & lmask;                  // forward strand
            xl[1] = xl[1] >> 2 | ((3 - c) as u64) << lshift;  // reverse strand
            xs[0] = (xs[0] << 2 | c as u64) & smask;                  // forward strand
            xs[1] = xs[1] >> 2 | ((3 - c) as u64) << sshift;  // reverse strand
            lp += 1;
            if params.s != 0 { //ksyncmer or kstrobemer
                if lp >= s {
                    let ys : u64 = match xs[0] < xs[1]{
                        true => xs[0],
                        false => xs[1]
                    };
                    let hash_s = hash(ys, smask);
                    if qs_size < l - s {
                        qs.push_back(hash_s);
                        qs_pos.push_back(i - s + 1);
                        qs_size += 1;
                    }
                    else if qs_size == l - s {
                        qs.push_back(hash_s);
                        qs_pos.push_back(i - s + 1);
                        qs_size += 1;
                        for j in 0..qs_size {
                            if qs[j] < qs_min_val {
                                qs_min_val = qs[j];
                                qs_min_pos = qs_pos[j] as i32;
                            }
                        }
                        if qs_min_pos == qs_pos[t-1] as i32 {
                            let yl : u64 = match xl[0] < xl[1]{
                                true => xl[0],
                                false => xl[1]
                            };
                            let hash_l = hash(yl, lmask);
                            seq_hashes.push(hash_l);
                            pos_to_seq_coord.push(i - l + 1);
                            hash_count += 1;
                        }
                    }
                    else {
                        let mut new_minimizer = false;
                        let tuple = update_window(&mut qs, &mut qs_pos, qs_min_val, qs_min_pos, hash_s, i - s + 1, new_minimizer);
                        qs_min_val = tuple.0; qs_min_pos = tuple.1; new_minimizer = tuple.2;
                        if qs_min_pos == qs_pos[t-1] as i32 { // occurs at t:th position in k-mer
                            let yl : u64 = match xl[0] < xl[1] {
                                true => xl[0],
                                false => xl[1]
                            };
                            let hash_l = hash(yl, lmask);
                            seq_hashes.push(hash_l);
                            pos_to_seq_coord.push(i - l + 1);
                            hash_count += 1;
                        }
                    }
                }
            }
            else { //kminmer
                let yl : u64 = match xl[0] < xl[1] {
                    true => xl[0],
                    false => xl[1]
                };
                let hash_l = hash(yl, lmask);
                if hash_l <= hash_bound {
                    seq_hashes.push(hash_l);
                    pos_to_seq_coord.push(i);
                    hash_count += 1;
                }
            }
        } else {
            qs_min_val = u64::max_value();
            qs_min_pos = -1;
            lp = 0; xs = [0; 2]; xl = [0; 2];
            qs_size = 0;
            qs.clear();
            qs_pos.clear();
        }
    }
    return (seq_hashes, pos_to_seq_coord)
} 

pub fn get_next_strobe(string_hashes: &Vec<u64>, strobe_hashval: u64, w_start: usize, w_end: usize, q: u64) -> (i32, u64) {
    let mut min_val = u64::max_value();
    let mut strobe_pos_next = -1;
    let mut strobe_hashval_next = u64::max_value();
    for i in w_start..w_end + 1 {
        let res : u64 = (strobe_hashval + string_hashes[i]) & q;
        if res < min_val {
            min_val = res;
            strobe_pos_next = i as i32;
            strobe_hashval_next = string_hashes[i];
        }
    }
    return (strobe_pos_next, strobe_hashval_next);
}

pub fn get_randstrobe(i: usize, wmin: usize, wmax: usize, num_hashes: usize, string_hashes: &Vec<u64>, pos_to_seq_coord: &Vec<usize>, q: u64) -> Option<(u64, usize, usize)> {
    let strobe_pos_next;
    let strobe_hashval_next;
    if i + wmax < num_hashes {
        let w_start = i + wmin;
        let w_end = i + wmax;
        let strobe_hash = string_hashes[i];
        let tuple = get_next_strobe(string_hashes, strobe_hash, w_start, w_end, q);
        strobe_pos_next = tuple.0; strobe_hashval_next = tuple.1;
    }
    else if (i + wmin + 1 < num_hashes) && (num_hashes <= i + wmax) {
        let w_start = i + wmin;
        let w_end = num_hashes - 1;
        let strobe_hash = string_hashes[i];
        let tuple = get_next_strobe(string_hashes, strobe_hash, w_start, w_end, q);
        strobe_pos_next = tuple.0; strobe_hashval_next = tuple.1;
    }
    else {return None;}
    let hash_randstrobe = string_hashes[i]/2 + strobe_hashval_next/3;
    let seq_pos_strobe1 = pos_to_seq_coord[i];
    let seq_pos_strobe2 = pos_to_seq_coord[strobe_pos_next as usize];
    let s = (hash_randstrobe, seq_pos_strobe1, seq_pos_strobe2);
    return Some(s);
}

pub fn seq_to_kmers(seq: &[u8], id: &str, params: &Params, read: bool) -> (usize, Vec<Mer>) {
    let k = params.k;
    let l = params.l;
    let s = params.s;
    let wmin = params.wmin;
    let wmax = params.wmax;
    let mut mers = Vec::<Mer>::new();
    let read_length = seq.len();
    if read_length < wmax {return (seq.len(), mers);}
    let lmask : u64 = ((1 as u64) << 2*l) - 1;
    let smask : u64 = ((1 as u64) << 2*s) - 1;
    let q = 2_u64.pow(16) - 1;
    let t = 3;
    let (mut string_hashes, mut pos_to_seq_coord) = extract_mers(seq, params);
    let num_hashes = string_hashes.len();
    for i in 0..num_hashes {
        let mut kmer_hash : u64 = 0;
        let mut kmer_hash_rev : u64 = 0;
        let mut kmer_start = 0;
        let mut kmer_end = 0;
        let mut reverse = false;
        if i + k - 1 < num_hashes {
            //if read && (string_hashes[i+k-1] < string_hashes[i]) {reverse = true;}
            for j in i..i+k {
                if params.use_strobe {
                    let s = get_randstrobe(j, wmin, wmax, num_hashes, &string_hashes, &pos_to_seq_coord, q);
                    if s.is_none() {break;}
                    let s_unw = s.unwrap();
                    kmer_hash += s_unw.0;
                    if j == i + k - 1 {kmer_end = s_unw.2;}
                }
                else {
                    kmer_hash += string_hashes[j];
                }
            }
                kmer_start = pos_to_seq_coord[i];
                kmer_end = pos_to_seq_coord[i + k - 1];
        }
        if kmer_end != 0 {
            mers.push((kmer_hash, kmer_start, kmer_end + l, i, false));
        }   
    }
    /*if read {
        string_hashes.reverse();
        pos_to_seq_coord.reverse();
        for i in 0..num_hashes {
            pos_to_seq_coord[i] = read_length - pos_to_seq_coord[i] - l;
        }
        for i in 0..num_hashes {
            let mut kmer_hash : u64 = 0;
            let mut kmer_start = 0;
            let mut kmer_end = 0;
            if i + k - 1 < num_hashes {
                for j in i..i+k {
                    if params.use_strobe {
                        let s = get_randstrobe(j, wmin, wmax, num_hashes, &string_hashes, &pos_to_seq_coord, q);
                        if s.is_none() {break;}
                        let s_unw = s.unwrap();
                        kmer_hash += s_unw.0;
                        if j == i + k - 1 {kmer_end = s_unw.2;}
                    }
                    else {kmer_hash += string_hashes[j];}
                }
                kmer_start = pos_to_seq_coord[i];
                kmer_end = pos_to_seq_coord[i + k - 1];
            }
            if kmer_end != 0 {mers.push((kmer_hash, kmer_start, kmer_end, i, true));}
        }
    }*/
    return (seq.len(), mers);
}

pub fn find_hits(query_id: &str, query_len: usize, query_mers: &Vec<Mer>, ref_lens: &DashMap<String, usize>, mers_index: &DashMap<String, DashMap<u64, (Mer, usize)>>, l: usize, k: usize) -> Vec<Match> {
    let mut hits_per_ref = HashMap::<String, Vec<Hit>>::new();
    let mut hit_count_all = 0;
    let mut i = 0;
    let mut i_rev = 0;
    println!("{}", query_id);
    while i < query_mers.len() {
        if i < query_mers.len() {
            let q = query_mers[i];
            let mut h = Hit {query_id: query_id.to_string(), ref_id: String::new(), query_s: q.1, query_e: q.2, ref_s: 0, ref_e: 0, hit_count: 0, is_rc: q.4, query_span: 0, ref_span: 0, query_offset: 0, ref_offset: 0};
            let mut max_extend_offset = 0;
            for index in mers_index.iter() {
                let ref_id = index.key().to_string();
                let ref_mer_count = index.value();
                let mer_entry = ref_mer_count.get(&q.0);
                if mer_entry.is_some() {
                    let tup = mer_entry.unwrap();
                    let (mer, count) = tup.value();
                    if *count > 1 {continue;}
                    // figure out if + or -
                    println!("{:?}\t{:?}", q, mer);
                    h.query_offset = q.3;
                    h.ref_offset = mer.3;
                    h.ref_id = ref_id.to_string();
                    h.ref_s = mer.1;
                    h.query_e = q.2;
                    h.ref_e = mer.2;
                    h.hit_count += 1;
                    h.is_rc = q.4;
                    h.query_span = h.query_e - h.query_s + 1; 
                    h.ref_span = h.ref_e - h.ref_s + 1; 
                    hits_per_ref.entry(ref_id.to_string()).or_insert(vec![]).push(h.clone());
                    hit_count_all += 1;
                }

            }
        }
        i += 1;
    }
    let mut final_matches = Vec::<Match>::new();
    let prev_key : u64 = 0;
    
    for (key, val) in hits_per_ref.into_iter() {
        let mut final_rc_mappings = Vec::<Vec::<Hit>>::new();
        let mut final_fwd_mappings = Vec::<Vec::<Hit>>::new();
        let mut max_fwd_mapping = 0;
        let mut max_fwd_mapping_score = 0;
        let mut max_rc_mapping = 0;
        let mut max_rc_mapping_score = 0;
        let ref_id = key;
        let mut hits_raw = val.to_vec();
        let ref_len = *ref_lens.get(&ref_id).unwrap().value();
            //println!("_____________________{} to {}_____________\nRC", query_id, ref_id);
        let prev_q_start = 0;
        let rc_mono = true;
        let fwd_hit_avg = 0.0;
        let rc_hit_avg = 0.0;
        let fwd_mono = true;
        let mut write_fwd = true;
        let mut write_rc = true;
        let mut s = 0;
        let mut e = 0;
        let mut span = 0;
        let mut score = 0;
        /*let rc_hits : Vec<&Hit> = hits.iter().filter(|hit| hit.is_rc == true && hit.query_id == query_id).collect();
        let rc_hits_cleared = Vec::<&Hit>::new();
        let fwd_hits : Vec<&Hit> = hits.iter().filter(|hit| hit.is_rc == false && hit.query_id == query_id).collect();
        let fwd_hits_cleared = Vec::<&Hit>::new();
        let rc_mappings = Vec::<Vec::<&Hit>>::new();
        let fwd_mappings = Vec::<Vec::<&Hit>>::new();*/
        let mut fwd_off = 0;
        let mut rc_off = 0;
        let mut filtered_hits_fwd = Vec::<&Hit>::new();
        let mut filtered_hits_rc = Vec::<&Hit>::new();
        let mut hits = Vec::<&Hit>::new();
        let mut max_fwd_off = 0;
        let mut max_rc_off = 0;
        if hits_raw.len() != 0 {
            hits_raw.sort_by(|a, b| a.ref_offset.cmp(&b.ref_offset));
            let mut dir = true;
            if hits_raw.len() > 1 {
                println!("FIRST OFFSET {} SECOND OFFSET {}, DIR {:?}", hits_raw[0].ref_offset, hits_raw[1].ref_offset, dir);
                dir = (hits_raw[1].query_offset >= hits_raw[0].query_offset);
            }
            let mut prev = &hits_raw[0];
            let mut prev_offset = hits_raw[0].query_offset;
            println!("FIRST OFFSET {} DIR {:?}", hits_raw[0].query_offset, dir);
            filtered_hits_fwd.push(&hits_raw[0]);
            filtered_hits_rc.push(&hits_raw[0]);
            for i in 1..hits_raw.len() {
                println!("{}", hits_raw[i].query_offset);
                if dir == true {
                    if hits_raw[i].query_offset >= prev_offset {
                        filtered_hits_fwd.push(&hits_raw[i]);
                        fwd_off += 1;
                        if max_fwd_off < fwd_off {max_fwd_off = fwd_off;}
                    }
                    else {
                        dir = false;
                        fwd_off = 0;
                    }
                }
                else {
                    if hits_raw[i].query_offset <= prev_offset {
                        filtered_hits_rc.push(&hits_raw[i]);
                        rc_off += 1;
                        if max_rc_off < rc_off {max_rc_off = rc_off;}
                    }
                    else {
                        dir = true;
                        rc_off = 0;
                    }
                }
                prev = &hits_raw[i];
                prev_offset = prev.query_offset;
            }
            if max_rc_off < max_fwd_off {hits = filtered_hits_fwd.clone();}
            else {hits = filtered_hits_rc.clone();}
        }
        println!("HITS LEN {}\tFWD OFF {}\tRC OFF {}\t", hits.len(), max_fwd_off, max_rc_off);
        s = hits[0].ref_s;
        e = hits[hits.len()-1].ref_e;
        span = e - s;
        score = hits.len() * span;
        let mut final_ref_s = hits[0].ref_s;
        let mut final_ref_e = hits[hits.len()-1].ref_e;
        let mut final_query_s = hits[0].query_s;
        let mut final_query_e = hits[hits.len()-1].query_e;
        if max_rc_off > max_fwd_off {
            final_query_s = hits[hits.len()-1].query_s;
            final_query_e = hits[0].query_e;
        }

        let ref_id = &hits[0].ref_id;
        if max_fwd_off > max_rc_off && hits.len() != 0 {
            let v = (query_id.to_string(), ref_id.to_string(), query_len, ref_len, final_query_s, final_query_e, final_ref_s, final_ref_e, score, false);
            final_matches.push(v);
        }
        else if max_rc_off > max_fwd_off && hits.len() != 0 {
            let v = (query_id.to_string(), ref_id.to_string(), query_len, ref_len, final_query_s, final_query_e, final_ref_s, final_ref_e, score, true);
            final_matches.push(v);
        }
    }
    return final_matches
}
/*


            // -----------------------FINDING RC MAPPINGS-----------------------
        if rc_hits.len() != 0 {
            let mut prev_rc_hit = rc_hits[0];
            let mut current_mapping = vec![prev_rc_hit];
            for i in 1..rc_hits.len() {
                let curr_rc_hit = rc_hits[i];
                if curr_rc_hit.ref_s < prev_rc_hit.ref_s || /*curr_rc_hit.hit_count == 1 ||*/ (prev_rc_hit.ref_e  + query_len < curr_rc_hit.ref_s) {
                    if current_mapping.len() != 0 {
                        final_rc_mappings.push(current_mapping.into_iter().cloned().collect());
                    }
                    current_mapping = Vec::<&Hit>::new();
                    continue;
                } 
                current_mapping.push(curr_rc_hit);
                prev_rc_hit = curr_rc_hit;               
            }
            final_rc_mappings.push(current_mapping.into_iter().cloned().collect());
        }
        // -----------------------FINDING FWD MAPPINGS-----------------------
        if fwd_hits.len() != 0 {
            let mut prev_fwd_hit = fwd_hits[0];
            let mut current_mapping = vec![prev_fwd_hit];
            for i in 1..fwd_hits.len() {
                let curr_fwd_hit = fwd_hits[i];
                if curr_fwd_hit.ref_s < prev_fwd_hit.ref_s || /* curr_fwd_hit.hit_count == 1  ||*/ (prev_fwd_hit.ref_e  + query_len < curr_fwd_hit.ref_s) {
                    if current_mapping.len() != 0 {
                        final_fwd_mappings.push(current_mapping.into_iter().cloned().collect());
                    }
                    current_mapping = Vec::<&Hit>::new();
                    continue;
                }
                current_mapping.push(curr_fwd_hit);
                prev_fwd_hit = curr_fwd_hit;
            }
            final_fwd_mappings.push(current_mapping.into_iter().cloned().collect());
        }
        for i in 0..final_fwd_mappings.len() {
            let mapping = &final_fwd_mappings[i];
            if mapping.is_empty() {continue;}
            let mapping_s = mapping.iter().min_by(|a, b| a.query_s.cmp(&b.query_s)).unwrap().query_s;
            let mapping_e = mapping.iter().max_by(|a, b| a.query_e.cmp(&b.query_e)).unwrap().query_e;
            let mapping_span = mapping_e - mapping_s;
            let mapping_score = mapping.len() * mapping_span;
            //for hit in mapping.iter() {mapping_score += (hit.query_e - hit.query_s);}
            if mapping_score > max_fwd_mapping_score {max_fwd_mapping_score = mapping_score; max_fwd_mapping = i;}
            /*println!("MAPPING {} SCORE {}", i+1, mapping_score);
            for hit in mapping.iter() {
                println!("{:?}", hit);
            }*/
        }
        for i in 0..final_rc_mappings.len() {
            let mapping = &final_rc_mappings[i];
            if mapping.is_empty() {continue;}
            let mapping_s = mapping.iter().min_by(|a, b| a.query_s.cmp(&b.query_s)).unwrap().query_s;
            let mapping_e = mapping.iter().max_by(|a, b| a.query_e.cmp(&b.query_e)).unwrap().query_e;
            let mapping_span = mapping_e - mapping_s;
            let mapping_score = mapping.len() * mapping_span;
            /*for hit in mapping.iter() {
                mapping_score += (hit.query_e - hit.query_s);
            }*/
            if mapping_score > max_rc_mapping_score {max_rc_mapping_score = mapping_score; max_rc_mapping = i;}
            /*println!("MAPPING {} SCORE {}", i+1, mapping_score);
            for hit in mapping.iter() {
                println!("{:?}", hit);
            }*/
        }  
        let mut final_ref_s = 0;
        let mut final_ref_e = 0;
        let mut final_query_s = 0;
        let mut final_query_e = 0;
        if max_fwd_mapping_score > max_rc_mapping_score {
            let max_mapping = final_fwd_mappings[max_fwd_mapping].to_vec();
            if max_mapping[0].ref_s < max_mapping[0].query_s {
                final_ref_s = 0;
                final_ref_e = query_len - 1;
                final_query_s = max_mapping[0].query_s - max_mapping[0].ref_s;
                final_query_e = max_mapping[0].query_s - max_mapping[0].ref_s + query_len - 1;

            }
            else {
                final_query_s = 0;
                final_query_e = query_len - 1;
                final_ref_s = max_mapping[0].ref_s - max_mapping[0].query_s;
                final_ref_e = max_mapping[0].ref_s - max_mapping[0].query_s + query_len - 1; 
            }
            let ref_id = &max_mapping[0].ref_id;
            let v = (query_id.to_string(), ref_id.to_string(), query_len, ref_len, final_query_s, final_query_e, final_ref_s, final_ref_e, max_fwd_mapping_score, false);
            final_matches.push(v);
        }
        else if max_fwd_mapping_score < max_rc_mapping_score {
            let max_mapping = final_rc_mappings[max_rc_mapping].to_vec();
            if max_mapping[0].ref_s < max_mapping[0].query_s {
                final_ref_s = 0;
                final_ref_e = query_len - 1;
                final_query_s = max_mapping[0].query_s - max_mapping[0].ref_s;
                final_query_e = max_mapping[0].query_s - max_mapping[0].ref_s + query_len - 1;

            }
            else {
                final_query_s = 0;
                final_query_e = query_len - 1;
                final_ref_s = max_mapping[0].ref_s - max_mapping[0].query_s;
                final_ref_e = max_mapping[0].ref_s - max_mapping[0].query_s + query_len - 1; 
            }
            let ref_id = &max_mapping[0].ref_id;
            let v = (query_id.to_string(), ref_id.to_string(), query_len, ref_len, final_query_s, final_query_e, final_ref_s, final_ref_e, max_rc_mapping_score, true);
            //println!("{:?}", v);
            final_matches.push(v);
        }
    }
    return final_matches
}*/

pub fn output_paf(all_matches: &Vec<(Vec<Match>, String)>, paf_file: &mut File) {
    for (matches, id) in all_matches.iter() {
        let v = matches.iter().max_by(|a, b| a.8.cmp(&b.8)).unwrap();
        let query_id = v.0.to_string();
        let ref_id = v.1.to_string();
        let query_len = v.2;
        let ref_len = v.3;
        let query_s = v.4;
        let query_e = v.5;
        let ref_s = v.6;
        let ref_e = v.7;
        let score = v.8;
        let rc : String = match v.9 {true => "-".to_string(), false => "+".to_string()};
        let paf_line = format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n", query_id, query_len, query_s, query_e, rc, ref_id, ref_len, ref_s, ref_e, score, ref_len, "255");
        write!(paf_file, "{}", paf_line).expect("Error writing line.");
    }
}

