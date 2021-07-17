use std::collections::VecDeque;
use crate::DashMap;
use crate::File;
use dashmap::ReadOnlyView;
use std::io::Write;
use std::cmp::Ordering;
use super::{MersVector, MersVectorReduced, MersVectorRead, PosIndex, KmerLookup, KmerLookupMod};
use super::Params;

#[derive(Clone, Debug)]
pub struct Hit {
    query_id: u64,
    ref_id: u64,
    query_s: usize,
    query_e: usize,
    ref_s: usize,
    ref_e: usize,
    hit_count: usize,
    is_rc: bool,
}
#[derive(Clone, Debug)]
pub struct NAM {
    query_id: u64,
    ref_id: u64,
    query_s: usize,
    query_e: usize,
    query_last_hit_pos: usize,
    ref_s: usize,
    ref_e: usize,
    ref_last_hit_pos: usize,
    n_hits: usize,
    previous_query_start: usize,
    previous_ref_start: usize,
    is_rc: bool,
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

pub fn update_window(mut q: &mut VecDeque<u64>, mut q_pos: &mut VecDeque<usize>, q_min_val: u64, q_min_pos: i32, new_strobe_hashval: u64, i: usize, new_minimizer: bool) -> (u64, i32, bool) {
    q.pop_front();
    let mut popped_index = q_pos.pop_front();
    q.push_back(new_strobe_hashval);
    q_pos.push_back(i);
    let mut min_val = q_min_val;
    let mut min_pos = q_min_pos;
    let mut new_minim = new_minimizer;
    if min_pos == popped_index.unwrap() as i32 {
        min_val = u64::max_value();
        min_pos = i as i32;
        for j in (0..q.len()).rev() {
            if (q[j] < min_val) {
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

pub fn extract_syncmers(seq: &[u8], params: &Params) -> (Vec<u64>, Vec<usize>) {
    let k = params.l;
    let s = params.s;
    let wmin = params.wmin;
    let wmax = params.wmax;
    let smask : u64 = ((1 as u64) << 2*s) - 1;
    let kmask : u64 = ((1 as u64) << 2*k) - 1;
    let t = 3;
    let mut hash_count = 0;
    let mut seq_hashes = Vec::new();
    let mut pos_to_seq_coord = Vec::new();
    let mut qs = VecDeque::<u64>::new();
    let mut qs_pos = VecDeque::<usize>::new();
    let mut seq_len = seq.len();
    let mut qs_size = 0;
    let mut qs_min_val = u64::max_value();
    let mut qs_min_pos : i32 = -1;
    let mut l = 0;
    let mut xk : [u64; 2] = [0; 2];
    let mut xs : [u64; 2] = [0; 2];
    let mut kshift : u64 = (k as u64 - 1) * 2;
    let mut sshift : u64 = (s as u64 - 1) * 2;
    for i in 0..seq.len() {
        let mut c = seq_nt4_table[seq[i] as usize];
        if c < 4 {
            xk[0] = (xk[0] << 2 | c as u64) & kmask;                  // forward strand
            xk[1] = xk[1] >> 2 | ((3 - c) as u64) << kshift;  // reverse strand
            xs[0] = (xs[0] << 2 | c as u64) & smask;                  // forward strand
            xs[1] = xs[1] >> 2 | ((3 - c) as u64) << sshift;  // reverse strand
            l += 1;
            if l >= s {
                let mut ys : u64 = match xs[0] < xs[1]{
                    true => xs[0],
                    false => xs[1]
                };
                let hash_s = hash(ys, smask);
                if qs_size < k - s {
                    qs.push_back(hash_s);
                    qs_pos.push_back(i - s+ 1);
                    qs_size += 1;
                }
                else if qs_size == k - s {
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
                        let mut yk : u64 = match xk[0] < xk[1]{
                            true => xk[0],
                            false => xk[1]
                        };
                        let mut hash_k = hash(yk, kmask);
                        seq_hashes.push(hash_k);
                        pos_to_seq_coord.push(i - k + 1);
                        hash_count += 1;
                    }
                }
                else {
                    let mut new_minimizer = false;
                    let tuple = update_window(&mut qs, &mut qs_pos, qs_min_val, qs_min_pos, hash_s, i - s + 1, new_minimizer);
                    qs_min_val = tuple.0; qs_min_pos = tuple.1; new_minimizer = tuple.2;
                    if qs_min_pos == qs_pos[t-1] as i32 { // occurs at t:th position in k-mer
                        let mut yk : u64 = match xk[0] < xk[1] {
                            true => xk[0],
                            false => xk[1]
                        };
                        let mut hash_k = hash(yk, kmask);
                        seq_hashes.push(hash_k);
                        pos_to_seq_coord.push(i - k + 1);
                        hash_count += 1;
                    }
                }
            }
        } else {
            qs_min_val = u64::max_value();
            qs_min_pos = -1;
            l = 0; xs = [0; 2]; xk = [0; 2];
            qs_size = 0;
            qs.clear();
            qs_pos.clear();

        }
    }
    return (seq_hashes, pos_to_seq_coord)
} 

pub fn get_ksyncmer(i: usize, k: usize, num_hashes: usize, ref_id: u64, string_hashes: &Vec<u64>, pos_to_seq_coord: &Vec<usize>, reverse: bool) -> Option<(u64, u64, usize, usize, bool)> {
    if i + k - 1 < num_hashes {
        let mut ksyncmer_pos_end = pos_to_seq_coord[i];
        let mut ksyncmer_hash : u64 = 0;
        let mut hash_frac = 2;
        for j in i..i+k {
            let mut lmer_hash = string_hashes[j];
            ksyncmer_pos_end = pos_to_seq_coord[j];
            ksyncmer_hash += lmer_hash;
            hash_frac += 1;
        }
        let mut ksyncmer_pos_start = pos_to_seq_coord[i];
        let mut s = (ksyncmer_hash, ref_id, ksyncmer_pos_start, ksyncmer_pos_end, reverse);
        return Some(s);
    }
    else {return None;}
}

pub fn seq_to_ksyncmers_read(seq: &[u8], id: u64, params: &Params) -> (Vec<u64>, Vec<usize>, MersVectorRead) {
    let l = params.l;
    let k = params.k;
    let mut ksyncmers = MersVectorRead::new();
    let read_length = seq.len();
    let lmask : u64 = ((1 as u64) << 2*l) - 1;
    let mut string_hashes = Vec::<u64>::new();
    let mut pos_to_seq_coord = Vec::<usize>::new();
    let (mut string_hashes, mut pos_to_seq_coord) = extract_syncmers(seq, params);
    let num_hashes = string_hashes.len();
    for i in 0..num_hashes + 1 {
        let s = get_ksyncmer(i, k, num_hashes, id, &string_hashes, &pos_to_seq_coord, false);
        if s.is_none() {continue;}
        ksyncmers.push(s.unwrap());
    }
    string_hashes.reverse();
    pos_to_seq_coord.reverse();
    for i in 0..num_hashes {
        pos_to_seq_coord[i] = read_length - pos_to_seq_coord[i] - l;
    }
    for i in 0..num_hashes + 1 {
        let s = get_ksyncmer(i, k, num_hashes, id, &string_hashes, &pos_to_seq_coord, true);
        if s.is_none() {continue;}
        ksyncmers.push(s.unwrap());
    }
    return (string_hashes, pos_to_seq_coord, ksyncmers);
}

pub fn seq_to_ksyncmers_ref(seq: &[u8], id: u64, params: &Params) -> (Vec<u64>, Vec<usize>, MersVectorRead) {
    let l = params.l;
    let k = params.k;
    let mut ksyncmers = MersVectorRead::new();
    let read_length = seq.len();
    let lmask : u64 = ((1 as u64) << 2*l) - 1;
    let mut string_hashes = Vec::<u64>::new();
    let mut pos_to_seq_coord = Vec::<usize>::new();
    let (mut string_hashes, mut pos_to_seq_coord) = extract_syncmers(seq, params);
    let num_hashes = string_hashes.len();
    for i in 0..num_hashes + 1 {
        let s = get_ksyncmer(i, k, num_hashes, id, &string_hashes, &pos_to_seq_coord, false);
        if s.is_none() {continue;}
        ksyncmers.push(s.unwrap());
    }
    return (string_hashes, pos_to_seq_coord, ksyncmers);
}

pub fn remove_kmer_hash_from_flat_vector(flat_vector: &MersVectorRead) -> MersVectorReduced {
    let mut mers_vector_reduced = MersVectorReduced::new();
    for entry in flat_vector.iter() {
        let s = (entry.1, entry.2, entry.3);
        mers_vector_reduced.push(s);
    }
    return mers_vector_reduced;
}

pub fn index(ref_seq: &[u8], ref_id: u64, params: &Params, read: bool) -> (Vec<u64>, Vec<usize>, MersVectorRead) {
    if read {return seq_to_ksyncmers_read(ref_seq, ref_id, params);}
    else {return seq_to_ksyncmers_ref(ref_seq, ref_id, params);}
}

pub fn find_nams(query_mers: &MersVectorRead, ref_index: &PosIndex, mers_index: &KmerLookupMod, hash_counts: &DashMap<u64, usize>, l: usize, read: &[u8], filter_cutoff: usize, k: usize, id_hashes: &DashMap<u64, String>, id_lengths: &DashMap<u64, usize>) -> Vec<NAM> {
    let mut hits_per_ref = DashMap::<u64, Vec<Hit>>::new();
    let mut hit_count_reduced = 0;
    let mut hit_count_all = 0;
    let read_length = read.len();
    let mut i = 0;
    let mut prev_hit_end = 0;
    let mut prev_rc = false;
    while i < query_mers.len() {
        let q = query_mers[i];
        let mut curr_rc = q.4;
        let mut query_new_start = 0;
        if prev_hit_end != 0 && prev_rc == curr_rc {query_new_start = prev_hit_end + 1;}
        let mut h = Hit {query_id: q.1, ref_id: 0, query_s: query_new_start, query_e: q.3, ref_s: 0, ref_e: 0, hit_count: 0, is_rc: q.4};
        let mer_hashv = q.0;
        let mut mer_entry = mers_index.get_mut(&mer_hashv);
        let mut extend_offset = 0;
        if mer_entry.is_some() {
            h.hit_count = k;
            let mer = mer_entry.unwrap();
            let count = hash_counts.get(&q.0).unwrap();
            let ref_id = mer.0;
            let offset = mer.1;
            let ref_mers = ref_index.get(&ref_id).unwrap();
            if *count < filter_cutoff || filter_cutoff == 0 {
                let mut j = offset;
                let r = ref_mers[j];
                h.ref_s = r.2 - (q.2 - query_new_start);
                //println!("Query: {:?}\tRef: {:?}",  q, r);
                h.ref_id = ref_id;
                while query_mers[i+extend_offset].0 == ref_mers[j+extend_offset].0 {
                    let r_next = ref_mers[j+extend_offset];
                    let q_next = query_mers[i+extend_offset];
                    h.query_e = q_next.3 + l - 1;
                    h.ref_e = r_next.3 + l - 1;
                    h.hit_count += 1;
                    prev_rc = h.is_rc;
                    //println!("Extend: {}\tQuery next: {}\tRef next: {}\tQuery end: {}\tRef end: {}",  extend_offset, q_next.0, r_next.0, h.query_e, h.ref_e);
                    extend_offset += 1;
                    if i + extend_offset == query_mers.len() {break;}
                    if j + extend_offset == ref_mers.len() {break;}
                }
                let ref_entry = hits_per_ref.get_mut(&r.1);
                if ref_entry.is_some() {
                    ref_entry.unwrap().push(h.clone());
                }
                else {hits_per_ref.insert(ref_id, vec![h.clone()]);}
                hit_count_all += 1;
                prev_hit_end = h.query_e;
                    /*let r = ref_mers[j];
                    let ref_id = r.1;
                    let ref_s = r.2;
                    let ref_e = r.3 + l;
                    h.ref_id = ref_id;
                    h.ref_s = ref_s;
                    h.ref_e = ref_e;
                    h.hit_count = count;*/
                    
                //}
           }
        }
        if extend_offset == 0 {i += 1;} else {i += extend_offset;}
    }
    let mut final_nams = Vec::<NAM>::new();
    let mut prev_key : u64 = 0;
    for (key, val) in hits_per_ref.into_iter() {
        let ref_id = key;
        let hits = val;
        let mut hit_ids : Vec<u64> = hits.iter().map(|hit| hit.query_id).collect();
        hit_ids.sort_unstable();
        hit_ids.dedup();
        for id in hit_ids.iter() {
            let query_id_str = id_hashes.get(&id).unwrap();
            let query_len = *id_lengths.get(&id).unwrap();
            let ref_id_str = id_hashes.get(&ref_id).unwrap();
            println!("_____________________{} to {}_____________\nRC", *query_id_str, *ref_id_str);
            let mut prev_q_start = 0;
            let mut mono_rc = true;
            let mut fwd_hit_avg = 0.0;
            let mut rc_hit_avg = 0.0;
            let mut mono_fwd = true;
            let mut fwd_counts = 0;
            let mut rc_counts = 0;
            let mut rc_ref_end = 0;
            let mut rc_query_end = 0;
            let mut fwd_ref_end = 0;
            let mut fwd_query_end = 0;
            let mut rc_hits : Vec<&Hit> = hits.iter().filter(|hit| hit.is_rc == true && hit.query_id == *id).collect();
            let mut rc_hits_cleared = Vec::<&Hit>::new();
            //rc_hits.retain(|x| x.hit_count > 1);
            if rc_hits.len() != 0 {
                rc_hits.sort_by(|a, b| ((a.query_s)).cmp(&((b.query_s))));
                let mut prev_hit = rc_hits[0];
                let mut prev_ref_s = prev_hit.ref_s;
                let query_hit_id_str = id_hashes.get(&prev_hit.query_id).unwrap();
                let ref_hit_id_str = id_hashes.get(&prev_hit.ref_id).unwrap();
                rc_hits_cleared.push(prev_hit);
                println!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", *query_hit_id_str, query_len, prev_hit.query_s, prev_hit.query_e, *ref_hit_id_str, prev_hit.ref_s, prev_hit.ref_e, prev_hit.hit_count, prev_hit.is_rc);
                for i in 1..rc_hits.len() {
                    let mut curr_hit = rc_hits[i];
                    let mut curr_ref_s = curr_hit.ref_s;
                    let query_hit_len = id_lengths.get(&curr_hit.query_id).unwrap();
                    let query_hit_id_str = id_hashes.get(&curr_hit.query_id).unwrap();
                    let ref_hit_id_str = id_hashes.get(&curr_hit.ref_id).unwrap();
                    if curr_ref_s < prev_ref_s || curr_ref_s > prev_ref_s + *query_hit_len {
                        continue;
                    }
                    else {
                        println!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", *query_hit_id_str, query_len, curr_hit.query_s, curr_hit.query_e, *ref_hit_id_str, curr_hit.ref_s, curr_hit.ref_e, curr_hit.hit_count, curr_hit.is_rc);
                        rc_hits_cleared.push(curr_hit);
                        prev_ref_s = curr_ref_s;
                    }
                }
                rc_counts = rc_hits_cleared.iter().map(|hit| hit.hit_count).sum();
                if rc_hits_cleared.len() != 0 {rc_hit_avg = rc_counts as f64 / rc_hits_cleared.len() as f64;}
                let mut rc_end = rc_hits_cleared.iter().max_by(|a, b| ((a.query_e)).cmp(&((b.query_e)))).unwrap();
                rc_ref_end = rc_end.ref_e; rc_query_end = rc_end.query_e;
            }
            let mut fwd_hits : Vec<&Hit> = hits.iter().filter(|hit| hit.is_rc == false && hit.query_id == *id).collect();
            let mut fwd_hits_cleared = Vec::<&Hit>::new();
            //fwd_hits.retain(|x| x.hit_count > 1);
            println!("FW");
            if fwd_hits.len() != 0 {
                fwd_hits.sort_by(|a, b| ((a.query_s)).cmp(&((b.query_s))));
                let mut prev_hit = fwd_hits[0];
                let mut prev_ref_s = fwd_hits[0].ref_s;
                let query_hit_id_str = id_hashes.get(&prev_hit.query_id).unwrap();
                let ref_hit_id_str = id_hashes.get(&prev_hit.ref_id).unwrap();
                fwd_hits_cleared.push(prev_hit);
                println!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", *query_hit_id_str, query_len, prev_hit.query_s, prev_hit.query_e, *ref_hit_id_str, prev_hit.ref_s, prev_hit.ref_e, prev_hit.hit_count, prev_hit.is_rc);
                for i in 1..fwd_hits.len() {
                    let mut curr_hit = fwd_hits[i];
                    let mut curr_ref_s = curr_hit.ref_s;
                    let query_hit_len = id_lengths.get(&curr_hit.query_id).unwrap();
                    let query_hit_id_str = id_hashes.get(&curr_hit.query_id).unwrap();
                    let ref_hit_id_str = id_hashes.get(&curr_hit.ref_id).unwrap();
                    if curr_ref_s < prev_ref_s || curr_ref_s > prev_ref_s + *query_hit_len {
                        continue;
                    }
                    else {
                        println!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", *query_hit_id_str, query_len, curr_hit.query_s, curr_hit.query_e, *ref_hit_id_str, curr_hit.ref_s, curr_hit.ref_e, curr_hit.hit_count, curr_hit.is_rc);
                        fwd_hits_cleared.push(curr_hit);
                        prev_ref_s = curr_ref_s;
                    }
                }
                let mut fwd_end = fwd_hits_cleared.iter().max_by(|a, b| ((a.query_e)).cmp(&((b.query_e)))).unwrap();
                fwd_ref_end = fwd_end.ref_e; fwd_query_end = fwd_end.query_e;
                fwd_counts = fwd_hits_cleared.iter().map(|hit| hit.hit_count).sum();
                if fwd_hits_cleared.len() != 0 {fwd_hit_avg = fwd_counts as f64 / fwd_hits_cleared.len() as f64;}
                    if fwd_hit_avg > rc_hit_avg {
                        let mut n = NAM {query_id: *id, ref_id: ref_id, query_s: 0, query_e: query_len - 1, query_last_hit_pos: 0, ref_s: fwd_hits_cleared[0].ref_s , ref_e: fwd_hits_cleared[0].ref_s  + query_len - 1, ref_last_hit_pos: 0, n_hits: fwd_counts, previous_query_start: 0, previous_ref_start: 0, is_rc: false};
                        //println!("{:?}", n);
                        
                        final_nams.push(n.clone());
                    }
                    else {
                        let mut n = NAM {query_id: *id, ref_id: ref_id, query_s: 0, query_e: query_len - 1, query_last_hit_pos: 0, ref_s: rc_hits_cleared[0].ref_s, ref_e: rc_hits_cleared[0].ref_s  + query_len - 1, ref_last_hit_pos: 0, n_hits: rc_counts, previous_query_start: 0, previous_ref_start: 0, is_rc: true};
                        //println!("{:?}", n);
                        final_nams.push(n.clone());
                    }
            }
           // println!("_____________________");
        }     
    }
    //println!("NAMs found:\t{}", final_nams.len());
    return final_nams;
}

pub fn output_paf(all_nams: &mut Vec<(Vec<NAM>, u64)>, l: usize, id_hashes: ReadOnlyView<u64, String>, id_lengths: ReadOnlyView<u64, usize>, paf_file: &mut File, k: usize) {
    let mut mapped_ref_nams = DashMap::<(u64, usize, usize), Vec<NAM>>::new();
    for (nams, id) in all_nams.iter_mut() {
        for nam in nams.iter() {
            println!("{}\t{:?}", id, nam);
        }
        if nams.len() == 0 {continue;}
        nams.retain(|x| x.n_hits > k);
        let max_nam = nams.iter().max_by(|a, b| ((a.query_e - a.query_s)).cmp(&((b.query_e - b.query_s))));
        if max_nam.is_none() {continue;}
        let n = max_nam.unwrap();
        let ref_entry = mapped_ref_nams.get_mut(&(n.ref_id, n.ref_s, n.ref_e));
        if ref_entry.is_some() {
            ref_entry.unwrap().push(n.clone());
        }
        else {mapped_ref_nams.insert((n.ref_id, n.ref_s, n.ref_e), vec![n.clone()]);} 
        
    }
    for mapped_nams in mapped_ref_nams.iter() {
        if mapped_nams.len() != 1 {
           /* for nam in mapped_nams.iter() {
                println!("Ref:\t{:?}\t{:?}", mapped_nams.key(), nam);
            }*/
        }
        let max_hit_nam = mapped_nams.iter().max_by(|a, b| a.n_hits.cmp(&(b.n_hits)));
        if max_hit_nam.is_none() {continue;}
        let n = max_hit_nam.unwrap();
        //let n = nams.iter().max_by(|a, b| (a.n_hits).cmp(&(b.n_hits))).unwrap();
        let mut o = String::new();
        if n.is_rc {o = "-".to_string();}
        else {o = "+".to_string();}
        let query_id = id_hashes.get(&n.query_id).unwrap();
        let ref_id = id_hashes.get(&n.ref_id).unwrap();
        let query_len = id_lengths.get(&n.query_id).unwrap();
        let ref_len = id_lengths.get(&n.ref_id).unwrap();
        let paf_line = format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n", query_id, query_len, n.query_s, n.query_e, o, ref_id, ref_len, n.ref_s, n.ref_e, n.n_hits, n.ref_e - n.ref_s + 1, "255");
        write!(paf_file, "{}", paf_line).expect("Error writing line.");
    }
}

