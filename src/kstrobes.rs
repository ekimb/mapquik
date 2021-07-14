use std::collections::VecDeque;
use crate::DashMap;
use crate::File;
use dashmap::ReadOnlyView;
use std::io::Write;
use std::cmp::Ordering;

use super::Params;
use super::{MersVector, MersVectorReduced, MersVectorRead, PosIndex, KmerLookup};

#[derive(Clone, Debug)]
pub struct Hit {
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
    let s = params.syncmer;
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

pub fn extend_kstrobe(string_hashes: &Vec<u64>, strobe_hashval: u64, w_start: usize, w_end: usize, q: u64) -> (i32, u64) {
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

/*get_kstrobe
while iter != k {
    let res = get_kstrobe(k, i, wmin, wmax, num_hashes, ref_id, string_hashes, pos_to_seq_coord, q, reverse)
    while res.is_some() {
        res = get_kstrobe
    }
    let mut strobe_pos_next;
    let mut strobe_hashval_next;
    if i + wmax < num_hashes {
        let mut w_start = i + wmin;
        let mut w_end = i + wmax;
        let mut strobe_hash = string_hashes[i];
        let tuple = get_next_strobe(string_hashes, strobe_hash, w_start, w_end, q);
        strobe_pos_next = tuple.0; strobe_hashval_next = tuple.1;
    }
    else if (i + wmin + 1 < num_hashes) && (num_hashes <= i + wmax) {
        let mut w_start = i + wmin;
        let mut w_end = num_hashes - 1;
        let mut strobe_hash = string_hashes[i];
        let tuple = get_next_strobe(string_hashes, strobe_hash, w_start, w_end, q);
        strobe_pos_next = tuple.0; strobe_hashval_next = tuple.1;
    }
    else {return None;}

    let mut hash_randstrobe = string_hashes[i]/2 + strobe_hashval_next/3;
    let mut seq_pos_strobe1 = pos_to_seq_coord[i];
    let mut seq_pos_strobe2 = pos_to_seq_coord[strobe_pos_next as usize];
    let mut s = (hash_randstrobe, ref_id, seq_pos_strobe1, seq_pos_strobe2, reverse);
}
*/
pub fn get_kstrobe(k: usize, i: usize, wmin: usize, wmax: usize, num_hashes: usize, ref_id: u64, string_hashes: &Vec<u64>, pos_to_seq_coord: &Vec<usize>, q: u64, reverse: bool) -> Option<(u64, u64, usize, usize, bool)> {
    let mut iterations = 1;
    let mut kstrobe_hashval = 0;
    let res = get_kstrobe_util(k, 1, i, wmin, wmax, num_hashes, ref_id, string_hashes, pos_to_seq_coord, q, reverse);
    if res.is_none() {return None;}
    else {
        let (k_fin, strobe_pos_next, strobe_hashval_next) = res.unwrap();
        let mut hash_randstrobe = string_hashes[i]/2 + strobe_hashval_next/3;
        let mut seq_pos_strobe1 = pos_to_seq_coord[i];
        let mut seq_pos_strobe2 = pos_to_seq_coord[strobe_pos_next as usize];
        let mut s = (hash_randstrobe, ref_id, seq_pos_strobe1, seq_pos_strobe2, reverse);
        if k != k_fin {return None;}
        else{return Some(s);}
    }
}

pub fn get_kstrobe_util(k_bound: usize, curr_k: usize, i: usize, wmin: usize, wmax: usize, num_hashes: usize, ref_id: u64, string_hashes: &Vec<u64>, pos_to_seq_coord: &Vec<usize>, q: u64, reverse: bool) -> Option<(usize, usize, u64)> {
        if i + wmax < num_hashes {
            let mut strobe_hash = string_hashes[i];
            let mut w_start = i + wmin;
            let mut w_end = i + wmax;
            let tuple = extend_kstrobe(string_hashes, strobe_hash, w_start, w_end, q);
            let mut k = curr_k + 1;
            if k == k_bound {return Some((k,tuple.0 as usize, tuple.1));}
            return get_kstrobe_util(k_bound, k, tuple.0 as usize, wmin, wmax, num_hashes, ref_id, string_hashes, pos_to_seq_coord, q, reverse)
        }
        else if (i + wmin + 1 < num_hashes) && (num_hashes <= i + wmax) {
            let mut strobe_hash = string_hashes[i];
            let mut w_start = i + wmin;
            let mut w_end = num_hashes - 1;
            let tuple = extend_kstrobe(string_hashes, strobe_hash, w_start, w_end, q);
            let mut k = curr_k + 1;
            if k == k_bound {return Some((k,tuple.0 as usize, tuple.1));}
            return get_kstrobe_util(k_bound, k, tuple.0 as usize, wmin, wmax, num_hashes, ref_id, string_hashes, pos_to_seq_coord, q, reverse)
        }
        else {return None;}
}

pub fn seq_to_kstrobes_read(seq: &[u8], id: u64, params: &Params) -> MersVectorRead {
    let k = params.l;
    let s = params.syncmer;
    let wmin = params.wmin;
    let wmax = params.wmax;
    let mut kstrobes = MersVectorRead::new();
    let read_length = seq.len();
    if read_length < wmax {return kstrobes;}
    let kmask : u64 = ((1 as u64) << 2*k) - 1;
    let smask : u64 = ((1 as u64) << 2*s) - 1;
    let q = 2_u64.pow(16) - 1;
    let mut string_hashes = Vec::<u64>::new();
    let mut pos_to_seq_coord = Vec::<usize>::new();
    let t = 3;
    let (mut string_hashes, mut pos_to_seq_coord) = extract_syncmers(seq, params);
    let num_hashes = string_hashes.len();
    for i in 0..num_hashes + 1 {
        let s = get_kstrobe(params.k, i, wmin, wmax, num_hashes, id, &string_hashes, &pos_to_seq_coord, q, false);
        if s.is_none() {continue;}
        kstrobes.push(s.unwrap());
    }
    string_hashes.reverse();
    pos_to_seq_coord.reverse();
    for i in 0..num_hashes {
        pos_to_seq_coord[i] = read_length - pos_to_seq_coord[i] - k;
    }
    for i in 0..num_hashes + 1 {
        let s = get_kstrobe(params.k, i, wmin, wmax, num_hashes, id, &string_hashes, &pos_to_seq_coord, q, true);
        if s.is_none() {continue;}
        kstrobes.push(s.unwrap());
    }
    return kstrobes;
}

pub fn seq_to_kstrobes_ref(seq: &[u8], id: u64, params: &Params) -> MersVectorRead {
    let k = params.l;
    let s = params.syncmer;
    let wmin = params.wmin;
    let wmax = params.wmax;
    let mut kstrobes = MersVectorRead::new();
    let read_length = seq.len();
    if read_length < wmax {return kstrobes;}
    let kmask : u64 = ((1 as u64) << 2*k) - 1;
    let smask : u64 = ((1 as u64) << 2*s) - 1;
    let q = 2_u64.pow(16) - 1;
    let mut string_hashes = Vec::<u64>::new();
    let mut pos_to_seq_coord = Vec::<usize>::new();
    let t = 3;
    let (mut string_hashes, mut pos_to_seq_coord) = extract_syncmers(seq, params);
    let num_hashes = string_hashes.len();
    for i in 0..num_hashes + 1 {
        let s = get_kstrobe(params.k, i, wmin, wmax, num_hashes, id, &string_hashes, &pos_to_seq_coord, q, false);
        if s.is_none() {continue;}
        kstrobes.push(s.unwrap());
    }
    return kstrobes;
}

pub fn remove_kmer_hash_from_flat_vector(flat_vector: &MersVectorRead) -> MersVectorReduced {
    let mut mers_vector_reduced = MersVectorReduced::new();
    for entry in flat_vector.iter() {
        let s = (entry.1, entry.2, entry.3);
        mers_vector_reduced.push(s);
    }
    return mers_vector_reduced;
}

pub fn index(ref_seq: &[u8], ref_id: u64, params: &Params, read: bool) -> MersVectorRead {
    if read {return seq_to_kstrobes_read(ref_seq, ref_id, params);}
    else {return seq_to_kstrobes_ref(ref_seq, ref_id, params);}
}

pub fn find_nams(query_mers: &MersVectorRead, ref_mers: &MersVectorReduced, mers_index: &KmerLookup, k: usize, read: &[u8], hit_upper_window_lim: usize, filter_cutoff: usize) -> Vec<NAM> {
    let mut hits_per_ref = DashMap::<u64, Vec<Hit>>::new();
    let mut hit_count_reduced = 0;
    let mut hit_count_all = 0;
    let read_length = read.len();
    for q in query_mers.iter() {
        let mut h = Hit {ref_id: 0, query_s: q.2, query_e: q.3 + k, ref_s: 0, ref_e: 0, hit_count: 0, is_rc: q.4};
        let mer_hashv = q.0;
        let mut mer_entry = mers_index.get_mut(&mer_hashv);
        if mer_entry.is_some() {
            let mer = mer_entry.unwrap();
            let offset = mer.0;
            let count = mer.1;
            if count <= filter_cutoff || filter_cutoff == 0 {
                for j in offset..offset+count {
                    let r = ref_mers[j];
                    let ref_id = r.0;
                    let ref_s = r.1;
                    let ref_e = r.2 + k;
                    h.ref_id = ref_id;
                    h.ref_s = ref_s;
                    h.ref_e = ref_e;
                    h.hit_count = count;
                    let ref_entry = hits_per_ref.get_mut(&ref_id);
                    if ref_entry.is_some() {
                        ref_entry.unwrap().push(h.clone());
                    }
                    else {hits_per_ref.insert(ref_id, vec![h.clone()]);}
                    hit_count_all += 1;
                }
           }
        }
    }
    let mut open_nams = Vec::<NAM>::new();
    let mut final_nams = Vec::<NAM>::new();
    for (key, val) in hits_per_ref.into_iter() {
        let ref_id = key;
        let hits = val;
        let mut prev_q_start = 0;
        for h in hits.iter() {
            let mut is_added = false;
            for o in open_nams.iter_mut() {
                if (o.is_rc == h.is_rc) && (o.previous_query_start < h.query_s) && (h.query_s <= o.query_e) && (o.previous_ref_start < h.ref_s) && (h.ref_s <= o.ref_e) {
                    if (h.query_e > o.query_e) && (h.ref_e > o.ref_e) {
                        o.query_e = h.query_e;
                        o.ref_e = h.ref_e;
                        o.previous_query_start = h.query_s;
                        o.previous_ref_start = h.ref_s; // keeping track so that we don't . Can be caused by interleaved repeats.
                        o.query_last_hit_pos = h.query_s; // log the last strobemer hit in case of outputting paf
                        o.ref_last_hit_pos = h.ref_s; // log the last strobemer hit in case of outputting paf
                        o.n_hits += 1;
//                        o.score += (float)1/ (float)h.hit_count;
                        is_added = true;
                        break;
                    }  
                    else if (h.query_e <= o.query_e) && (h.ref_e <= o.ref_e) {
                        o.previous_query_start = h.query_s;
                        o.previous_ref_start = h.ref_s; // keeping track so that we don't . Can be caused by interleaved repeats.
                        o.query_last_hit_pos = h.query_s; // log the last strobemer hit in case of outputting paf
                        o.ref_last_hit_pos = h.ref_s; // log the last strobemer hit in case of outputting paf
                        o.n_hits += 1;
//                        o.score += (float)1/ (float)h.hit_count;
                        is_added = true;
                        break;
                    }
                }
            }
            if !is_added {
                let mut n = NAM {ref_id: ref_id, query_s: h.query_s, query_e: h.query_e, query_last_hit_pos: h.query_s, ref_s: h.ref_s, ref_e: h.ref_e, ref_last_hit_pos: h.ref_s, n_hits: 1, previous_query_start: h.query_s, previous_ref_start: h.ref_s, is_rc: h.is_rc};
//                n.score += (float)1/ (float)h.hit_count;
                open_nams.push(n);
            }
            if h.query_s > prev_q_start + k {
                // Output all NAMs from open_matches to final_nams that the current hit have passed
                for n in open_nams.iter() {
                    if n.query_e < h.query_s {final_nams.push(n.clone());}
                }
                let c = h.query_s;
                open_nams.retain(|nam| nam.query_e >= c);
                prev_q_start = h.query_s;
            }

        }
        for n in open_nams.iter() {final_nams.push(n.clone());}
    }
    //println!("NAMs found:\t{}", final_nams.len());
    return final_nams;
}

pub fn output_paf(all_nams: &mut Vec<(Vec<NAM>, u64)>, k: usize, id_hashes: ReadOnlyView<u64, String>, id_lengths: ReadOnlyView<u64, usize>, paf_file: &mut File) {
    for (nams, id) in all_nams.iter_mut() {
        if nams.len() == 0 {continue;}
        let n = nams.iter().max_by(|a, b| (a.n_hits * (a.query_e - a.query_s)).cmp(&(b.n_hits * (b.query_e - b.query_s)))).unwrap();
        let mut o = String::new();
        if n.is_rc {o = "-".to_string();}
        else {o = "+".to_string();}
        let query_id = id_hashes.get(&id).unwrap();
        let ref_id = id_hashes.get(&n.ref_id).unwrap();
        let query_len = id_lengths.get(&id).unwrap();
        let ref_len = id_lengths.get(&n.ref_id).unwrap();
        let paf_line = format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n", query_id, query_len, n.query_s, n.query_last_hit_pos + k, o, ref_id, ref_len, n.ref_s, n.ref_last_hit_pos + k, n.n_hits, n.ref_last_hit_pos + k - n.ref_s, "255");
        write!(paf_file, "{}", paf_line).expect("Error writing line.");
    }
}

