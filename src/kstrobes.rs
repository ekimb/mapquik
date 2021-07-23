use std::collections::VecDeque;
use crate::DashMap;
use crate::File;
use crate::hash_id;
use dashmap::ReadOnlyView;
use std::io::Write;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::collections::HashSet;
use super::Params;
pub type Match = (String, String, usize, usize, usize, usize, usize, usize, f64, bool);

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
    offset: usize,
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

/*pub fn extend_kstrobe(string_hashes: &Vec<u64>, strobe_hashval: u64, w_start: usize, w_end: usize, q: u64) -> (i32, u64) {
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
}*/

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

pub fn get_randstrobe(i: usize, wmin: usize, wmax: usize, num_hashes: usize, ref_id: &str, string_hashes: &Vec<u64>, pos_to_seq_coord: &Vec<usize>, q: u64, reverse: bool) -> Option<(u64, usize, usize, bool)> {
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
    let mut s = (hash_randstrobe, seq_pos_strobe1, seq_pos_strobe2, reverse);
    return Some(s);
}

pub fn seq_to_kstrobes_read(seq: &[u8], id: &str, params: &Params) -> (usize, Vec<(u64, usize, usize, bool, usize)>) {
    let k = params.k;
    let l = params.l;
    let s = params.s;
    let wmin = params.wmin;
    let wmax = params.wmax;
    let mut kstrobes = Vec::<(u64, usize, usize, bool, usize)>::new();
    let read_length = seq.len();
    if read_length < wmax {return (seq.len(), kstrobes);}
    let kmask : u64 = ((1 as u64) << 2*l) - 1;
    let smask : u64 = ((1 as u64) << 2*s) - 1;
    let q = 2_u64.pow(16) - 1;
    let mut string_hashes = Vec::<u64>::new();
    let mut pos_to_seq_coord = Vec::<usize>::new();
    let t = 3;
    let (mut string_hashes, mut pos_to_seq_coord) = extract_syncmers(seq, params);
    let num_hashes = string_hashes.len();
    for i in 0..num_hashes + 1 {
        let mut kstrobe_hash = Vec::<u64>::new();
        let mut kstrobe_id : u64 = 0;
        let mut kstrobe_start = 0;
        let mut kstrobe_end = 0;
        let mut kstrobe_reverse = false;
        for j in i..i+k {
            let mut s = get_randstrobe(j, wmin, wmax, num_hashes, id, &string_hashes, &pos_to_seq_coord, q, false);
            if s.is_none() {break;}
            let s_unw = s.unwrap();
            kstrobe_hash.push(s_unw.0);
            kstrobe_reverse = s_unw.3;
            if j == i {kstrobe_start = s_unw.1;}
            if j == i + k - 1 {kstrobe_end = s_unw.2;}
        }
        let mut kstrobe_u64_hash = hash_id(&kstrobe_hash);
        if kstrobe_end != 0 {kstrobes.push((kstrobe_u64_hash, kstrobe_start, kstrobe_end, kstrobe_reverse, i));}
    }
    string_hashes.reverse();
    pos_to_seq_coord.reverse();
    for i in 0..num_hashes {
        pos_to_seq_coord[i] = read_length - pos_to_seq_coord[i] - k;
    }
    for i in 0..num_hashes + 1 {
        let mut kstrobe_hash = Vec::<u64>::new();
        let mut kstrobe_id : u64 = 0;
        let mut kstrobe_start = 0;
        let mut kstrobe_end = 0;
        let mut kstrobe_reverse = true;
        for j in i..i+k {
            let mut s = get_randstrobe(j, wmin, wmax, num_hashes, id, &string_hashes, &pos_to_seq_coord, q, true);
            if s.is_none() {break;}
            let s_unw = s.unwrap();
            kstrobe_hash.push(s_unw.0);
            kstrobe_reverse = s_unw.3;
            if j == i {kstrobe_start = s_unw.1;}
            if j == i + k - 1 {kstrobe_end = s_unw.2;}
        }
        let mut kstrobe_u64_hash = hash_id(&kstrobe_hash);
        if kstrobe_end != 0 {kstrobes.push((kstrobe_u64_hash, kstrobe_start, kstrobe_end, kstrobe_reverse, i));}
    }
    return (seq.len(), kstrobes);
}

pub fn seq_to_kstrobes_ref(seq: &[u8], id: &str, params: &Params) -> (usize, Vec<(u64, usize, usize, bool, usize)>) {
    let k = params.k;
    let l = params.l;
    let s = params.s;
    let wmin = params.wmin;
    let wmax = params.wmax;
    let mut kstrobes = Vec::<(u64, usize, usize, bool, usize)>::new();
    let read_length = seq.len();
    if read_length < wmax {return (seq.len(), kstrobes);}
    let lmask : u64 = ((1 as u64) << 2*l) - 1;
    let smask : u64 = ((1 as u64) << 2*s) - 1;
    let q = 2_u64.pow(16) - 1;
    let mut string_hashes = Vec::<u64>::new();
    let mut pos_to_seq_coord = Vec::<usize>::new();
    let t = 3;
    let (mut string_hashes, mut pos_to_seq_coord) = extract_syncmers(seq, params);
    let num_hashes = string_hashes.len();
    for i in 0..num_hashes + 1 {
        let mut kstrobe_hash = Vec::<u64>::new();
        let mut kstrobe_id : u64 = 0;
        let mut kstrobe_start = 0;
        let mut kstrobe_end = 0;
        let mut kstrobe_reverse = false;
        for j in i..i+k {
            let mut s = get_randstrobe(j, wmin, wmax, num_hashes, id, &string_hashes, &pos_to_seq_coord, q, false);
            if s.is_none() {break;}
            let s_unw = s.unwrap();
            kstrobe_hash.push(s_unw.0);
            kstrobe_reverse = s_unw.3;
            if j == i {kstrobe_start = s_unw.1;}
            if j == i + k - 1 {kstrobe_end = s_unw.2;}
        }
        let mut kstrobe_u64_hash = hash_id(&kstrobe_hash);
        if kstrobe_end != 0 {kstrobes.push((kstrobe_u64_hash, kstrobe_start, kstrobe_end, kstrobe_reverse, i));}
    }
    return (seq.len(), kstrobes);
}

/*pub fn get_kstrobe(k: usize, i: usize, wmin: usize, wmax: usize, num_hashes: usize, ref_id: u64, string_hashes: &Vec<u64>, pos_to_seq_coord: &Vec<usize>, q: u64, reverse: bool) -> Option<(u64, u64, usize, usize, bool)> {
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
}*/

/*pub fn get_kstrobe_util(k_bound: usize, curr_k: usize, i: usize, wmin: usize, wmax: usize, num_hashes: usize, ref_id: u64, string_hashes: &Vec<u64>, pos_to_seq_coord: &Vec<usize>, q: u64, reverse: bool) -> Option<(usize, usize, u64)> {
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
    let s = params.s;
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
    let s = params.s;
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
}*/

/*pub fn remove_kmer_hash_from_flat_vector(flat_vector: &MersVectorRead) -> MersVectorReduced {
    let mut mers_vector_reduced = MersVectorReduced::new();
    for entry in flat_vector.iter() {
        let s = (entry.1, entry.2, entry.3);
        mers_vector_reduced.push(s);
    }
    return mers_vector_reduced;
}*/

pub fn index(ref_seq: &[u8], ref_id: &str, params: &Params, read: bool) -> (usize, Vec<(u64, usize, usize, bool, usize)>) {
    if read {return seq_to_kstrobes_read(ref_seq, ref_id, params);}
    else {return seq_to_kstrobes_ref(ref_seq, ref_id, params);}
}

/*pub fn find_nams(query_mers: &MersVectorRead, ref_mers: &MersVectorReduced, mers_index: &KmerLookup, k: usize, read: &[u8], hit_upper_window_lim: usize, filter_cutoff: usize) -> Vec<NAM> {
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
}*/

pub fn find_hits(query_id: &str, query_len: usize, query_mers: &Vec<(u64, usize, usize, bool, usize)>, ref_index: &DashMap<String, (usize, Vec<(u64, usize, usize, bool, usize)>)>, mers_index: &DashMap<u64, Vec<(String, usize)>>, l: usize, k: usize) -> Vec<Match> {
    let mut hits_per_ref = HashMap::<String, Vec<Hit>>::new();
    let mut hit_count_reduced = 0;
    let mut hit_count_all = 0;
    let mut i = 0;
    let mut prev_hit_end = 0;
    let mut prev_rc = 2;
    while i < query_mers.len() {
        let q = query_mers[i];
        let mut curr_rc = 0;
        if q.3 == true {curr_rc = 1;}
        /*let mut query_new_start = 0;
        if prev_rc == curr_rc {query_new_start = prev_hit_end + 1;}
        else {prev_hit_end = 0;}*/
        let mut h = Hit {query_id: query_id.to_string(), ref_id: String::new(), query_s: q.1, query_e: q.2, ref_s: 0, ref_e: 0, hit_count: 0, is_rc: q.3, query_span: 0, ref_span: 0, offset: 0};
        let mer_hashv = q.0;
        let mut mer_entry = mers_index.get(&mer_hashv);
        let mut extend_offset = 0;
        if mer_entry.is_some() {
            let mer_vec = mer_entry.unwrap();
            let mer = &mer_vec[0];
            let ref_id = mer.0.to_string();
            let offset = mer.1;
            let ref_mers_len = ref_index.get(&ref_id).unwrap();
            let ref_len = ref_mers_len.0;
            let ref_mers = &ref_mers_len.1;
            let mut j = offset;
            let r = ref_mers[j];
            h.offset = j;
            //println!("Query RC {}\tref RC {}", q.3, r.3);
            //println!("Query: {:?}\tRef: {:?}",  q, r);
            h.ref_id = ref_id.to_string();
            h.ref_s = r.1;
            while query_mers[i+extend_offset].0 == ref_mers[j+extend_offset].0 {
                let r_next = ref_mers[j+extend_offset];
                let q_next = query_mers[i+extend_offset];
                h.query_e = q_next.2;
                h.ref_e = r_next.2;
                h.hit_count += 1;
                h.query_span = h.query_e - h.query_s + 1; 
                h.ref_span = h.ref_e - h.ref_s + 1; 
                if h.is_rc == true {prev_rc = 1;} else {prev_rc = 0;}
                //println!("Extend: {}\tQuery next: {}\tRef next: {}\tQuery end: {}\tRef end: {}",  extend_offset, q_next.0, r_next.0, h.query_e, h.ref_e);
                extend_offset += 1;
                if i + extend_offset == query_mers.len() {break;}
                if j + extend_offset == ref_mers.len() {break;}
            }
            hits_per_ref.entry(ref_id).or_insert(vec![h.clone()]).push(h.clone());
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
        if extend_offset == 0 {i += 1;} else {i += extend_offset;}
    }
    let mut final_matches = Vec::<Match>::new();
    let mut prev_key : u64 = 0;
    for (key, val) in hits_per_ref.into_iter() {
        let ref_id = key;
        let hits = val;
        let mut ref_len = ref_index.get(&ref_id).unwrap().0;
        let mut hit_ids : HashSet<String> = hits.iter().map(|hit| hit.query_id.to_string()).collect();
        for id in hit_ids.iter() {
            //println!("_____________________{} to {}_____________\nRC", query_id, ref_id);
            let mut prev_q_start = 0;
            let mut rc_mono = true;
            let mut fwd_hit_avg = 0.0;
            let mut rc_hit_avg = 0.0;
            let mut fwd_mono = true;
            let mut write_fwd = false;
            let mut write_rc = false;
            let mut fwd_counts = 0;
            let mut rc_counts = 0;
            let mut rc_ref_end = 0;
            let mut rc_query_end = 0;
            let mut fwd_ref_end = 0;
            let mut fwd_query_end = 0;
            let mut rc_hits : Vec<&Hit> = hits.iter().filter(|hit| hit.is_rc == true && hit.query_id == query_id).collect();
            let mut rc_hits_cleared = Vec::<&Hit>::new();
            let mut fwd_hits : Vec<&Hit> = hits.iter().filter(|hit| hit.is_rc == false && hit.query_id == query_id).collect();
            let mut fwd_hits_cleared = Vec::<&Hit>::new();
            let mut rc_mappings = Vec::<Vec::<&Hit>>::new();
            let mut fwd_mappings = Vec::<Vec::<&Hit>>::new();

            // -----------------------FINDING RC MAPPINGS-----------------------
            if rc_hits.len() != 0 {
                let mut prev_rc_hit = rc_hits[0];
                let mut current_mapping = vec![prev_rc_hit];
                for i in 1..rc_hits.len() {
                    let curr_rc_hit = rc_hits[i];
                    if (curr_rc_hit.ref_s < prev_rc_hit.ref_s || /*curr_rc_hit.hit_count == 1 ||*/ (prev_rc_hit.ref_e  + query_len < curr_rc_hit.ref_s)) {
                        /*if current_mapping.len() != 0 {
                            rc_mappings.push(current_mapping);
                        }*/
                        current_mapping = Vec::<&Hit>::new();
                        continue;
                    } 
                    current_mapping.push(curr_rc_hit);
                    prev_rc_hit = curr_rc_hit;               
                }
                rc_mappings.push(current_mapping.to_vec());
            }
            // -----------------------FINDING FWD MAPPINGS-----------------------
            if fwd_hits.len() != 0 {
                let mut prev_fwd_hit = fwd_hits[0];
                let mut current_mapping = vec![prev_fwd_hit];
                for i in 1..fwd_hits.len() {
                    let curr_fwd_hit = fwd_hits[i];
                    if (curr_fwd_hit.ref_s < prev_fwd_hit.ref_s || /* curr_fwd_hit.hit_count == 1  ||*/ (prev_fwd_hit.ref_e  + query_len < curr_fwd_hit.ref_s)) {
                        /*if current_mapping.len() != 0 {
                            fwd_mappings.push(current_mapping);
                        }*/
                        current_mapping = Vec::<&Hit>::new();
                        continue;
                    }
                    current_mapping.push(curr_fwd_hit);
                    prev_fwd_hit = curr_fwd_hit;
                }
                fwd_mappings.push(current_mapping.to_vec());
            }
           
            let mut idf = HashMap::<u64, usize>::new();
            // -----------------------RECORDING RC IDF-----------------------
            for i in 0..rc_mappings.len() {
                let mapping = &rc_mappings[i];
                for hit in mapping.iter() {
                    let offset = hit.offset;
                    let ref_mers_len = ref_index.get(&ref_id).unwrap();
                    let ref_mers = &ref_mers_len.1;
                    for j in offset..offset+hit.hit_count {
                        *idf.entry(ref_mers[j].0).or_insert(1) += 1;
                    }
                }
            }
            // -----------------------RECORDING FWD IDF-----------------------
            for i in 0..fwd_mappings.len() {
                let mapping = &fwd_mappings[i];
                for hit in mapping.iter() {
                    let offset = hit.offset;
                    let ref_mers_len = ref_index.get(&ref_id).unwrap();
                    let ref_mers = &ref_mers_len.1;
                    for j in offset..offset+hit.hit_count {
                        *idf.entry(ref_mers[j].0).or_insert(1) += 1;
                    }
                }
            }

            let tot = fwd_mappings.len() + rc_mappings.len();
            // -----------------------CALCULATING FWD TF-IDF-----------------------
            let mut max_fwd_mapping = 0;
            let mut max_fwd_tf_idf = 0.0;
            for i in 0..fwd_mappings.len() {
                let mapping = &fwd_mappings[i];
                let mut mapping_tf_idf = 0.0;
                for hit in mapping.iter() {
                    let offset = hit.offset;
                    let ref_mers_len = ref_index.get(&ref_id).unwrap();
                    let ref_mers = &ref_mers_len.1;
                    for j in offset..offset+hit.hit_count {
                        let tf_idf = 1.0/(hit.hit_count as f64) * (idf[&ref_mers[j].0] as f64).log10();
                        mapping_tf_idf += tf_idf;
                    }
                }
                if mapping_tf_idf > max_fwd_tf_idf {max_fwd_tf_idf = mapping_tf_idf; max_fwd_mapping = i;}
               /* if mapping.len() > 1 {
                    println!("MAPPING {} TF-IDF {}", i+1, mapping_tf_idf);
                    for hit in mapping.iter() {
                        println!("{:?}", hit);
                    }
                }*/
            }
            // -----------------------CALCULATING RC TF-IDF-----------------------
            let mut max_rc_tf_idf = 0.0;
            let mut max_rc_mapping = 0;
            for i in 0..rc_mappings.len() {
                let mapping = &rc_mappings[i];
                let mut mapping_tf_idf = 0.0;
                for hit in mapping.iter() {
                    let offset = hit.offset;
                    let ref_mers_len = ref_index.get(&ref_id).unwrap();
                    let ref_mers = &ref_mers_len.1;
                    for j in offset..offset+hit.hit_count {
                        let tf_idf = 1.0/(hit.hit_count as f64) * (idf[&ref_mers[j].0] as f64).log10();
                        mapping_tf_idf += tf_idf;
                    }
                }
                if mapping_tf_idf > max_rc_tf_idf {max_rc_tf_idf = mapping_tf_idf; max_rc_mapping = i;}
               /* if mapping.len() > 1 {
                    println!("MAPPING {} TF-IDF {}", i+1, mapping_tf_idf);
                    for hit in mapping.iter() {
                        println!("{:?}", hit);
                    }
                }*/
            }
            let mut final_ref_s = 0;
            let mut final_ref_e = 0;
            let mut final_query_s = 0;
            let mut final_query_e = 0;
            if max_fwd_tf_idf > max_rc_tf_idf {
                let max_mapping = fwd_mappings[max_fwd_mapping].to_vec();
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
                let mut v = (query_id.to_string(), ref_id.to_string(), query_len, ref_len, final_query_s, final_query_e, final_ref_s, final_ref_e, max_fwd_tf_idf, false);
                println!("FWD:\n{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", query_id, final_query_s, final_query_e, ref_id, final_ref_s, final_ref_e, max_fwd_tf_idf, false);
                final_matches.push(v);
            }
            else if max_fwd_tf_idf < max_rc_tf_idf {
                let max_mapping = rc_mappings[max_rc_mapping].to_vec();
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
                let mut v = (query_id.to_string(), ref_id.to_string(), query_len, ref_len, final_query_s, final_query_e, final_ref_s, final_ref_e, max_rc_tf_idf, true);
                println!("RC:\n{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", query_id, final_query_s, final_query_e, ref_id, final_ref_s, final_ref_e, max_rc_tf_idf, true);
                final_matches.push(v);
            }
            /*
            //rc_hits.retain(|x| x.hit_count > k);
            /*if rc_hits.len() != 0 {
                //rc_hits.sort_by(|a, b| ((a.query_s)).cmp(&((b.query_s))));
                let mut mappings = Vec::<Vec::<&Hit>>::new();
                let mut mapping_ends = Vec::<&Hit>::new();
                let mut prev_hit = rc_hits[0];
                let mut prev_ref_s = prev_hit.ref_s;
                let mut prev_ref_e = prev_hit.ref_e;
                println!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", query_id, prev_hit.query_span, prev_hit.query_s, prev_hit.query_e, ref_id, prev_hit.ref_span, prev_hit.ref_s, prev_hit.ref_e, prev_hit.hit_count, prev_hit.is_rc);
                let mut current_mapping = Vec::<&Hit>::new();
                //current_mapping.push(prev_hit);
                for i in 1..rc_hits.len() {
                    let mut curr_hit = rc_hits[i];
                    let mut curr_ref_s = curr_hit.ref_s;
                    let mut curr_ref_e = curr_hit.ref_e;
                    let mut continue_sparse : bool = false;
                    if !((prev_ref_e + l - 1 >= curr_ref_s) && (prev_ref_e < curr_ref_s)) && !((curr_ref_e <= prev_ref_s + query_len) && (prev_ref_s < curr_ref_e)) {
                        for j in 0..mapping_ends.len() {
                            let mut mapping_end = mapping_ends[j];
                            let mut prev_end = mapping_end.ref_e;
                            let mut prev_start = mapping_end.ref_s;
                            //println!("Prev end {}\tprev start {}\tcurr end {}\tcurr start {}\tpredicate 1 {}\tpredicate 2 {}", prev_end, prev_start, curr_ref_e, curr_ref_s, ((prev_end + l - 1 >= curr_ref_s) && (prev_end < curr_ref_s)), ((curr_ref_e <= prev_start + query_len) && (prev_start < curr_ref_e)));  
                            if ((prev_end + l - 1 >= curr_ref_s) && (prev_end < curr_ref_s)) || ((curr_ref_e <= prev_start + query_len) && (prev_start < curr_ref_e)) {
                                println!("NEW MAPPING JOINED WITH {}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", prev_end, query_id, curr_hit.query_span, curr_hit.query_s, curr_hit.query_e, ref_id, curr_hit.ref_span, curr_hit.ref_s, curr_hit.ref_e, curr_hit.hit_count, curr_hit.is_rc);
                                current_mapping = mappings[j].to_vec();
                                current_mapping.push(curr_hit);
                                prev_ref_s = curr_ref_s;
                                prev_ref_e = curr_ref_e;
                                continue_sparse = true;
                                break;
                            }
                        }
                        if continue_sparse == false {
                            println!("NEW MAPPING STARTED {}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", query_id, curr_hit.query_span, curr_hit.query_s, curr_hit.query_e, ref_id, curr_hit.ref_span, curr_hit.ref_s, curr_hit.ref_e, curr_hit.hit_count, curr_hit.is_rc);
                            mappings.push(current_mapping);
                            mapping_ends.push(prev_hit);
                            current_mapping = Vec::<&Hit>::new();
                            current_mapping.push(curr_hit);
                            prev_ref_s = curr_ref_s;
                            prev_ref_e = curr_ref_e;
                        }
                        prev_hit = curr_hit;
                    }
                    else {
                        current_mapping.push(curr_hit);
                        prev_ref_s = curr_ref_s;
                        prev_ref_e = curr_ref_e;
                        println!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", query_id, curr_hit.query_span, curr_hit.query_s, curr_hit.query_e, ref_id, curr_hit.ref_span, curr_hit.ref_s, curr_hit.ref_e, curr_hit.hit_count, curr_hit.is_rc);
                        prev_hit = curr_hit;
                    }
                }
                mappings.push(current_mapping);
                let mut max_score = 0;
                for mapping in mappings.iter() {
                    let score = mapping.iter().map(|hit| (hit.hit_count - k) * hit.query_span).sum();
                    if score > max_score {max_score = score; rc_hits_cleared = mapping.to_vec();}
                }
            }
            
            //fwd_hits.retain(|x| x.hit_count > k);
            //println!("FW");
            if fwd_hits.len() != 0 {
                //fwd_hits.sort_by(|a, b| ((a.query_s)).cmp(&((b.query_s))));
                let mut mappings = Vec::<Vec::<&Hit>>::new();
                let mut mapping_ends = Vec::<&Hit>::new();
                let mut prev_hit = fwd_hits[0];
                let mut prev_ref_s = prev_hit.ref_s;
                let mut prev_ref_e = prev_hit.ref_e;
                println!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", query_id, prev_hit.query_span, prev_hit.query_s, prev_hit.query_e, ref_id, prev_hit.ref_span, prev_hit.ref_s, prev_hit.ref_e, prev_hit.hit_count, prev_hit.is_rc);
                let mut current_mapping = Vec::<&Hit>::new();
                //current_mapping.push(prev_hit);
                for i in 1..fwd_hits.len() {
                    let mut curr_hit = fwd_hits[i];
                    let mut curr_ref_s = curr_hit.ref_s;
                    let mut curr_ref_e = curr_hit.ref_e;
                    let mut continue_sparse : bool = false;
                    if !((prev_ref_e + l - 1 >= curr_ref_s) && (prev_ref_e < curr_ref_s)) && !((curr_ref_e <= prev_ref_s + query_len) && (prev_ref_s < curr_ref_e)) {
                        for j in 0..mapping_ends.len() {
                            let mut mapping_end = mapping_ends[j];
                            let mut prev_end = mapping_end.ref_e;
                            let mut prev_start = mapping_end.ref_s;
                            //println!("Prev end {}\tprev start {}\tcurr end {}\tcurr start {}\tpredicate 1 {}\tpredicate 2 {}", prev_end, prev_start, curr_ref_e, curr_ref_s, ((prev_end + l - 1 >= curr_ref_s) && (prev_end < curr_ref_s)), ((curr_ref_e <= prev_start + query_len) && (prev_start < curr_ref_e)));  
                            if ((prev_end + l - 1 >= curr_ref_s) && (prev_end < curr_ref_s)) || ((curr_ref_e <= prev_start + query_len) && (prev_start < curr_ref_e)) {
                                println!("NEW MAPPING JOINED WITH {}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", prev_end, query_id, curr_hit.query_span, curr_hit.query_s, curr_hit.query_e, ref_id, curr_hit.ref_span, curr_hit.ref_s, curr_hit.ref_e, curr_hit.hit_count, curr_hit.is_rc);
                                current_mapping = mappings[j].to_vec();
                                current_mapping.push(curr_hit);
                                prev_ref_s = curr_ref_s;
                                prev_ref_e = curr_ref_e;
                                continue_sparse = true;
                                break;
                            }
                        }
                        if continue_sparse == false {
                            println!("NEW MAPPING STARTED {}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", query_id, curr_hit.query_span, curr_hit.query_s, curr_hit.query_e, ref_id, curr_hit.ref_span, curr_hit.ref_s, curr_hit.ref_e, curr_hit.hit_count, curr_hit.is_rc);
                            mappings.push(current_mapping);
                            mapping_ends.push(prev_hit);
                            current_mapping = Vec::<&Hit>::new();
                            current_mapping.push(curr_hit);
                            prev_ref_s = curr_ref_s;
                            prev_ref_e = curr_ref_e;
                        }
                        prev_hit = curr_hit;
                    }
                    else {
                        current_mapping.push(curr_hit);
                        prev_ref_s = curr_ref_s;
                        prev_ref_e = curr_ref_e;
                        println!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", query_id, curr_hit.query_span, curr_hit.query_s, curr_hit.query_e, ref_id, curr_hit.ref_span, curr_hit.ref_s, curr_hit.ref_e, curr_hit.hit_count, curr_hit.is_rc);
                        prev_hit = curr_hit;
                    }
                }
                mappings.push(current_mapping);
                let mut max_score = 0;
                for mapping in mappings.iter() {
                    let score = mapping.iter().map(|hit| (hit.hit_count - k) * hit.query_span).sum();
                    if score > max_score {max_score = score; fwd_hits_cleared = mapping.to_vec();}
                }
            }*/
            /*if fwd_hits_cleared.len() != 0 {
                fwd_counts = fwd_hits_clearedlet mut v = (query_id.to_string(), ref_id.to_string(), query_len, ref_len, 0, query_len - 1, fwd_hits[0].ref_s, fwd_hits[0].ref_s  + query_len - 1, fwd_score, false);
                println!("FWD:\n{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", query_id, 0, query_len - 1, ref_id, fwd_hits[0].ref_s, fwd_hits[0].ref_s  + query_len - 1, fwd_score, false);
                final_matches.push(v);.iter().map(|hit| hit.hit_count).sum();
                fwd_hit_avg = fwd_counts as f64 / fwd_hits_cleared.len() as f64;
            }
            if rc_hits_cleared.len() != 0 {
                rc_counts = rc_hits_cleared.iter().map(|hit| hit.hit_count).sum();
                rc_hit_avg = rc_counts as f64 / rc_hits_cleared.len() as f64;
            }*/
            let mut fwd_score : usize = 0;
            let mut rc_score : usize = 0;
            if fwd_hits.len() != 0 {fwd_score = fwd_hits.iter().map(|hit| (hit.hit_count - k) * hit.query_span).sum();}
            if rc_hits.len() != 0 {rc_score = rc_hits.iter().map(|hit| (hit.hit_count - k) * hit.query_span).sum();}
            //println!("FWD mapping: {:?} RC mapping: {:?}", fwd_hits, rc_hits);
            println!("FWD score: {:?} RC score: {:?}", fwd_score, rc_score);

            if fwd_score > rc_score {write_fwd = true; write_rc = false;}
            else if rc_score > fwd_score {write_fwd = false; write_rc = true;}
            else {continue;}
            if fwd_mono == false {write_fwd = false;}
            if rc_mono == false {write_rc = false;}
            //println!("FWD:");

            if write_fwd == true {
                //println!("Selected FW");
                let mut v = (query_id.to_string(), ref_id.to_string(), query_len, ref_len, 0, query_len - 1, fwd_hits[0].ref_s, fwd_hits[0].ref_s  + query_len - 1, fwd_score, false);
                println!("FWD:\n{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", query_id, 0, query_len - 1, ref_id, fwd_hits[0].ref_s, fwd_hits[0].ref_s  + query_len - 1, fwd_score, false);
                final_matches.push(v);

            }
            else if write_rc == true {
                //println!("Selected RC\n______________________");
                let mut v = (query_id.to_string(), ref_id.to_string(), query_len, ref_len, 0, query_len - 1, rc_hits[0].ref_s, rc_hits[0].ref_s  + query_len - 1, rc_score, true);
                println!("RC:\n{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", query_id, 0, query_len - 1, ref_id, rc_hits[0].ref_s, rc_hits[0].ref_s  + query_len - 1, rc_score, true);
                final_matches.push(v);

            }
            else {continue;}
               println!("_____________________");
            */
        }     
    }
    return final_matches
}

pub fn output_paf(all_matches: &Vec<(Vec<Match>, String)>, paf_file: &mut File) {
    for (matches, id) in all_matches.iter() {
        let mut max_score = 0.0;
        let mut max_i = 0;
        for i in 0..matches.len() {
            let score = matches[i].8;
            //if score as usize <= 1 {continue;}
            //if score > max_score {max_score = score; max_i = i;}
           // if max_score == 0.0 {continue;}
            let v = &matches[i];    //let v = &matches[max_i];
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
}

