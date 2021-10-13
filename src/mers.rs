use std::collections::VecDeque;
use crate::DashMap;
use crate::File;
use std::io::Write;
use std::collections::HashMap;
use std::collections::HashSet;
use super::Params;
use dashmap::DashSet;
use crate::graph::DbgEntry;
use crate::graph::DbgIndex;
use crate::graph::Kmer;
use crate::kmer_vec::KmerVec;
use std::sync::Arc;

pub type Match = (String, String, usize, usize, usize, usize, usize, usize, usize, bool);
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
    count: usize,

}
#[derive(Clone, Debug)]
pub struct NAM {
    ref_id: String,
    query_s: usize,
    query_e: usize,
    query_last_hit_pos: usize,
    ref_s: usize,
    ref_e: usize,
    ref_last_hit_pos: usize,
    hit_count: usize,
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

fn try_mer(seq: &[u8], params: &Params, i: usize) -> Option<(u64, usize)> {
    let mut pos = i;
    let l = params.l;
    let hash_bound = ((params.d as f64) * 4_usize.pow(l as u32) as f64) as u64; 
    let lmask : u64 = ((1 as u64) << 2*l) - 1;
    let lshift : u64 = (l as u64 - 1) * 2;
    let mut xl : [u64; 2] = [0; 2];
    for i in pos..seq.len() {
        let c = seq_nt4_table[seq[i] as usize];
        if c < 4 {
            xl[0] = (xl[0] << 2 | c as u64) & lmask;
            xl[1] = xl[1] >> 2 | ((3 - c) as u64) << lshift;
            let yl : u64 = match xl[0] < xl[1] {
                true => xl[0],
                false => xl[1]
            };
            let hash_l = hash(yl, lmask);
            if hash_l <= hash_bound {
                return Some((hash_l, i))
            }
        } 
        else {xl = [0; 2];}
    }
    return None;
}

fn generate_kmer<'a>(seq: &'a [u8], params: &'a Params, s_pos: usize, offset: usize) -> Option<Kmer> {
    let mut k = params.k;
    let mut i = s_pos;
    let mut num = 0;
    let mut start = 0;
    let mut end = 0;
    let mut hashes = Vec::<u64>::new();
    let mut result = try_mer(seq, params, i);
    while num < k {
        if result.is_none() {return None;}
        let (hash, pos) = result.unwrap();
        if num == 0 {start = pos;}
        hashes.push(hash);
        num += 1;
        if num == k {end = pos;}
        i = pos + 1;
        result = try_mer(seq, params, i);
    }
    let kmer = Kmer::make_from(&hashes, start, end, offset);
    return Some(kmer)

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

pub fn extract_ref_kmers(seq: &[u8], params: &Params, rc: bool) -> Vec<Kmer> {
    let l = params.l;
    let hash_bound = ((params.d as f64) * 4_usize.pow(l as u32) as f64) as u64; 
    let lmask : u64 = ((1 as u64) << 2*l) - 1;
    let mut kmers = Vec::new();
    let seq_len = seq.len();
    let mut xl : [u64; 2] = [0; 2];
    let lshift : u64 = (l as u64 - 1) * 2;
    let mut offset = 0;
    let mut start = 0;
    let mut kmer_opt = generate_kmer(seq, params, start, offset);
    while kmer_opt.is_some() {
        let mut kmer = kmer_opt.unwrap();
        offset += 1;
        start = kmer.start + 1;
        kmers.push(kmer.clone());
        kmer_opt = generate_kmer(seq, params, start, offset); 
    } 
    if rc {
        let mut kmers_rev = kmers.iter().map(|k| k.reverse()).collect::<Vec<Kmer>>();
        kmers_rev.reverse();
        return kmers_rev
    }

    return kmers
} 

pub fn extract_mers(seq: &[u8], params: &Params) -> (Vec<u64>, Vec<usize>) {
    let l = params.l;
    let hash_bound = ((params.d as f64) * 4_usize.pow(l as u32) as f64) as u64; 
    let s = params.s;
    let wmin = params.wmin;
    let wmax = params.wmax;
    let smask : u64 = ((1 as u64) << 2*s) - 1;
    let lmask : u64 = ((1 as u64) << 2*l) - 1;
    let t = f64::ceil((params.l - params.s + 1) as f64 / 2.0) as usize;
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
            xl[0] = (xl[0] << 2 | c as u64) & lmask;
            xl[1] = xl[1] >> 2 | ((3 - c) as u64) << lshift;
            xs[0] = (xs[0] << 2 | c as u64) & smask;
            xs[1] = xs[1] >> 2 | ((3 - c) as u64) << sshift;
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
                        }
                    }
                    else {
                        let mut new_minimizer = false;
                        let tuple = update_window(&mut qs, &mut qs_pos, qs_min_val, qs_min_pos, hash_s, i - s + 1, new_minimizer);
                        qs_min_val = tuple.0; qs_min_pos = tuple.1; new_minimizer = tuple.2;
                        if qs_min_pos == qs_pos[t-1] as i32 {
                            let yl : u64 = match xl[0] < xl[1] {
                                true => xl[0],
                                false => xl[1]
                            };
                            let hash_l = hash(yl, lmask);
                            seq_hashes.push(hash_l);
                            pos_to_seq_coord.push(i - l + 1);
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

pub fn seq_to_kmers(seq: &[u8], id: &str, params: &Params, read: bool, query_mers_index: &DashMap<u64, usize>, ftup: (usize, usize)) -> (usize, Vec<Kmer>, Vec<Kmer>) {
    let mut mers = Vec::<Kmer>::new();
    let mut mers_rev = Vec::<Kmer>::new();
    let k = params.k;
    let l = params.l;
    let s = params.s;
    let wmin = params.wmin;
    let wmax = params.wmax;
    let read_length = seq.len();
    if read_length < wmax {return (seq.len(), mers, mers_rev);}
    let lmask : u64 = ((1 as u64) << 2*l) - 1;
    let smask : u64 = ((1 as u64) << 2*s) - 1;
    let q = 2_u64.pow(16) - 1;
    let t = 3;
    let (mut string_hashes_raw, mut pos_to_seq_coord_raw) = extract_mers(seq, params);
    let mut string_hashes = Vec::new();
    let mut pos_to_seq_coord = Vec::new();
    if !read {
        for i in 0..string_hashes_raw.len() {
            if query_mers_index.get(&string_hashes_raw[i]).is_some() {
                string_hashes.push(string_hashes_raw[i]);
                pos_to_seq_coord.push(pos_to_seq_coord_raw[i]);
            }
        }
    }
    else {
        for i in 0..string_hashes_raw.len() {
            if (*query_mers_index.get(&string_hashes_raw[i]).unwrap().value() >= ftup.0) && (*query_mers_index.get(&string_hashes_raw[i]).unwrap().value() <= ftup.1) {
                string_hashes.push(string_hashes_raw[i]);
                pos_to_seq_coord.push(pos_to_seq_coord_raw[i]);
            }
        }
    }
    let num_hashes = string_hashes.len();
    let mut string_hashes_rev = string_hashes.clone();
    string_hashes_rev.reverse();
    let mut pos_to_seq_coord_rev = pos_to_seq_coord.clone();
    pos_to_seq_coord_rev.reverse();
    for i in 0..num_hashes {
        let mut kmer_hashes : Vec<u64> = Vec::new();
        let mut kmer_hashes_rev : Vec<u64> = Vec::new();
        let mut kmer_start = 0;
        let mut kmer_end = 0;
        let mut kmer_start_rev = 0;
        let mut kmer_end_rev = 0;
        let mut lmer_counts = Vec::<usize>::new();
        if i + k - 1 < num_hashes {
            for j in i..i+k {
                /*if read {lmer_counts.push(*query_mers_index.get(&string_hashes[j]).unwrap().value());}*/
                if params.strobe {
                    let s = get_randstrobe(j, wmin, wmax, num_hashes, &string_hashes, &pos_to_seq_coord, q);
                    if s.is_none() {break;}
                    let s_unw = s.unwrap();
                    kmer_hashes.push(s_unw.0);
                    if j == i + k - 1 {kmer_end = s_unw.2;}
                }
                else {
                    kmer_hashes.push(string_hashes[j]);
                    kmer_hashes_rev.push(string_hashes_rev[j]);
                }
            }
            kmer_start_rev = pos_to_seq_coord_rev[i];
            kmer_end_rev = pos_to_seq_coord_rev[i + k - 1];
            kmer_start = pos_to_seq_coord[i];
            kmer_end = pos_to_seq_coord[i + k - 1];
        }
        if kmer_end != 0 {
            let kmer = KmerVec::make_from(&kmer_hashes, kmer_start, kmer_end, i);
            let kmer_rev = KmerVec::make_from(&kmer_hashes_rev,  kmer_end_rev, kmer_start_rev, i);
           // let (norm, rev) = kmer.normalize();
            mers.push(kmer.clone());

            if read {
                mers_rev.push(kmer_rev);
            }
        }   
    }
    return (seq.len(), mers, mers_rev);
}

pub fn edge_util(node: &DbgEntry, dbg_nodes: &DashMap<Vec<u64>, DbgEntry>, dbg_edges: &DashMap<DbgIndex, Vec<Kmer>>, mut count: usize, visited: &mut Vec<DbgIndex>, kmers: &Vec<Kmer>, kmers_rev: &Vec<Kmer>, rc: bool, mut depth: usize) -> (DbgEntry, usize){
    visited.push(node.index);
    //println!("{}", count);
    if count == kmers.len() - 1 {return (node.clone(), depth)}
    let mut origin_eq = &node.origin.clone();
    let edges = dbg_edges.get(&node.index);
    if edges.is_some() {
        for out in edges.unwrap().iter() {
            let out_node = dbg_nodes.get(&out.normalize().0.hashes);
            if count == kmers.len() - 1 {return (node.clone(), depth)}
            let mut kmer = kmers[count+1].clone();
            if rc {kmer = kmers_rev[count + 1].clone();}
            if &out_node.as_ref().unwrap().origin == origin_eq && out_node.as_ref().unwrap().abundance == 1 && !visited.contains(&out_node.as_ref().unwrap().index) && out_node.as_ref().unwrap().mers[0].normalize().0.hashes == kmer.normalize().0.hashes {
                count += 1;
                depth += 1;
                return edge_util(out_node.unwrap().value(), dbg_nodes, dbg_edges, count, visited, &kmers, &kmers_rev, rc, depth);
            }           
        }
    }
    return (node.clone(), depth);
}

pub fn new_query_graph(seq_id: &str, seq: &[u8],  params: &Params, dbg_nodes: &DashMap<Vec<u64>, DbgEntry>, dbg_edges: &DashMap<DbgIndex, Vec<Kmer>>, lens: &DashMap<String, usize>, rc: bool, query_mers_index: &DashMap<u64, usize>, ftup: (usize, usize)) -> Match {
    let mut count = 0;
    let mut pos = 0;
    let mut all_paths = Vec::<Vec<u64>>::new();
    let mut end_pos = 0;
    let mut out_rc = false;
    let mut q_start = 0; let mut q_end = 0;
    let mut r_start = 0; let mut r_end = 0;
    let mut ori = String::new();
    let mut first_done = false;
    let mut lengths = Vec::<usize>::new();
    let (seq_len, kmers, kmers_rev) = seq_to_kmers(seq, seq_id, params, true, query_mers_index, ftup);
    let mut visited = Vec::<DbgIndex>::new();
    while count < kmers.len() {
        let mut kmer = kmers[count].clone();
        if rc {kmer = kmers_rev[count].clone();}
        let mut node = dbg_nodes.get(&kmer.normalize().0.hashes);
        if node.is_some() && node.as_ref().unwrap().abundance == 1 {
            //println!("{}\t{:?}", seq_id, kmer);
            ori = node.as_ref().unwrap().origin[0].to_string();
            let (mut last_node, mut depth) = edge_util(&node.as_ref().unwrap()
            , dbg_nodes, dbg_edges, count, &mut visited, &kmers, &kmers_rev, rc, 1);
            lengths.push(depth);
            if !first_done {
                q_start = kmer.start;
                q_end = kmer.end;
                r_start = node.as_ref().unwrap().mers[0].start;
                r_end = last_node.mers[0].end;
                first_done = true;
            }
            else if ((last_node.mers[0].end as i32 - r_end as i32).abs() > 0 && (last_node.mers[0].end as i32 - r_end as i32).abs() < 100) {
                if rc {
                    r_end = last_node.mers[0].end; 
                    q_start = kmers_rev[count].start;
                }
                else {
                    r_end = last_node.mers[0].end;
                    q_end = kmers[count].end;
                }
            } 
            count += 1;
        }
        else {
            count += 1;
            continue;
        }
    }
    if ori.len() > 0 {
        if rc {
            r_end += (q_start - kmers_rev[kmers_rev.len() - 1].start);
            r_start -= (kmers_rev[0].end - q_end);
            q_start = kmers_rev[kmers_rev.len() - 1].start;
            q_end = kmers_rev[0].end;
        }
        else {
            r_start -= (q_start - kmers[0].start);
            r_end += (kmers[kmers.len() - 1].end - q_end);
            q_start = kmers[0].start;
            q_end = kmers[kmers.len() - 1].end; 
        }
        let mut score = lengths.iter().max().unwrap();
        let v = (seq_id.to_string(), ori.to_string(), seq.len(), *lens.get(&ori).unwrap().value(), q_start, q_end, r_start, r_end, *score, rc);
        println!("{:?}", v);
        return v;
    }
    else {
        let v = (seq_id.to_string(), ori.to_string(), seq.len(), 0, 0, q_end, 0, 0, 0, rc);
        return v;
    }
    
}

pub fn query_graph(query_id: &str, query_len: usize, query_mers: &Vec<Kmer>, dbg_nodes: &DashMap<Vec<u64>, DbgEntry>, dbg_edges: &DashMap<DbgIndex, Vec<Kmer>>, lens: &DashMap<String, usize>, rc: bool) -> Vec<Match> {
    let mut matches = Vec::<Match>::new();
    let mut visited : Arc::<DashMap::<Kmer, usize>> = Arc::new(DashMap::new());
    if query_mers.len() == 0 {return matches;}
    let mut start = dbg_nodes.get(&query_mers[0].hashes);
    let mut i = 1;
    while start.is_none() {
        if i == query_mers.len() {
            return vec![]
        }
        start = dbg_nodes.get(&query_mers[i].hashes);
        i += 1;
    }
    let start_n = start.unwrap();
    let (it, ref_start, ref_end, ref_ori) = dfs_util(i, &query_mers, &dbg_edges, &dbg_nodes, &visited, &mut 0, &mut 0, &mut Vec::new());
    if ref_ori.len() > 1 || ref_ori.len() == 0 {return Vec::new();}
    else {
        let v = (query_id.to_string(), ref_ori[0].to_string(), query_len, *lens.get(&ref_ori[0]).unwrap().value(), query_mers[i].start, query_mers[it-1].end, ref_start, ref_end, it-i, rc);
       // println!("{:?}", v);
        matches.push(v);
    }
    return matches
}

pub fn dfs_util(it: usize, query_mers: &Vec<Kmer>, dbg_edges: &DashMap<DbgIndex, Vec<Kmer>>, dbg_nodes: &DashMap<Vec<u64>, DbgEntry>, visited: &DashMap<Kmer, usize>, ref_start: &mut usize, ref_end: &mut usize, ref_ori: &mut Vec<String>) -> (usize, usize, usize, Vec<String>)  {
    if it == query_mers.len() {
        return (it, *ref_start, *ref_end, ref_ori.clone())
    }
    visited.insert(query_mers[it].clone(), it);
    let mut start = &query_mers[it];
    let n = &dbg_nodes.get(&start.hashes);
    if n.is_none() {
        return dfs_util(it+1, &query_mers, &dbg_edges, &dbg_nodes, &visited, ref_start, ref_end, ref_ori);
    }
    let idx = n.as_ref().unwrap().index;
    if dbg_edges.get(&idx).is_none() {
        return dfs_util(it+1, &query_mers, &dbg_edges, &dbg_nodes, &visited, ref_start, ref_end, ref_ori);   
    }
    else if it < query_mers.len() - 1 {
        for entry in dbg_edges.get(&idx).unwrap().iter() {
            let mut next = dbg_nodes.get(&entry.hashes).unwrap();
            if ref_ori.is_empty() {*ref_ori = next.origin.clone();}
            if !visited.contains_key(entry) && start.suffix() == entry.prefix() && next.origin == *ref_ori {
                if *ref_start == 0 {*ref_start = entry.start.clone();}
                *ref_end = entry.end.clone();
                //println!("{:?}\t{}\t{}\t{}\t{:?}\t{}\t{}\t{}", query_mers[it+1], query_mers[it+1].start, query_mers[it+1].end, query_mers[it+1].offset, entry.hashes, entry.start, entry.end, entry.offset);
                return dfs_util(it+1, &query_mers, &dbg_edges, &dbg_nodes, &visited, ref_start, ref_end, ref_ori);
            }
        }
    }
    return (it, *ref_start, *ref_end, ref_ori.clone()) 
}
/*pub fn find_hits(query_id: &str, query_len: usize, query_mers: &Vec<Mer>, ref_lens: &DashMap<String, usize>, query_mers_index: &DashMap<u64, usize>, mers_index: &DashMap<u64, Vec<(Mer, String)>>, params: &Params) -> Vec<Match> {
    let k = params.k;
    let l = params.l;
    let nam = params.nam;
    let z = params.z;
    let mut hits_per_ref = HashMap::<String, Vec<Hit>>::new();
    let mut hit_count_all = 0;
    let mut i = 0;
    let mut i_rev = 0;
    while i < query_mers.len() {
        let q = &query_mers[i];
        let mut h = Hit {query_id: query_id.to_string(), ref_id: String::new(), query_s: q.start, query_e: q.end, ref_s: 0, ref_e: 0, hit_count: 0, is_rc: q.rc, query_span: 0, ref_span: 0, query_offset: q.offset, ref_offset: 0, count: q.sum};
        let mut max_extend_offset = 0;
        let mer_entry = mers_index.get(&q.hash);
        if mer_entry.is_some() {       
            let vec = mer_entry.unwrap().clone();     
            let (mer, ref_id) = &vec[0];
            h.ref_offset = mer.offset;
            h.ref_id = ref_id.to_string();
            h.ref_s = mer.start;
            h.ref_e = mer.end;
            h.hit_count = 1;
            h.query_span = h.query_e - h.query_s + 1; 
            h.ref_span = h.ref_e - h.ref_s + 1;
            hits_per_ref.entry(ref_id.to_string()).or_insert(vec![]).push(h.clone());
            hit_count_all += 1; 
        }
        i += 1;
    }
    let mut final_matches = Vec::<Match>::new();
    if nam {
        for (key, val) in hits_per_ref.into_iter() {
            let mut open_nams = Vec::<NAM>::new();
            let ref_id = key;
            let ref_len = *ref_lens.get(&ref_id).unwrap().value();
            let mut hits_raw = val.clone();
            let mut prev_q_start = 0;
            let mut prev_rc = false;
            let mut hits_rc : Vec<&Hit> = hits_raw.iter().filter(|h| h.is_rc).collect();
            let mut hits_fwd : Vec<&Hit> = hits_raw.iter().filter(|h| !h.is_rc).collect();
            hits_rc.sort_by(|a, b| a.query_s.cmp(&b.query_s));
            hits_fwd.sort_by(|a, b| a.query_s.cmp(&b.query_s));
            let mut hits_part = vec![hits_rc, hits_fwd];
            for hits in hits_part.iter() {
                for i in 0..hits.len() {
                    let mut is_added = false;
                    let mut h = hits[i];
                    for o in open_nams.iter_mut() {
                        //println!("{:?}", o);
                        if (h.is_rc == o.is_rc) {
                            if (h.query_e > o.query_e) && (h.ref_e > o.ref_e) && h.query_e < (o.query_e + query_len) {
                                o.query_e = h.query_e;
                                o.ref_e = h.ref_e;
                                o.hit_count += 1;
                                is_added = true;
                                break;
                            }  
                        }
                    }
                    if !is_added {
                        let mut n = NAM {ref_id: ref_id.clone(), query_s: h.query_s, query_e: h.query_e, query_last_hit_pos: h.query_s, ref_s: h.ref_s, ref_e: h.ref_e, ref_last_hit_pos: h.ref_s, hit_count: 1, previous_query_start: h.query_s, previous_ref_start: h.ref_s, is_rc: h.is_rc};
        //                n.score += (float)1/ (float)h.hit_count;
                        open_nams.push(n);
                    }
                }
            }
            for n in open_nams.iter() {
                let v = (query_id.to_string(), n.ref_id.to_string(), query_len, ref_len, n.query_s, n.query_e + l, n.ref_s, n.ref_e + l, n.hit_count, n.is_rc);
                final_matches.push(v);
            }
        }
    }
    else {
        for (key, val) in hits_per_ref.into_iter() {
            let ref_id = key;
            let mut rc = false;
            let ref_len = *ref_lens.get(&ref_id).unwrap().value();
            let mut score = 0;
            let mut hits = val.to_vec();
            let mut sum : usize = hits.iter().map(|h| h.ref_offset).sum();
            let mut m = sum as f64 / hits.len() as f64;
            let mut sd = f64::sqrt(hits.iter().map(|c| f64::abs(c.ref_offset as f64 - m).powf(2.0)).sum::<f64>() / (hits.len() as f64));
            let mut sc_sum : usize = hits.iter().map(|h| h.count).sum();
            let mut sm = sc_sum as f64 / hits.len() as f64;
            let mut ssd = f64::sqrt(hits.iter().map(|c| f64::abs(c.count as f64 - sm).powf(2.0)).sum::<f64>() / (hits.len() as f64));
            hits.reverse();
            let mut lhits = hits.clone();
            for i in 0..lhits.len() {
                let hit = &lhits[i];
                let s = f64::abs((hit.count as f64 - sm) / ssd);
                let z = f64::abs((hit.ref_offset as f64 - m) / sd);
                if s >= params.z && z >= params.z {
                    //println!("S {}\tZ {}\tHit {:?}", s, z, hit);
                    hits.retain(|h| h.ref_offset != hit.ref_offset);
                }
            }
            let mut rc_count = hits.iter().filter(|h| h.is_rc).collect::<Vec<&Hit>>().len();
            let mut fwd_count = hits.iter().filter(|h| !h.is_rc).collect::<Vec<&Hit>>().len();
            if rc_count > fwd_count {hits.retain(|h| h.is_rc); rc = true;}
            else {hits.retain(|h| !h.is_rc); rc = false;}
            //for hit in hits.iter() {println!("{:?}", hit);}
            if hits.is_empty() {continue;}
            let mut final_ref_s = hits[0].ref_s;
            let mut final_ref_e = hits[hits.len()-1].ref_e;
            if final_ref_s > final_ref_e {let mut temp = final_ref_e; final_ref_e = final_ref_s; final_ref_s = temp;}
            let mut final_query_s = hits[0].query_s;
            let mut final_query_e = hits[hits.len()-1].query_e;
            let mut score = hits.len();
            let ref_id = &hits[0].ref_id;
            if !rc && hits.len() > 1 {
                let v = (query_id.to_string(), ref_id.to_string(), query_len, ref_len, final_query_s, final_query_e, final_ref_s, final_ref_e, score, false);
                final_matches.push(v);
            }
            else if rc && hits.len() > 1 {
                let v = (query_id.to_string(), ref_id.to_string(), query_len, ref_len, final_query_s, final_query_e, final_ref_s, final_ref_e, score, true);
                final_matches.push(v);
            }
        }
    }
    return final_matches

}
*/
pub fn output_paf(all_matches: &Vec<(Match, String)>, paf_file: &mut File) {
    for (v, id) in all_matches.iter() {
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

