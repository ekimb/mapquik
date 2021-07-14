use std::collections::VecDeque;
use crate::DashMap;
use crate::File;
use dashmap::ReadOnlyView;
use std::io::Write;
use std::cmp::Ordering;
use super::{MersVector, MersVectorReduced, MersVectorRead, PosIndex, KmerLookup};
use super::Params;

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

pub fn extract_lmers(seq: &[u8], params: &Params) -> (Vec<u64>, Vec<usize>) {
    let l = params.l;
    let density = params.density;
    let lmask : u64 = ((1 as u64) << 2*l) - 1;
    let mut hash_count = 0;
    let mut seq_hashes = Vec::new();
    let mut pos_to_seq_coord = Vec::new();
    let mut seq_len = seq.len();
    let mut lp = 0;
    let mut xl : [u64; 2] = [0; 2];
    let mut lshift : u64 = (l as u64 - 1) * 2;
    let hash_bound = ((density as f64) * u64::max_value() as f64) as u64; 
    for i in 0..seq.len() {
        let mut c = seq_nt4_table[seq[i] as usize];
        if c < 4 {
            xl[0] = (xl[0] << 2 | c as u64) & lmask;                  // forward strand
            xl[1] = xl[1] >> 2 | ((3 - c) as u64) << lshift;  // reverse strand
            lp += 1;
            let mut yl : u64 = match xl[0] < xl[1] {
                true => xl[0],
                false => xl[1]
            };
            let mut hash_l = hash(yl, lmask);
            if hash_l <= hash_bound {
                seq_hashes.push(hash_l);
                pos_to_seq_coord.push(i - l + 1);
                hash_count += 1;
            }
        } else {
            lp = 0;  xl = [0; 2];
        }
    }
    return (seq_hashes, pos_to_seq_coord)
} 

pub fn get_kminmer(i: usize, k: usize, num_hashes: usize, ref_id: u64, string_hashes: &Vec<u64>, pos_to_seq_coord: &Vec<usize>, reverse: bool) -> Option<(u64, u64, usize, usize, bool)> {
    if i + k - 1 < num_hashes {
        let mut kminmer_pos_end = pos_to_seq_coord[i];
        let mut kminmer_hash : u64 = 0;
        let mut hash_frac = 2;
        for j in i..i+k {
            let mut lmer_hash = string_hashes[j];
            kminmer_pos_end = pos_to_seq_coord[j];
            kminmer_hash += lmer_hash / hash_frac;
            hash_frac += 1;
        }
        let mut kminmer_pos_start = pos_to_seq_coord[i];
        let mut s = (kminmer_hash, ref_id, kminmer_pos_start, kminmer_pos_end, reverse);
        return Some(s);
    }
    else {return None;}
}

pub fn seq_to_kminmers_read(seq: &[u8], id: u64, params: &Params) -> MersVectorRead {
    let l = params.l;
    let k = params.k;
    let mut kminmers = MersVectorRead::new();
    let read_length = seq.len();
    let lmask : u64 = ((1 as u64) << 2*l) - 1;
    let mut string_hashes = Vec::<u64>::new();
    let mut pos_to_seq_coord = Vec::<usize>::new();
    let (mut string_hashes, mut pos_to_seq_coord) = extract_lmers(seq, params);
    let num_hashes = string_hashes.len();
    for i in 0..num_hashes + 1 {
        let s = get_kminmer(i, k, num_hashes, id, &string_hashes, &pos_to_seq_coord, false);
        if s.is_none() {continue;}
        kminmers.push(s.unwrap());
    }
    string_hashes.reverse();
    pos_to_seq_coord.reverse();
    for i in 0..num_hashes {
        pos_to_seq_coord[i] = read_length - pos_to_seq_coord[i] - l;
    }
    for i in 0..num_hashes + 1 {
        let s = get_kminmer(i, k, num_hashes, id, &string_hashes, &pos_to_seq_coord, true);
        if s.is_none() {continue;}
        kminmers.push(s.unwrap());
    }
    return kminmers;
}

pub fn seq_to_kminmers_ref(seq: &[u8], id: u64, params: &Params) -> MersVectorRead {
    let l = params.l;
    let k = params.k;
    let mut kminmers = MersVectorRead::new();
    let read_length = seq.len();
    let lmask : u64 = ((1 as u64) << 2*l) - 1;
    let mut string_hashes = Vec::<u64>::new();
    let mut pos_to_seq_coord = Vec::<usize>::new();
    let (mut string_hashes, mut pos_to_seq_coord) = extract_lmers(seq, params);
    let num_hashes = string_hashes.len();
    for i in 0..num_hashes + 1 {
        let s = get_kminmer(i, k, num_hashes, id, &string_hashes, &pos_to_seq_coord, false);
        if s.is_none() {continue;}
        kminmers.push(s.unwrap());
    }
    return kminmers;
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
    if read {return seq_to_kminmers_read(ref_seq, ref_id, params);}
    else {return seq_to_kminmers_ref(ref_seq, ref_id, params);}
}

pub fn find_nams(query_mers: &MersVectorRead, ref_mers: &MersVectorReduced, mers_index: &KmerLookup, l: usize, read: &[u8], filter_cutoff: usize) -> Vec<NAM> {
    let mut hits_per_ref = DashMap::<u64, Vec<Hit>>::new();
    let mut hit_count_reduced = 0;
    let mut hit_count_all = 0;
    let read_length = read.len();
    for q in query_mers.iter() {
        let mut h = Hit {ref_id: 0, query_s: q.2, query_e: q.3 + l, ref_s: 0, ref_e: 0, hit_count: 0, is_rc: q.4};
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
                    let ref_e = r.2 + l;
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
            if h.query_s > prev_q_start + l {
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

pub fn output_paf(all_nams: &mut Vec<(Vec<NAM>, u64)>, l: usize, id_hashes: ReadOnlyView<u64, String>, id_lengths: ReadOnlyView<u64, usize>, paf_file: &mut File) {
    for (nams, id) in all_nams.iter_mut() {
        if nams.len() == 0 {continue;}
        nams.retain(|x| x.n_hits > 1);

        let max_nam = nams.iter().max_by(|a, b| (a.n_hits * (a.query_e - a.query_s)).cmp(&(b.n_hits * (b.query_e - b.query_s))));
        if max_nam.is_none() {continue;}
        let n = max_nam.unwrap();
        //let n = nams.iter().max_by(|a, b| (a.n_hits).cmp(&(b.n_hits))).unwrap();
        let mut o = String::new();
        if n.is_rc {o = "-".to_string();}
        else {o = "+".to_string();}
        let query_id = id_hashes.get(&id).unwrap();
        let ref_id = id_hashes.get(&n.ref_id).unwrap();
        let query_len = id_lengths.get(&id).unwrap();
        let ref_len = id_lengths.get(&n.ref_id).unwrap();
        let paf_line = format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n", query_id, query_len, n.query_s, n.query_last_hit_pos + l, o, ref_id, ref_len, n.ref_s, n.ref_last_hit_pos + l, n.n_hits, n.ref_last_hit_pos + l - n.ref_s, "255");
        write!(paf_file, "{}", paf_line).expect("Error writing line.");
    }
}

