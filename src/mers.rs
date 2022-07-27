use std::collections::VecDeque;
use crate::DashMap;
use crate::File;
use std::io::Write;
use crate::kminmer::Kminmer;
use std::collections::HashMap;
use std::collections::HashSet;
use nthash::NtHashIterator;
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

pub fn encode_rle(inp_seq: &str) -> (String, Vec<usize>) {
    let mut prev_char = '#';
    let mut hpc_seq = String::new();
    let mut pos_vec = Vec::<usize>::new();
    let mut prev_i = 0;
    for (i, c) in inp_seq.chars().enumerate() {
        if c == prev_char && "ACTGactgNn".contains(c) {continue;}
        if prev_char != '#' {
            hpc_seq.push(prev_char);
            pos_vec.push(prev_i);
            prev_i = i;
        }
        prev_char = c;
    }
    hpc_seq.push(prev_char);
    pos_vec.push(prev_i);
    (hpc_seq, pos_vec)
}

pub fn extract(inp_seq_raw: &[u8], params: &Params) -> (Vec<u64>, Vec<usize>) {
    let density = params.density;
    let l = params.l;
    let mut read_minimizers_pos = Vec::<usize>::new();
    let mut read_transformed = Vec::<u64>::new();
    let hash_bound = ((density as f64) * (u64::max_value() as f64)) as u64;
    let mut tup = (String::new(), Vec::<usize>::new());
    let inp_seq = String::from_utf8(inp_seq_raw.to_vec()).unwrap();
    let mut seq;
    if !params.use_hpc {
        tup = encode_rle(&inp_seq); //get HPC sequence and positions in the raw nonHPCd sequence
        seq = tup.0; //assign new HPCd sequence as input
    }
    else {
        seq = inp_seq;  //already HPCd before so get the raw sequence
    }
    if seq.len() < l {
        return (read_transformed, read_minimizers_pos)
    }
    let iter = NtHashIterator::new(seq.as_bytes(), l).unwrap().enumerate().filter(|(_i, x)| *x <= hash_bound);
    for (i, hash) in iter {
        if !params.use_hpc {read_minimizers_pos.push(tup.1[i]);} //if not HPCd need raw sequence positions
        else {read_minimizers_pos.push(i);} //already HPCd so positions are the same
        read_transformed.push(hash);
    }
    return (read_transformed, read_minimizers_pos)

}

pub fn kminmers(sk: &Vec<u64>, pos: &Vec<usize>, params: &Params) -> Vec<Kminmer> {
    let k = params.k;
    let mut kminmers = Vec::<Kminmer>::new();
    if sk.len() >= k {
        for i in 0..(sk.len() - k + 1) {
            let kminmer = Kminmer::new(&sk[i..i+k], pos[i], pos[i + k - 1], i);
            kminmers.push(kminmer);
        }
    }
    kminmers
}

pub fn filter_hits(hits: &Vec<Hit>) -> (Vec<&Hit>, bool) {
    let mut fwd = hits.iter().filter(|h| !h.is_rc).collect::<Vec<&Hit>>();
    let mut rc = hits.iter().filter(|h| h.is_rc).collect::<Vec<&Hit>>();
    let fin = match fwd.len() >= rc.len() {
        true => (fwd, false),
        false => (rc, true),
    };
    fin

}

pub fn partition(hits: &mut Vec<&Hit>, g: usize) {
    if hits.len() > 1 {
        let mut prev_offset = hits[0].ref_offset;
        let mut partitions = HashMap::<usize, Vec<&Hit>>::new();
        let mut prev_i = 0;
        partitions.insert(0, vec![&hits[0]]);
        for i in 1..hits.len() {
            let mut curr_offset = hits[i].ref_offset;
            if ((curr_offset as i32 - prev_offset as i32)).abs() <= g as i32 {
                partitions.entry(prev_i).or_insert(vec![]).push(&hits[i]);
            }
            else {
                partitions.insert(i, vec![&hits[i]]);
                prev_i = i;
            }
            prev_offset = curr_offset;
        }
        let mut hits_f = partitions.iter().max_by(|a, b| a.1.len().cmp(&b.1.len())).unwrap().1.clone();
        *hits = hits_f;
    }
}
pub fn find_coords(hits: &Vec<&Hit>, rc: bool, ref_id: &str, ref_len: usize, query_id: &str, query_len: usize) -> Match {
    let mut final_ref_s = hits[0].ref_s;
    let mut final_ref_e = hits[hits.len()-1].ref_e;
    let mut final_query_s = hits[0].query_s;
    let mut final_query_e = hits[hits.len()-1].query_e;
    let mut score = hits.len();
    if final_ref_s > final_query_s {
        final_ref_s -= final_query_s;
        final_query_s = 0;
    }
    /*if final_ref_e - final_query_s > query_len {
        let excess_add = query_len / 2;
        if final_ref_s > excess_add {final_ref_s -= excess_add;}
        else {final_ref_s = 0;}
        if ref_len - final_ref_e > excess_add {final_ref_e += excess_add;}
        else {final_ref_e = ref_len;}
        if final_query_s > excess_add {final_query_s -= excess_add;}
        else {final_query_s = 0;}
        if query_len - final_query_e > excess_add {final_query_e += excess_add;}
        else {final_query_e = query_len;}
    }*/
    (query_id.to_string(), ref_id.to_string(), query_len, ref_len, final_query_s, final_query_e, final_ref_s, final_ref_e, score, rc)
}

pub fn find_hits(query_id: &str, query_len: usize, query_mers: &Vec<Kminmer>, ref_lens: &DashMap<String, usize>, mers_index: &DashMap<Kminmer, String>, params: &Params) -> Vec<Match> {
    let l = params.l;
    let k = params.k;
    let g = params.g;
    let mut hits_per_ref = HashMap::<String, Vec<Hit>>::new();
    let mut hit_count_all = 0;
    let mut i = 0;
    let mut i_rev = 0;
    for i in 0..query_mers.len() {
        let q = &query_mers[i];
        let mut h = Hit {query_id: query_id.to_string(), ref_id: String::new(), query_s: q.start, query_e: q.end, ref_s: 0, ref_e: 0, hit_count: 0, is_rc: q.rev, query_span: q.end - q.start + 1, ref_span: 0, query_offset: q.offset, ref_offset: 0};
        let mut max_extend_offset = 0;
        let mer_entry = mers_index.get(&q);
        if mer_entry.is_some() {
            let tup = mer_entry.unwrap();
            let mer = tup.key();
            let ref_id = tup.value();
            h.ref_offset = mer.offset;
            h.ref_id = ref_id.to_string();
            h.ref_s = mer.start;
            h.ref_e = mer.end;
            h.is_rc = (q.rev != mer.rev);
            h.hit_count += 1;
            h.ref_span = h.ref_e - h.ref_s + 1; 
            hits_per_ref.entry(ref_id.to_string()).or_insert(vec![]).push(h.clone());
            hit_count_all += 1; 
        }
    }
    let mut final_matches = Vec::<Match>::new();    
    for (key, val) in hits_per_ref.into_iter() {
        let ref_id = key;
        let mut hits_raw = val.to_vec();
        let ref_len = *ref_lens.get(&ref_id).unwrap().value();
        if !hits_raw.is_empty() {
            let (mut hits, rc) = filter_hits(&hits_raw);
            if rc { hits.reverse() }
            partition(&mut hits, g);
            let v = find_coords(&hits, rc, &ref_id, ref_len, query_id, query_len);
            final_matches.push(v);
        }
    }
    return final_matches
}

pub fn output_paf(all_matches: &Vec<(Vec<Match>, String)>, paf_file: &mut File, params: &Params) {
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
        if query_id == ref_id && query_s == ref_s && query_e == ref_e {continue;}
        let paf_line = format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n", query_id, query_len, query_s, query_e, rc, ref_id, ref_len, ref_s, ref_e, score, ref_len, "255");
        write!(paf_file, "{}", paf_line).expect("Error writing line.");
    }
}

