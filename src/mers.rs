use std::collections::VecDeque;
use crate::DashMap;
use crate::File;
use std::io::Write;
use crate::kminmer::Kminmer;
use crate::graph::DAG;
use std::collections::HashMap;
use std::collections::HashSet;
use nthash::NtHashIterator;
use std::cmp;
use super::Params;
pub type Match = (String, String, usize, usize, usize, usize, usize, usize, usize, bool, usize);
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
    match_score: usize,
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

pub fn filter_hits(hits: &Vec<Hit>) -> (Vec<&Hit>, Vec<&Hit>) {
    let mut fwd = hits.iter().filter(|h| !h.is_rc).collect::<Vec<&Hit>>();
    let mut rc = hits.iter().filter(|h| h.is_rc).collect::<Vec<&Hit>>();
    (fwd, rc)
}

pub fn score(h: &Vec<&Hit>) -> usize {
    return h.iter().map(|h| h.match_score).sum::<usize>();
}
pub fn hits(h: &Vec<&Hit>) -> usize {
    return h.iter().map(|h| h.hit_count).sum::<usize>();
}

pub fn mapq_calc(h1: &Vec<&Hit>, h2: &Vec<&Hit>) -> usize {
    let f1 = score(&h1) as f64;
    let f2 = score(&h1) as f64;
    return (40.0 * (1.0 - (f2 / f1)) * 1.0_f64.min((hits(&h1) as f64 /10.0)) * f1.ln()) as usize;
}

pub fn get_edges(hits: &Vec<&Hit>, params: &Params) -> Vec<(usize, usize, usize)> {
    let mut info = Vec::<(usize, usize, usize)>::new();
    let k = params.k; let l = params.l;
    for u in 0..hits.len() - 1 {
        for v in u + 1..hits.len() {
            if hits[u].is_rc == hits[v].is_rc && !hits[u].is_rc {
                let offset_difference = ((hits[v].query_offset - hits[u].query_offset) as i32 - (hits[v].ref_offset - hits[u].ref_offset) as i32).abs();
                if hits[u].ref_s < hits[v].ref_s &&
                hits[u].ref_e < hits[v].ref_e &&
                hits[u].ref_offset < hits[v].ref_offset &&
                hits[u].query_s < hits[v].query_s &&
                hits[u].query_e < hits[v].query_e &&
                hits[u].query_offset < hits[v].query_offset &&
                (offset_difference as usize * l) <= hits[u].match_score {
                    info.push((u, v, (hits[u].match_score - (offset_difference as usize * l))));
                }
            }
            else if hits[u].is_rc == hits[v].is_rc && hits[u].is_rc {
                let offset_difference =  ((hits[v].query_offset - hits[u].query_offset) as i32 - (hits[u].ref_offset - hits[v].ref_offset) as i32).abs();
                if hits[u].ref_s > hits[v].ref_s &&
                hits[u].ref_e > hits[v].ref_e &&
                hits[u].ref_offset > hits[v].ref_offset &&
                hits[u].query_s < hits[v].query_s &&
                hits[u].query_e < hits[v].query_e &&
                hits[u].query_offset < hits[v].query_offset &&
                (offset_difference as usize * l) <= hits[u].match_score {
                    info.push((u, v, (hits[u].match_score - (offset_difference as usize * l))));
                }
            }
        }
    }
    info
}

pub fn partition(hits: &mut Vec<&Hit>, params: &Params) -> (usize, usize) {
    let mut mapq = 60;
    let mut score = 0;
    if hits.len() > 1 {
        hits.sort_by(|a, b| a.query_offset.cmp(&b.query_offset));
        let edges = get_edges(hits, params);
        let dag = DAG::new(edges, hits.len());
        let mut longest_per_node = dag.longest_paths();
        if longest_per_node.is_empty() {return (0, 0);}
        longest_per_node.sort_by(|a, b| a.2.cmp(&b.2));
        longest_per_node.reverse();
        let (best_u, best_path, best_score) = &longest_per_node[0];
        if best_path.len() == 0  {return (0, 0);}
        let mut best_hits = Vec::<&Hit>::new();
        if params.debug {
            println!("BEST PATH");
            println!("-----------");
        }
        let hit = hits[*best_u];
        best_hits.push(hit);
        if params.debug {
            println!("SOURCE NODE:\nQID {}\tQSTART {}\tQEND {}\tQOFF {}\tRID {}\tRSTART {}\tREND {}\tROFF {}\tRC {}\tCOUNT {}\tSCORE {}", hit.query_id, hit.query_s, hit.query_e, hit.query_offset, hit.ref_id, hit.ref_s, hit.ref_e, hit.ref_offset, hit.is_rc, hit.hit_count, hit.match_score);
            println!("-----------");
        }
        for p in best_path.iter() {
            let hitp = hits[*p];
            best_hits.push(hitp);
            if params.debug {
                println!("PATH NODE:\nQID {}\tQSTART {}\tQEND {}\tQOFF {}\tRID {}\tRSTART {}\tREND {}\tROFF {}\tRC {}\tCOUNT {}\tSCORE {}", hitp.query_id, hitp.query_s, hitp.query_e, hitp.query_offset, hitp.ref_id, hitp.ref_s, hitp.ref_e, hitp.ref_offset, hitp.is_rc, hitp.hit_count, hitp.match_score);
            }
        }
        if params.debug {println!("-----------");
        println!("PATH SCORE {}", best_score);}
        *hits = best_hits;
        score = *best_score;

        /*let mut prev_offset = hits[0].ref_offset;
        let mut partitions = HashMap::<usize, Vec<&Hit>>::new();
        let mut prev_i = 0;
        partitions.insert(0, vec![&hits[0]]);
        for i in 1..hits.len() {
            let mut curr_offset = hits[i].ref_offset;
            if ((curr_offset as i32 - prev_offset as i32)).abs() <= 10 as i32 {
                partitions.entry(prev_i).or_insert(vec![]).push(&hits[i]);
            }
            else {
                partitions.insert(i, vec![&hits[i]]);
                prev_i = i;
            }
            prev_offset = curr_offset;
        }
        let mut v: Vec<_> = partitions.into_iter().collect();
        v.sort_by(|a, b| score(&a.1).cmp(&score(&b.1)));
        v.reverse();
        //println!("{:?}", v[0].1);
        mapq = match v.len() > 1 {
            true => mapq_calc(&v[0].1, &v[1].1),
            false => 60,
        };
        //println!("MAX SCORE {}\tMAPQ {}", score(&v[0].1), mapq);
        *hits = v[0].1.clone();
        */
    }
    (mapq, score)
}
pub fn find_coords(hits: &Vec<&Hit>, rc: bool, ref_id: &str, ref_len: usize, query_id: &str, query_len: usize, mapq: usize) -> Match {
    let mut final_ref_s = hits[0].ref_s;
    let mut final_ref_e = hits[hits.len()-1].ref_e;
    let mut final_query_s = hits[0].query_s;
    let mut final_query_e = hits[hits.len()-1].query_e;
    let mut score = score(hits);
    if rc {
        final_query_s = hits[hits.len()-1].query_s;
        final_query_e = hits[0].query_e;
    }
    if !rc {
        if query_len > final_ref_e - final_ref_s {
            final_ref_e += query_len - final_query_e;
            final_query_e = query_len - 1;
        }
        if final_ref_s > final_query_s {
            final_ref_s -= final_query_s;
            final_query_s = 0;
        }
    }
    else {
        if query_len - final_query_e < ref_len - final_ref_e {
            final_ref_s -= query_len - final_query_e;
            final_query_e = query_len - 1;
        }
        if ref_len - final_ref_e > final_query_s {
            final_ref_e += final_query_s;
            final_query_s = 0;
        }
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
    (query_id.to_string(), ref_id.to_string(), query_len, ref_len, final_query_s, final_query_e, final_ref_s, final_ref_e, score, rc, mapq)
}
pub fn extend_hit(index: usize, h: &mut Hit, query_mers: &Vec<Kminmer>, mers_index: &DashMap<Vec<u64>, Vec<(String, Kminmer)>>, prev_mer: &Kminmer, params: &Params) -> usize {
    if index == query_mers.len() - 1 {return h.hit_count;}
    let q = &query_mers[index + 1];
    let mer_entry = mers_index.get(&q.mers());
    if mer_entry.is_some() {
        let tup = mer_entry.unwrap();
        let mer_hash = tup.key();
        let mer_id_vec = tup.value();
        for (ref_id, mer) in mer_id_vec.iter() {
            if ref_id == &h.ref_id && (mer.offset as i32 - prev_mer.offset as i32).abs() == 1 as i32 && ((q.rev != mer.rev) == h.is_rc) {
                //println!("Q\t{}\t{}\tR\t{}\t{}\t{}", q.start, q.end, mer.start, mer.end, (q.rev != mer.rev));
                if h.is_rc {
                    h.ref_s = mer.start;
                }
                else {
                    h.ref_e = mer.end;
                }
                h.query_e = q.end;
                h.query_span = q.end - h.query_s + 1;
                h.ref_span = mer.end - h.ref_s + 1;
                h.hit_count += 1;
                let new_score = match_score_extend(q, &query_mers[index], mer, prev_mer, h.match_score, params);
                h.match_score += new_score;
                //println!("QUERY ID: {}\tstart {}\tend {}\toffset {}\nREF ID: {}\tstart {}\tend {}\toffset {}\nSCORE: {}", h.query_id.to_string(), q.start, q.end, q.offset, ref_id.to_string(), mer.start, mer.end, mer.offset, h.match_score);
                extend_hit(index + 1, h, query_mers, mers_index, mer, params);
            }
        }
    }

    return h.hit_count;
}

pub fn match_score(q: &Kminmer, r: &Kminmer, params: &Params) -> usize {
    cmp::min(cmp::min((q.end - q.start), (r.end - r.start)), (params.k * params.l))
}
pub fn match_score_extend(curr_q: &Kminmer, prev_q: &Kminmer, curr_r: &Kminmer, prev_r: &Kminmer, prev_score: usize, params: &Params) -> usize {
    cmp::min(cmp::min(curr_q.end - prev_q.end, curr_r.end - prev_r.end), params.l)
}

pub fn chain_hits(query_id: &str, query_len: usize, query_mers: &Vec<Kminmer>, mers_index: &DashMap<Vec<u64>, Vec<(String, Kminmer)>>, params: &Params) -> HashMap<String, Vec<Hit>> {
    let mut hits_per_ref = HashMap::<String, Vec<Hit>>::new();
    let mut hit_count_all = 0;
    let l = params.l;
    let k = params.k;
    let g = params.g;
    let mut i = 0;
    while i < query_mers.len() {
        let q = &query_mers[i];
        let mer_entry = mers_index.get(&q.mers());
        let mut count = 1;
        if mer_entry.is_some() {
            let tup = mer_entry.unwrap();
            let mer_hash = tup.key();
            let mer_id_vec = tup.value();
            if mer_id_vec.len() > params.f && params.f != 0 {i += 1; continue;}
            //if mer_id_vec.len() > 1 
            for (ref_id, mer) in mer_id_vec.iter() {
                let mut match_score = match_score(q, mer, params);
                //println!("QUERY ID: {}\tstart {}\tend {}\toffset {}\nREF ID: {}\tstart {}\tend {}\toffset {}\nSCORE: {}", query_id.to_string(), q.start, q.end, q.offset, ref_id.to_string(), mer.start, mer.end, mer.offset, match_score);
                let mut h = Hit {query_id: query_id.to_string(), ref_id: ref_id.to_string(), query_s: q.start, query_e: q.end, ref_s: mer.start, ref_e: mer.end, hit_count: 1, match_score: match_score, is_rc: (q.rev != mer.rev), query_span: q.end - q.start + 1, ref_span: mer.start - mer.end + 1, query_offset: q.offset, ref_offset: mer.offset};
                let mut prev_offset = mer.offset;
                count = extend_hit(i, &mut h, query_mers, mers_index, mer, params);
                hits_per_ref.entry(ref_id.to_string()).or_insert(vec![]).push(h.clone());
                hit_count_all += 1; 
            }
            
        }
        i += count + 1;
    }
    hits_per_ref
}

pub fn find_hits(query_id: &str, query_len: usize, query_mers: &Vec<Kminmer>, ref_lens: &DashMap<String, usize>, mers_index: &DashMap<Vec<u64>, Vec<(String, Kminmer)>>, params: &Params) -> Vec<Match> {
    let g = params.g;
    let mut hits_per_ref = chain_hits(query_id, query_len, query_mers, mers_index, params);
    let mut final_matches = Vec::<Match>::new();    
    for (key, val) in hits_per_ref.into_iter() {
        let ref_id = key;
        let mut hits_raw = val.to_vec();
        let ref_len = *ref_lens.get(&ref_id).unwrap().value();
        if !hits_raw.is_empty() {
            let (mut fwd, mut rc) = filter_hits(&hits_raw);
            rc.reverse();
            let (fwd_mapq, fwd_score) = partition(&mut fwd, params);
            let (rc_mapq, rc_score) = partition(&mut rc, params);
            if params.debug{println!("FWD HIT LEN {} SCORE {}\tRC HIT LEN {}\tSCORE {}", fwd.len(), fwd_score, rc.len(), rc_score);}
            let mut hits = match rc_score > fwd_score {
                true => rc,
                false => fwd,
            };
            if hits.is_empty() {continue;}
            let mut mapq = match rc_score > fwd_score {
                true => rc_mapq,
                false => fwd_mapq,
            };
            let v = find_coords(&hits, (rc_score > fwd_score), &ref_id, ref_len, query_id, query_len, mapq);
            final_matches.push(v);
        }
    }
    
    return final_matches
}

pub fn output_paf(all_matches: &Vec<(Vec<Match>, String)>, paf_file: &mut File, unmap_file: &mut File, params: &Params) {
    for (matches, id) in all_matches.iter() {
        if matches.len() == 0 {
            write!(unmap_file, "{}\n", id).expect("Error writing line.");
            continue;
        }
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
        let mapq = v.10;
        if mapq == 0 {
            write!(unmap_file, "{}\n", id).expect("Error writing line.");
        }
        if query_id == ref_id && query_s == ref_s && query_e == ref_e {continue;}
        let paf_line = format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n", query_id, query_len, query_s, query_e, rc, ref_id, ref_len, ref_s, ref_e, score, ref_len, mapq);
        write!(paf_file, "{}", paf_line).expect("Error writing line.");
    }
}

