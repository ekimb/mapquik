// chain.rs
// Contains the "Chain" struct, which is a Vec of Hits, and various necessary operations to filter out bad Hits.

use crate::{Entry, Hit, Index, Kminmer, Match, Params};
use std::collections::HashMap;
use std::fmt;

pub struct Chain {
    hits: Vec<Hit>,
}
impl Chain {

    // New Chain from &Vec of Hits.
    pub fn new(hits: &Vec<Hit>) -> Self {
        Chain {hits: hits.to_vec()}
    }

    // An empty Chain object.
    pub fn empty() -> Self {
        Chain {hits: Vec::new()}
    }

    // Generate a new Chain from Hits that are only in the forward direction.
    pub fn get_fwd(&self) -> Chain {
        let fwd = self.hits.iter().filter(|h| !h.rc).cloned().collect::<Vec<Hit>>();
        Chain::new(&fwd)
    }

    // Generate a new Chain from Hits that are only in the reverse direction.
    pub fn get_rc(&self) -> Chain {
        let rc = self.hits.iter().filter(|h| h.rc).cloned().collect::<Vec<Hit>>();
        Chain::new(&rc)
    }

    // Get the sum of all scores from the Hits in the Chain.
    pub fn get_score(&self) -> usize {
        return self.hits.iter().map(|h| h.score).sum::<usize>();
    }

    // Get the Hit that has the maximum score in the Chain.
    pub fn get_max_score(&self) -> &Hit {
        return self.hits.iter().max_by(|a, b| a.score.cmp(&b.score)).unwrap();
    }

    pub fn filter_hit_span(&mut self) {
        self.hits.retain(|h| h.span_diff() < (h.r_span()));
    }

    // Get the total difference in span (difference between length of query covered by a Hit and length of reference covered by the Hit) of the Chain.
    pub fn span_diff(&self) -> usize {
        let mut s = 0;
        for i in 0..self.len() {
            let h = self.nth(i);
            s += (h.q_span() as i32 - h.r_span() as i32).abs() as usize;
        }
        s
    }

    // Get total number of k-min-mer matches in the Chain (a Hit can have multiple consecutive k-min-mer matches).
    pub fn get_count(&self) -> usize {
        return self.hits.iter().map(|h| h.count).sum::<usize>();
    }

    // Get the Hit that has the maximum number of consecutive k-min-mer matches.
    pub fn get_max_count(&self) -> &Hit {
        return self.hits.iter().max_by(|a, b| a.count.cmp(&b.count)).unwrap();
    }

    // Get the number of Hits in the Chain.
    pub fn len(&self) -> usize {
        return self.hits.len()
    } 

    // Check if the Chain is empty.
    pub fn is_empty(&self) -> bool {
        self.hits.is_empty()
    }

    // Check if the Chain contains the Hit.
    pub fn contains(&self, h: &Hit) -> bool {
        self.hits.contains(h)
    }

    // Return a raw Vec of Hits.
    pub fn hits(&self) -> &Vec<Hit> {
        &self.hits
    }

    // Get the first Hit in the Chain.
    pub fn first(&self) -> &Hit {
        &self.hits[0]
    }

    // Get the last Hit in the Chain.
    pub fn last(&self) -> &Hit {
        &self.hits[self.len() - 1]
    }

    // Get the nth Hit in the Chain.
    pub fn nth(&self, i: usize) -> &Hit {
        &self.hits[i]
    }

    // Remove the ith Hit in the Chain.
    pub fn remove(&mut self, i: usize) {
        self.hits.remove(i);
    }

    // Reverse the order of Hits in the Chain.
    pub fn reverse(&mut self) {
        self.hits.reverse();
    }

    pub fn check_hit_compatible(&self, h1: &Hit, h2: &Hit) -> bool {
        let (u, v) = match h1.q_offset < h2.q_offset {
            true => (h1, h2),
            false => (h2, h1),
        };
        let mut u_q_s = u.q_start;
        let mut u_q_e = u.q_end;
        let mut v_q_s = v.q_start;
        let mut v_q_e = v.q_end;
        let mut u_r_s = u.r_start;
        let mut u_r_e = u.r_end;
        let mut v_r_s = v.r_start;
        let mut v_r_e = v.r_end;
        if u.rc {
            if (self.rc_inconsistent(u_r_s, v_r_s, u_r_e, v_r_e) || self.rc_gap_too_long(u_q_s, u_r_s, u_q_e, u_r_e, v_q_s, v_r_s, v_q_e, v_r_e)) {return false;}
        }
        else {
            if (self.fwd_inconsistent(u_r_s, v_r_s, u_r_e, v_r_e) || self.fwd_gap_too_long(u_q_s, u_r_s, u_q_e, u_r_e, v_q_s, v_r_s, v_q_e, v_r_e)) {return false;} 
        }
        return true;
    }

    pub fn retain_unique(&mut self) {
        let len = self.len();
        let mut remove_hit_indices = vec![false; len];
        for i in 0..self.len() {
            let h_i = self.nth(i);
            let h_c = self.hits.iter().filter(|hp| hp.r_start == h_i.r_start && hp.r_end == h_i.r_end).count();
            if h_c > 1 {remove_hit_indices[i] = true;}
        }
        self.hits = (0..len).filter(|i| !remove_hit_indices[*i]).map(|i| self.nth(i).clone()).collect::<Vec<Hit>>();

    }

    pub fn filter_hits(&mut self) {
        let mut hits_to_remove_per_hit = Vec::<(Vec<bool>, usize, usize)>::new();
        let len = self.len();
        for i in 0..len {
            let mut h_i = self.nth(i);
            let mut hits_to_remove = vec![true; len];
            hits_to_remove[i] = false;
            let mut hit_remove_count = len - 1;
            let mut hits_count = h_i.count;
            for j in 0..len {
                if i == j {continue;}
                let mut h_j = self.nth(j);
                let b = self.check_hit_compatible(h_i, h_j);
                if b {
                    hits_to_remove[j] = false;
                    hit_remove_count -= 1;
                    hits_count += h_j.count;
                }
            }
            hits_to_remove_per_hit.push((hits_to_remove, hit_remove_count, hits_count));
        }
        if hits_to_remove_per_hit.len() == 0 {self.hits = Vec::new(); return;}
        let mut idx_to_remove = &hits_to_remove_per_hit.iter().max_by(|a, b| a.2.cmp(&b.2)).unwrap().0;
        self.hits = (0..self.len()).filter(|i| !idx_to_remove[*i]).map(|i| self.nth(i).clone()).collect::<Vec<Hit>>();
        self.sort_by_q_offset();
    }


    pub fn fwd_gap_too_long(&self, u_q_s: usize, u_r_s: usize, u_q_e: usize, u_r_e: usize, v_q_s: usize, v_r_s: usize, v_q_e: usize, v_r_e: usize) -> bool {
       // (self.fwd_gap_start_too_long(u_q_s, u_r_s, v_q_s, v_r_s) ||self.fwd_gap_end_too_long(u_q_e, u_r_e, v_q_e, v_r_e) )

       let g_1 = v_q_s - u_q_e;
       let g_2 = v_r_s - u_r_e;
       (g_1 as i32 - g_2 as i32).abs() > 2000
    }

    pub fn rc_gap_too_long(&self, u_q_s: usize, u_r_s: usize, u_q_e: usize, u_r_e: usize, v_q_s: usize, v_r_s: usize, v_q_e: usize, v_r_e: usize) -> bool {
        //(self.rc_gap_start_too_long(u_q_s, u_r_s, v_q_s, v_r_s) ||self.rc_gap_end_too_long(u_q_e, u_r_e, v_q_e, v_r_e))
        let g_1 = v_q_s - u_q_e;
        let g_2 = u_r_s - v_r_e;
        (g_1 as i32 - g_2 as i32).abs() > 2000
    }


    pub fn fwd_start_inconsistent(&self, u_r_s: usize, v_r_s: usize) -> bool {
        (v_r_s <= u_r_s)
    }

    pub fn rc_start_inconsistent(&self, u_r_s: usize, v_r_s: usize) -> bool {
        (u_r_s <= v_r_s)
    }

    pub fn fwd_end_inconsistent(&self, u_r_e: usize, v_r_e: usize) -> bool {
        (v_r_e <= u_r_e)
    }

    pub fn rc_end_inconsistent(&self, u_r_e: usize, v_r_e: usize) -> bool {
        (u_r_e <= v_r_e)
    }

    pub fn fwd_inconsistent(&self, u_r_s: usize, v_r_s: usize, u_r_e: usize, v_r_e: usize) -> bool {
        (self.fwd_start_inconsistent(u_r_s, v_r_s) || self.fwd_end_inconsistent(u_r_e, v_r_e))
    }
    pub fn rc_inconsistent(&self, u_r_s: usize, v_r_s: usize, u_r_e: usize, v_r_e: usize) -> bool {
        (self.rc_start_inconsistent(u_r_s, v_r_s) || self.rc_end_inconsistent(u_r_e, v_r_e))
    }

    // Sort the Chain by the query k-min-mer offsets of the Hits.
    pub fn sort_by_q_offset(&mut self) {
        self.hits.sort_by(|a, b| a.q_offset.cmp(&b.q_offset));
    }

    // Sort the Chain by the reference k-min-mer offsets of the Hits.
    pub fn sort_by_r_offset(&mut self) {
        self.hits.sort_by(|a, b| a.q_offset.cmp(&b.r_offset));
    }


    // Find an "equivalence class" of a Hit in the Chain:

    // A Hit h_p is in the "equivalence class" E(h) of a Hit h if the absolute difference of the k-min-mer offset differences of h and h_p are below some threshold g (see hit.rs for a definition of offset difference in the context of a Hit).
    pub fn find_eq_class(&self, h: &Hit, g: usize) -> Chain {
        let c = self.hits.iter().filter(|h_p| (h_p.offset_diff() as i32 - h.offset_diff() as i32).abs() < g as i32 && self.check_hit_compatible(h, h_p)).cloned().collect::<Vec<Hit>>();
        Chain::new(&c)
    }

    // Generate a new Chain that has no elements that also occur in another Chain c.
    pub fn find_complement(&self, c: &Chain) -> Chain {
        let mut c_p = self.hits.iter().filter(|h| !c.contains(h)).cloned().collect::<Vec<Hit>>();
        Chain::new(&c_p)
    }

    // Modify this particular Chain so that it has no elements that also occur in another Chain c.
    pub fn retain_complement(&mut self, c: &Chain) {
        self.hits.retain(|h| !c.contains(h));
    }

    // Replace this Chain object with another Chain c.
    pub fn replace(&mut self, c: &Chain) {
        self.hits = c.hits().to_vec();
    }

    // Find the longest equivalence class: argmax(len(E(h))) for all Hits h in the Chain.
    pub fn find_best_eqc(&self, g: usize) -> Chain {
        let mut best_eqc = Chain::empty();
        for i in 0..self.len() {
            let mut eqc = self.find_eq_class(self.nth(i), g);
            if eqc.len() > best_eqc.len() {best_eqc = eqc;}
            else if eqc.len() == best_eqc.len() {
                if eqc.get_count() >= best_eqc.get_count() {
                    best_eqc = eqc;
                }
            }
        }
        best_eqc
    }

    // Find the second longest equivalence class in the Chain (this is useful for MAPQ calculations).
    pub fn find_best_2eqc(&self, best_eqc: &Chain, g: usize) -> Chain {
        let mut rest = self.find_complement(best_eqc);
        rest.find_best_eqc(g)
    }

    // Obtain a MAPQ score, given the two longest equivalence classes. 
    
    // MAPQ score defaults to 60 if the second longest equivalence class is empty, or if the longest equivalence class has length > 3.

    // For all pairs of equivalence classes with MAPQ < 60, MAPQ score defaults to 1 if the longest equivalence class has length <= 3.
    pub fn get_mapq(&mut self, g: usize, k: usize) -> usize {
        let mut mapq = 1;
        let mut best_eqc = self.find_best_eqc(g);
        let mut second_eqc = self.find_best_2eqc(&best_eqc, g);
        //println!("EQC1!CNT!{}!{}", best_eqc, best_eqc.get_count());
        //println!("EQC2!CNT!{}!{}", second_eqc, second_eqc.get_count());
        if !second_eqc.is_empty() {
            mapq = kminmer_mapq(best_eqc.get_count(), second_eqc.get_count());
        }
        else if best_eqc.len() > 3 || best_eqc.get_count() > (k + 1) * 3 {mapq = 60;}
        best_eqc.sort_by_q_offset();
        self.replace(&best_eqc);
        mapq
    }

    // Calculates the MAPQ score for the Chain (by computing the two longest equivalence classes), sorts the Chain by the query k-min-mer offsets, and checks for consistency.

    // MAPQ score defaults to 60 if the resulting Chain (the longest equivalence class) has length > 3, and 1 otherwise.
    pub fn partition(&mut self, params: &Params, q_len: usize) -> (usize, usize) {
        let mut k = params.k; // This could be added to Params later.
        let mut mapq = 0;
        let mut count = 0;
        let mut len = self.len();
        if len >= 1 {
            if len > 1 {
                self.retain_unique();
                if self.len() == 0 {return (mapq, count);}
                self.sort_by_q_offset();
                self.filter_hits();
            }
            count = self.get_count(); 
            if self.len() > 3 || count > k + 4 {mapq = 60;}
        }
        (mapq, count)
    }
    
    // Extends the first and last Hit locations in the Chain to the length of the query.
    pub fn find_coords(&mut self, rc: bool, r_id: &str, r_len: usize, q_id: &str, q_len: usize, mapq: usize, first: &Hit, last: &Hit, score: usize) -> Match {
        let mut q_start = first.q_start;
        let mut q_end = last.q_end;
        let mut r_start = first.r_start;
        let mut r_end = last.r_end;
        if rc {
            r_start = last.r_end;
            r_end = first.r_start;
        }
        let mut final_r_start = 0;
        let mut final_r_end = 0;
        let mut final_q_start = 0;
        let mut final_q_end = 0;
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
        let m = (q_id.to_string(), r_id.to_string(), q_len, r_len, final_q_start, final_q_end, final_r_start, final_r_end, score, rc, mapq);
        m
    }

    // Wrapper function that filters bad Hits, checks for consistency, and obtains final coordinates.

    // Outputs a Match object (see mers.rs for a definition).
    pub fn get_match(&mut self, r_id: &str, r_len: usize, q_id: &str, q_len: usize, params: &Params) -> Option<Match> {
        if self.is_empty() {return None;}
        let mut fwd = self.get_fwd();
        let mut rc = self.get_rc();
        let (fwd_mapq, fwd_score) = fwd.partition(&params, q_len);
        let (rc_mapq, rc_score) = rc.partition(&params, q_len);
        if fwd_score == rc_score && fwd_score == 0 {return None;}
        let mut is_rc = (rc_score > fwd_score);
        let (mapq, score, first, last) = match is_rc {
            true => (rc_mapq, rc_score, rc.first(), rc.last()),
            false => (fwd_mapq, fwd_score, fwd.first(), fwd.last()),
        };

        Some(self.find_coords(is_rc, r_id, r_len, q_id, q_len, mapq, first, last, score))
    }

    // Obtains query and reference intervals that are not covered by a Hit in the final Chain object (only for base-level alignment).
    pub fn get_remaining_seqs(&self, m: &Match) -> (Vec<(usize, usize)>, Vec<(usize, usize)>) {
        let mut q_coords = Vec::<(usize, usize)>::new();
        let mut r_coords = Vec::<(usize, usize)>::new();

        let (q_id, r_id, q_len, r_len, q_start, q_end, r_start, r_end, score, rc, mapq) = m;
        if !rc {
            let mut prev_q_end = q_start;
            let mut prev_r_end = r_start;
            let rc_s = match rc {
                true => "-",
                false => "+",
            };
            for i in 0..self.len() {
                let hit = self.nth(i);
                let rc_s = match hit.rc {
                    true => "-",
                    false => "+",
                };
                //println!("QALN!{}!{}!RALN!{}!{}", prev_q_end, hit.q_start, prev_r_end, hit.r_start);
                //println!("QMMM!{}!{}!RMMM!{}!{}", hit.q_start, hit.q_end, hit.r_start, hit.r_end); 
                if prev_q_end < &hit.q_start && &hit.r_start > prev_r_end {
                    q_coords.push((*prev_q_end, hit.q_start));
                    r_coords.push((*prev_r_end, hit.r_start));
                    prev_q_end = &hit.q_end;
                    prev_r_end = &hit.r_end;
                }
            }
            //println!("QALN!{}!{}!RALN!{}!{}", prev_q_end, q_end, prev_r_end, r_end);
            if prev_q_end != q_end {
                q_coords.push((*prev_q_end, *q_end));
                r_coords.push((*prev_r_end, *r_end));
            }

        }
        if *rc {
            let mut prev_q_end = q_start;
            let mut prev_r_start = r_end;
            let rc_s = match rc {
                true => "-",
                false => "+",
            };
            for i in 0..self.len() {
                let hit = self.nth(i);
                let rc_s = match hit.rc {
                    true => "-",
                    false => "+",
                };
                if hit.r_end < *r_start {break;}
                //println!("QALN!{}!{}!RALN!{}!{}", prev_q_end, hit.q_start, prev_r_end, hit.r_start);
                //println!("QMMM!{}!{}!RMMM!{}!{}", hit.q_start, hit.q_end, hit.r_start, hit.r_end); 
                if prev_q_end < &hit.q_start && &hit.r_end < prev_r_start {
                    q_coords.push((*prev_q_end, hit.q_start));
                    r_coords.push((hit.r_end, *prev_r_start));
                    prev_q_end = &hit.q_end;
                    prev_r_start = &hit.r_start;
                }
            }
            //println!("QALN!{}!{}!RALN!{}!{}", prev_q_end, q_end, r_start, prev_r_start);
            if prev_q_end != q_end {
                q_coords.push((*prev_q_end, *q_end));
                r_coords.push((*r_start, *prev_r_start));
            }
        }
        (q_coords, r_coords)
    } 
}

// Pretty-prints a Chain.
impl fmt::Display for Chain {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut res = String::new();
        res = format!("{}\n", res);
        for i in 0..self.len() {
            let h = self.nth(i);
            res = format!("{}{}", res, h);
        }
        write!(f, "{}", res)
    }
}

// Calculates a "MAPQ" score, given two lengths of equivalence classes.
pub fn kminmer_mapq(best_len: usize, second_len: usize) -> usize {
    (60.0 * (1.0 - (second_len as f64/best_len as f64))) as usize
}
