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
    pub fn hits(&self) -> Vec<Hit> {
        self.hits.to_vec()
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

    // Check if two Hits in the Chain are "compatible": 

    // a tuple of Hits (u, v) are compatible if:

    // query start locations, 
    // query end locations,
    // query k-min-mer offsets,
    // reference start locations, 
    // reference end locations,
    // reference k-min-mer offsets,

    // are strictly increasing in the forward direction, or

    // query start locations, 
    // query end locations,
    // query k-min-mer offsets,

    // are strictly increasing, and

    // reference start locations, 
    // reference end locations,
    // reference k-min-mer offsets,

    // are strictly decreasing in the reverse direction.

    // plus, the difference between 

    // the distance between the query starting locations of u and v //// the distance between the reference starting locations of u and v 
    
    // needs to be small (less than twice the distance between the query starting locations of u and v).
    pub fn is_edge(&self, i: usize, j: usize, q_len: usize) -> bool {
        let u = &self.nth(i);
        let v = &self.nth(j);
        if !u.rc {
            if u.r_start < v.r_start &&
            u.r_end < v.r_end &&
            u.r_offset < v.r_offset &&
            u.q_start < v.q_start &&
            u.q_end < v.q_end &&
            u.q_offset < v.q_offset &&
            (v.r_start - u.r_start) < (v.q_start - u.q_start) * 2 {
                return true;
            }
        }
        else if u.rc {
            if u.r_start > v.r_start &&
            u.r_end > v.r_end &&
            u.r_offset > v.r_offset &&
            u.q_start < v.q_start &&
            u.q_end < v.q_end &&
            u.q_offset < v.q_offset &&
            (u.r_start - v.r_start) < (v.q_start - u.q_start) * 2 {
                return true;
            }
        }
        return false;
    }

    // Check if a given chain is "consistent": 
    
    // A chain c is "consistent" if all consecutive tuples (u, v) of Hits are compatible.

    // This function does a linear pass between all consecutive tuples (u, v), and if u and v are not compatible, removes v from the Chain. If there are no incompatible tuples, it returns true; otherwise it returns false.
    pub fn consistent(&mut self, q_len: usize) -> bool {
        let mut rem_self_v = Vec::<Hit>::new();
        let mut prev_i = 0;
        for i in 1..self.len() {
            if !self.is_edge(prev_i, i, q_len) {
                rem_self_v.push(self.nth(i).clone());
            }
            else {prev_i = i;}
        }
        let mut rem_self = Chain::new(&rem_self_v);
        if rem_self_v.is_empty() {return true;}
        self.retain_complement(&rem_self);  
        return false;
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
        let c = self.hits.iter().filter(|h_p| (h_p.offset_diff() as i32 - h.offset_diff() as i32).abs() < g as i32).cloned().collect::<Vec<Hit>>();
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
    pub fn get_mapq(&mut self, g: usize) -> usize {
        let mut mapq = 0;
        let mut best_eqc = self.find_best_eqc(g);
        let mut second_eqc = self.find_best_2eqc(&best_eqc, g);
        if !second_eqc.is_empty() {
            mapq = kminmer_mapq(best_eqc.len(), second_eqc.len());
        }
        else if best_eqc.len() > 3 {mapq = 60;}
        else {mapq = 1;}
        if mapq < 60 {
            if best_eqc.len() <= 3 {mapq = 1;}
        }
        self.replace(&best_eqc);
        mapq
    }

    // Calculates the MAPQ score for the Chain (by computing the two longest equivalence classes), sorts the Chain by the query k-min-mer offsets, and checks for consistency.

    // MAPQ score defaults to 60 if the resulting Chain (the longest equivalence class) has length > 3, and 1 otherwise.
    pub fn partition(&mut self, params: &Params, q_len: usize) -> (usize, usize) {
        let mut g = 1000; // This could be added to Params later.
        let mut mapq = 0;
        let mut score = self.get_score();
        if self.len() > 1 {
            mapq = self.get_mapq(g);
            self.sort_by_q_offset();
            if self.consistent(q_len) {score = self.get_score();} 
            else {
                if self.len() > 3 {mapq = 60;}
                else {mapq = 1;}
            }
        }
        (mapq, score)
    }
    
    // Extends the first and last Hit locations in the Chain to the length of the query.
    pub fn find_coords(&mut self, rc: bool, r_id: &str, r_len: usize, q_id: &str, q_len: usize, mapq: usize) -> Match {
        let first = self.first();
        let last = self.last();
        let mut q_start = first.q_start;
        let mut q_end = last.q_end;
        let mut r_start = first.r_start;
        let mut r_end = last.r_end;
        if rc {
            r_start = last.r_start;
            r_end = first.r_end;
        }
        let mut score = self.get_count();
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
        let mut fwd = self.get_fwd();
        let mut rc = self.get_rc();
        let (fwd_mapq, fwd_score) = fwd.partition(&params, q_len);
        let (rc_mapq, rc_score) = rc.partition(&params, q_len);
        let mut is_rc = (rc_score > fwd_score);
        if is_rc {self.replace(&rc);}
        else {self.replace(&fwd);}
        if self.is_empty() {return None;}
        let mut mapq = match is_rc {
            true => rc_mapq,
            false => fwd_mapq,
        };
        Some(self.find_coords(is_rc, r_id, r_len, q_id, q_len, mapq))
    }

    // Obtains query and reference intervals that are not covered by a Hit in the final Chain object (optional, for base-level alignment).
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