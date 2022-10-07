// chain.rs
// Contains the "Chain" struct, which is a Vec of Hits, and various necessary operations to filter out bad Hits.

use crate::{Entry, Hit, Index, Kminmer, Match, Params};
use std::collections::HashMap;
use std::fmt;

#[derive(Clone, Debug, PartialEq)]
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

    pub fn check_hit_compatible(&self, h1: &Hit, h2: &Hit, g: usize) -> bool {
        if h1 == h2 {return true;}
        if h1.rc != h2.rc {return false;}
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
            if (self.rc_inconsistent(u_r_s, v_r_s, u_r_e, v_r_e) || self.rc_gap_too_long(u_r_s, u_q_e, v_q_s, v_r_e, g)) {return false;}
        }
        else {
            if (self.fwd_inconsistent(u_r_s, v_r_s, u_r_e, v_r_e) || self.fwd_gap_too_long(u_q_e, u_r_e, v_q_s, v_r_s, g)) {return false;} 
        }
        return true;
    }

    pub fn retain_unique(&mut self) {
        let len = self.len();
        let mut remove_hit_indices = vec![false; len];
        for i in 0..self.len() {
            let h_i = self.nth(i);
            let h_c = self.hits.iter().filter(|hp| hp.r_start == h_i.r_start && hp.r_end == h_i.r_end && hp.rc == h_i.rc).count();
            if h_c > 1 {remove_hit_indices[i] = true;}
        }
        self.hits = (0..len).filter(|i| !remove_hit_indices[*i]).map(|i| self.nth(i).clone()).collect::<Vec<Hit>>();

    }

    pub fn colinear_idx_per_hit(&self, h: &Hit, g: usize) ->  Vec<usize> {
        let col = (0..self.len()).filter(|i| self.check_hit_compatible(h, self.nth(*i), g)).collect::<Vec<usize>>();
        col
    }
    pub fn colinear_hits_per_hit(&self, h: &Hit, g: usize) ->  Vec<&Hit> {
        let col = self.hits.iter().filter(|hp| self.check_hit_compatible(h, hp, g)).collect::<Vec<&Hit>>();
        col
    }

    pub fn check_colinear_hit_set(&self, col: &Vec<usize>, g: usize) -> (Vec<usize>, usize) {
        let mut remove_h = HashMap::<usize, Vec<usize>>::new();
        let mut remove_count : i32 = 0;
        let mut to_remove = vec![false; self.len()];

        for idx_i in 0..col.len() - 1 {
            let i = col[idx_i];
            let mut h_i = self.nth(i);
            for idx_j in idx_i..col.len() {
                let j = col[idx_j];
                let mut h_j = self.nth(j);
                if !self.check_hit_compatible(h_i, h_j, g) {
                    remove_h.entry(i).or_insert(Vec::new()).push(j);
                    remove_h.entry(j).or_insert(Vec::new()).push(i);
                    remove_count += 1;
                }
            }
        } 
        if remove_count > 0 {
            let mut remove_v = remove_h.into_iter().collect::<Vec<(usize, Vec<usize>)>>();
            remove_v.sort_by(|a, b| a.1.len().cmp(&b.1.len()));
            remove_v.reverse();
            for i in 0..remove_v.len() {
                let (max_to_remove, max_v) = &remove_v[i];
                //eprintln!("SAVED!{}!{}", remove_count, self.nth(*max_to_remove));
                to_remove[*max_to_remove] = true;
                remove_count -= max_v.len() as i32;
                if remove_count <= 0 {break;}
            }
        }
        let mut score = col.iter().map(|i| self.nth(*i).count).sum::<usize>();
        let col_f = (0..col.len()).map(|idx_i| col[idx_i]).filter(|i| !to_remove[*i]).collect::<Vec<usize>>();
        (col_f, score)
    }

    pub fn check_colinear(&mut self, g: usize) {
        let len = self.len();
        let mut max_col = Vec::<usize>::new();
        let mut max_score = 0;
        for i in 0..len {
            let h_i = self.nth(i);
            let col = self.colinear_idx_per_hit(h_i, g);
            let (col_f, score) = self.check_colinear_hit_set(&col, g);
            if score > max_score {
                max_col = col_f;
                max_score = score;
            }
        }
        self.hits = max_col.iter().map(|i| self.nth(*i).clone()).collect::<Vec<Hit>>();
        self.sort_by_q_offset();
    }

    pub fn filter_hits(&mut self, g: usize) {

        let len = self.len();
        if len == 0 {return;}
        let max_chain = self.colinear_hits_per_hit(self.hits.iter().max_by(|a, b| self.colinear_hits_per_hit(a, g).iter().map(|h| h.count).sum::<usize>().cmp(&self.colinear_hits_per_hit(b, g).iter().map(|hp| hp.count).sum::<usize>())).unwrap(), g);
        self.hits = max_chain.into_iter().cloned().collect();
        self.sort_by_q_offset();
    }

    pub fn find_largest_two_hits(&self) -> (usize, usize) {
        let mut max = 0;
        let mut max_count = 0;
        let mut second_max = 0;
        let mut second_max_count = 0;
        for i in 0..self.len() {
            let count = self.nth(i).count;
            if count > max_count {
                second_max = max;
                second_max_count = max_count;
                max = i;
                max_count = count;
            } else if count > second_max_count {
                second_max = i;
                second_max_count = count;
            }
        }
        (max, second_max)
    }

    pub fn filter_hits_max_two(&mut self, g: usize) {
        let len = self.len();
        let (max, second_max) = self.find_largest_two_hits();
        let max_chain = self.colinear_hits_per_hit(self.nth(max), g);
        let max_count = max_chain.iter().map(|h| h.count).sum::<usize>();
        self.hits = max_chain.into_iter().cloned().collect();
        let second_max_chain = self.colinear_hits_per_hit(self.nth(second_max), g);
        let second_max_count = second_max_chain.iter().map(|h| h.count).sum::<usize>();
        if second_max_count > max_count {
            self.hits = second_max_chain.into_iter().cloned().collect();
        }
        self.sort_by_q_offset();
    }

    pub fn filter_hits_max(&mut self, g: usize) {
        let len = self.len();
        let (max, second_max) = self.find_largest_two_hits();
        let max_chain = self.colinear_hits_per_hit(self.nth(max), g);
        self.hits = max_chain.into_iter().cloned().collect();
        self.sort_by_q_offset();
    }


    pub fn fwd_gap_too_long(&self, u_q_e: usize, u_r_e: usize, v_q_s: usize, v_r_s: usize, g: usize) -> bool {
       // (self.fwd_gap_start_too_long(u_q_s, u_r_s, v_q_s, v_r_s) ||self.fwd_gap_end_too_long(u_q_e, u_r_e, v_q_e, v_r_e) )

       let g_1 = v_q_s as i32 - u_q_e as i32;
       let g_2 = v_r_s as i32 - u_r_e as i32;
       (g_1 - g_2).abs() as usize > g
    }

    pub fn rc_gap_too_long(&self, u_r_s: usize, u_q_e: usize, v_q_s: usize, v_r_e: usize, g: usize) -> bool {
        //(self.rc_gap_start_too_long(u_q_s, u_r_s, v_q_s, v_r_s) ||self.rc_gap_end_too_long(u_q_e, u_r_e, v_q_e, v_r_e))
        let g_1 = v_q_s as i32 - u_q_e as i32;
        let g_2 = u_r_s as i32- v_r_e as i32;
        (g_1 - g_2).abs() as usize > g
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
    
    // Extends the first and last Hit locations in the Chain to the length of the query.
    pub fn find_coords(&self, rc: bool, r_id: &str, r_len: usize, q_id: &str, q_len: usize, mapq: usize, first: &Hit, last: &Hit, score: usize) -> Match {
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
        let mut len = self.len();
        if len > 1 {
            self.retain_unique();
            self.filter_hits(params.g);
        }
        let len_f = self.len();
        if len_f == 0 {return None;}
        let score = self.get_count(); 
        let mapq = match (params.s != 0 && params.c != 0) && ((len_f >= params.c) || (score >= params.s)) {
            true => 60,
            false => 0,
        };
        Some(self.find_coords(self.first().rc, r_id, r_len, q_id, q_len, mapq, self.first(), self.last(), score))
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
