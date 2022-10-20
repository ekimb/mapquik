// chain.rs
// Contains the "Chain" struct, which is a Vec of matches, and various necessary operations to filter out bad matches.

use crate::{Entry, Index, r#match::Match, Params, PseudoChainCoords, PseudoChainCoordsTuple};
use std::collections::HashMap;
use std::fmt;

#[derive(Clone, Debug, PartialEq)]
pub struct Chain {
    matches: Vec<Match>,
}
impl Chain {

    // New Chain from &Vec of matches.
    pub fn new(matches: &Vec<Match>) -> Self {
        Chain {matches: matches.to_vec()}
    }

    // An empty Chain object.
    pub fn empty() -> Self {
        Chain {matches: Vec::new()}
    }

    pub fn filter_match_span(&mut self) {
        self.matches.retain(|h| h.span_diff() < (h.r_span()));
    }

    // Get the total difference in span (difference between length of query covered by a Match and length of reference covered by the Match) of the Chain.
    pub fn span_diff(&self) -> usize {
        let mut s = 0;
        for i in 0..self.len() {
            let h = self.nth(i);
            s += (h.q_span() as i32 - h.r_span() as i32).abs() as usize;
        }
        s
    }

    // Get total number of k-min-mer matches in the Chain (a Match can have multiple consecutive k-min-mer matches).
    pub fn get_count(&self) -> usize {
        return self.matches.iter().map(|h| h.count).sum::<usize>();
    }

    // Get the number of matches in the Chain.
    pub fn len(&self) -> usize {
        return self.matches.len()
    } 

    // Check if the Chain is empty.
    pub fn is_empty(&self) -> bool {
        self.matches.is_empty()
    }

    // Check if the Chain contains the Match.
    pub fn contains(&self, h: &Match) -> bool {
        self.matches.contains(h)
    }

    // Return a raw Vec of matches.
    pub fn matches(&self) -> &Vec<Match> {
        &self.matches
    }

    // Get the first Match in the Chain.
    pub fn first(&self) -> &Match {
        &self.matches[0]
    }

    // Get the last Match in the Chain.
    pub fn last(&self) -> &Match {
        &self.matches[self.len() - 1]
    }

    // Get the nth Match in the Chain.
    pub fn nth(&self, i: usize) -> &Match {
        &self.matches[i]
    }

    // Remove the ith Match in the Chain.
    pub fn remove(&mut self, i: usize) {
        self.matches.remove(i);
    }

    // Reverse the order of matches in the Chain.
    pub fn reverse(&mut self) {
        self.matches.reverse();
    }

    pub fn check_match_compatible(&self, h1: &Match, h2: &Match, g: usize) -> bool {
        if h1 == h2 {return true;}
        if h1.rc != h2.rc {return false;}
        let (u, v) = match h1.q_start < h2.q_start {
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
        let mut remove_match_indices = vec![false; len];
        for i in 0..self.len() {
            let h_i = self.nth(i);
            let h_c = self.matches.iter().filter(|hp| hp.r_start == h_i.r_start && hp.r_end == h_i.r_end && hp.rc == h_i.rc).count();
            if h_c > 1 {remove_match_indices[i] = true;}
        }
        self.matches = (0..len).filter(|i| !remove_match_indices[*i]).map(|i| self.nth(i).clone()).collect::<Vec<Match>>();

    }

    pub fn colinear_idx_per_match(&self, h: &Match, g: usize) ->  Vec<usize> {
        let col = (0..self.len()).filter(|i| self.check_match_compatible(h, self.nth(*i), g)).collect::<Vec<usize>>();
        col
    }
    pub fn colinear_matches_per_match(&self, h: &Match, g: usize) ->  (Vec<&Match>, usize) {
        let mut col = Vec::new();
        let mut tot = 0;
        for hp in self.matches.iter() {
            if self.check_match_compatible(h, hp, g) {
                col.push(hp);
                tot += hp.count;
            }
        }
        (col, tot)
    }

    pub fn check_colinear_match_set(&self, col: &Vec<usize>, g: usize) -> (Vec<usize>, usize) {
        let mut remove_h = HashMap::<usize, Vec<usize>>::new();
        let mut remove_count : i32 = 0;
        let mut to_remove = vec![false; self.len()];

        for idx_i in 0..col.len() - 1 {
            let i = col[idx_i];
            let mut h_i = self.nth(i);
            for idx_j in idx_i..col.len() {
                let j = col[idx_j];
                let mut h_j = self.nth(j);
                if !self.check_match_compatible(h_i, h_j, g) {
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
            let col = self.colinear_idx_per_match(h_i, g);
            let (col_f, score) = self.check_colinear_match_set(&col, g);
            if score > max_score {
                max_col = col_f;
                max_score = score;
            }
        }
        self.matches = max_col.iter().map(|i| self.nth(*i).clone()).collect::<Vec<Match>>();
    }

    pub fn filter_matches(&mut self, g: usize) {

        let len = self.len();
        if len == 0 {return;}
        let mut max_chain = Vec::new();
        let mut max_count = 0;
        for i in 0..len {
            let h = self.nth(i);
            let (chain, count) = self.colinear_matches_per_match(h, g);
            if count > max_count {
                max_chain = chain;
                max_count = count;
            }
        }
        //let max_chain = self.colinear_matches_per_match(self.matches.iter().max_by(|a, b| self.colinear_matches_per_match(a, g).iter().map(|h| h.count).sum::<usize>().cmp(&self.colinear_matches_per_match(b, g).iter().map(|hp| hp.count).sum::<usize>())).unwrap(), g);
        self.matches = max_chain.into_iter().cloned().collect();
        //self.sort_by_q_start();
    }

    pub fn find_largest_two_matches(&self) -> (usize, usize) {
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

    pub fn find_largest_match(&self) -> usize {
        let mut max = 0;
        let mut max_count = 0;
        for i in 0..self.len() {
            let count = self.nth(i).count;
            if count > max_count {
                max = i;
                max_count = count;
            }
        }
        max
    }

    pub fn filter_matches_c(&mut self, g: usize, c: usize) {
        let len = self.len();
        let mut sorted_idx = (0..len).collect::<Vec<usize>>();
        sorted_idx.sort_by(|a, b| self.nth(*a).count.cmp(&self.nth(*b).count));
        let mut max_chain = Vec::new();
        let mut max_count = 0;
        for i in 0..c {
            let max_i = sorted_idx[len - 1 - i];
            let (chain, count) = self.colinear_matches_per_match(self.nth(max_i), g);
            if count > max_count {
                max_chain = chain;
                max_count = count;
            }
        }
        self.matches = max_chain.into_iter().cloned().collect();
    }

    pub fn filter_matches_max(&mut self, g: usize) {
        let len = self.len();
        if len <= 1 {return;}
        let max = self.find_largest_match();
        let (max_chain, tot) = self.colinear_matches_per_match(self.nth(max), g);
        self.matches = max_chain.into_iter().cloned().collect();
        //self.sort_by_q_start();
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
        (self.fwd_start_inconsistent(u_r_s, v_r_s)) //|| self.fwd_end_inconsistent(u_r_e, v_r_e))
    }
    pub fn rc_inconsistent(&self, u_r_s: usize, v_r_s: usize, u_r_e: usize, v_r_e: usize) -> bool {
        (self.rc_start_inconsistent(u_r_s, v_r_s)) //|| self.rc_end_inconsistent(u_r_e, v_r_e))
    }

    // Sort the Chain by the query k-min-mer offsets of the matches.
    pub fn sort_by_q_start(&mut self) {
        self.matches.sort_by(|a, b| a.q_start.cmp(&b.q_start));
    }

    // Generate a new Chain that has no elements that also occur in another Chain c.
    pub fn find_complement(&self, c: &Chain) -> Chain {
        let mut c_p = self.matches.iter().filter(|h| !c.contains(h)).cloned().collect::<Vec<Match>>();
        Chain::new(&c_p)
    }

    // Modify this particular Chain so that it has no elements that also occur in another Chain c.
    pub fn retain_complement(&mut self, c: &Chain) {
        self.matches.retain(|h| !c.contains(h));
    }

    // Replace this Chain object with another Chain c.
    pub fn replace(&mut self, c: &Chain) {
        self.matches = c.matches().to_vec();
    }

    // Wrapper function that filters bad matches, checks for consistency, and obtains final coordinates.

    // Outputs a Match object (see mers.rs for a definition).
    pub fn get_match(&mut self, params: &Params) -> Option<PseudoChainCoords> {
        let mut len = self.len();
        if len > 1 {
            //self.retain_unique();
            //self.filter_matches_c(params.g, len - 1);
            self.filter_matches_max(params.g);
            //self.check_colinear(params.g);
        }
        let len_f = self.len();
        if len_f == 0 {return None;}
        let score = self.get_count(); 
        let mapq = match (params.s != 0 && params.c != 0) && ((len_f >= params.c) || (score >= params.s)) {
            true => 60,
            false => 0,
        };
        let first = self.first();
        let last = self.last();
        let rc = first.rc;
        return match rc && self.len() > 1 {
            true => Some((rc, first.q_start, last.q_end - 1, last.r_start, first.r_end - 1, score, mapq)),
            false => Some((rc, first.q_start, last.q_end - 1, first.r_start, last.r_end - 1, score, mapq)),
        };
    }
    

    // Obtains query and reference intervals that are not covered by a Match in the final Chain object (only for base-level alignment).
    /*pub fn get_remaining_seqs(&self, m: &Match) -> (Vec<(usize, usize)>, Vec<(usize, usize)>) {
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
                let match = self.nth(i);
                let rc_s = match match.rc {
                    true => "-",
                    false => "+",
                };
                //println!("QALN!{}!{}!RALN!{}!{}", prev_q_end, match.q_start, prev_r_end, match.r_start);
                //println!("QMMM!{}!{}!RMMM!{}!{}", match.q_start, match.q_end, match.r_start, match.r_end); 
                if prev_q_end < &match.q_start && &match.r_start > prev_r_end {
                    q_coords.push((*prev_q_end, match.q_start));
                    r_coords.push((*prev_r_end, match.r_start));
                    prev_q_end = &match.q_end;
                    prev_r_end = &match.r_end;
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
                let match = self.nth(i);
                let rc_s = match match.rc {
                    true => "-",
                    false => "+",
                };
                if match.r_end < *r_start {break;}
                //println!("QALN!{}!{}!RALN!{}!{}", prev_q_end, match.q_start, prev_r_end, match.r_start);
                //println!("QMMM!{}!{}!RMMM!{}!{}", match.q_start, match.q_end, match.r_start, match.r_end); 
                if prev_q_end < &match.q_start && &match.r_end < prev_r_start {
                    q_coords.push((*prev_q_end, match.q_start));
                    r_coords.push((match.r_end, *prev_r_start));
                    prev_q_end = &match.q_end;
                    prev_r_start = &match.r_start;
                }
            }
            //println!("QALN!{}!{}!RALN!{}!{}", prev_q_end, q_end, r_start, prev_r_start);
            if prev_q_end != q_end {
                q_coords.push((*prev_q_end, *q_end));
                r_coords.push((*r_start, *prev_r_start));
            }
        }
        (q_coords, r_coords)
    } */
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
