// chain.rs
// Contains the "Chain" struct, which is a Vec of matches, and various necessary operations to filter out bad matches.

use crate::{r#match::Match, Params, PseudoChainCoords};
use std::fmt;

#[derive(Clone, Debug, PartialEq)]
pub struct Chain {
    matches: Vec<Match>,
}
impl Chain {

    // New Chain from &Vec of matches.
    pub fn new(matches: &[Match]) -> Self {
        Chain {matches: matches.to_vec()}
    }

    // Get total number of k-min-mer matches in the Chain (a Match can have multiple consecutive k-min-mer matches).
    pub fn get_count(&self) -> usize {
        return self.matches.iter().map(|h| h.count).sum::<usize>();
    }

    // Get the number of matches in the Chain.
    pub fn len(&self) -> usize {
        return self.matches.len()
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

    pub fn check_match_compatible(&self, h1: &Match, h2: &Match, g: usize) -> bool {
        if h1 == h2 {return true;}
        if h1.rc != h2.rc {return false;}
        let (u, v) = match h1.q_start < h2.q_start {
            true => (h1, h2),
            false => (h2, h1),
        };
        let u_q_e = u.q_end;
        let v_q_s = v.q_start;
        let u_r_s = u.r_start;
        let u_r_e = u.r_end;
        let v_r_s = v.r_start;
        let v_r_e = v.r_end;
        if u.rc {
            if u_r_s <= v_r_s || self.rc_gap_too_long(u_r_s, u_q_e, v_q_s, v_r_e, g) {return false;}
        }
        else if v_r_s <= u_r_s || self.fwd_gap_too_long(u_q_e, u_r_e, v_q_s, v_r_s, g) {
            return false;
        }
        return true;
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

    /*pub fn filter_matches(&mut self, g: usize) {
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
        self.matches = max_chain.into_iter().cloned().collect();
    }*/

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

    /*pub fn filter_matches_c(&mut self, g: usize, c: usize) {
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
    }*/

    pub fn filter_matches_max(&mut self, g: usize) {
        let len = self.len();
        if len <= 1 {return;}
        let max = self.find_largest_match();
        let (max_chain, _) = self.colinear_matches_per_match(self.nth(max), g);
        self.matches = max_chain.into_iter().cloned().collect();
    }


    pub fn fwd_gap_too_long(&self, u_q_e: usize, u_r_e: usize, v_q_s: usize, v_r_s: usize, g: usize) -> bool {
       let g_1 = v_q_s as i32 - u_q_e as i32;
       let g_2 = v_r_s as i32 - u_r_e as i32;
       (g_1 - g_2).abs() as usize > g
    }

    pub fn rc_gap_too_long(&self, u_r_s: usize, u_q_e: usize, v_q_s: usize, v_r_e: usize, g: usize) -> bool {
        let g_1 = v_q_s as i32 - u_q_e as i32;
        let g_2 = u_r_s as i32- v_r_e as i32;
        (g_1 - g_2).abs() as usize > g
    }

    // Wrapper function that filters bad matches, checks for consistency, and obtains final coordinates.

    // Outputs a Match object (see mers.rs for a definition).
    pub fn get_match(&mut self, params: &Params) -> Option<PseudoChainCoords> {
        let len = self.len();
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
