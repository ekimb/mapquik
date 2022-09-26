// hit.rs
// Contains the "Hit" struct, which represents a collection of consecutive k-min-mer matches from query to reference.

use crate::{Entry, Index, Kminmer, Params};
use std::cmp;
use std::fmt;

#[derive(Clone, Debug, PartialEq)]
pub struct Hit {
    pub q_id: String, // Query ID
    pub r_id: String, // Reference ID
    pub q_start: usize, // Query start location
    pub q_end: usize, // Query end location
    pub r_start: usize, // Reference start location
    pub r_end: usize, // Reference end location
    pub count: usize, // Number of k-min-mer matches
    pub score: usize, // Approximate number of bases matched
    pub rc: bool, // Strand direction
    pub q_offset: usize, // Query k-min-mer offset
    pub r_offset: usize, // Reference k-min-mer offset
}
impl Hit {
    
    // A new Hit object from a query Kminmer matching a reference Entry.
    pub fn new(q_id: &str, q: &Kminmer, r: &Entry, params: &Params) -> Self {
        let mut h = Hit {
            q_id: q_id.to_string(),
            r_id: r.id.to_string(),
            q_start: q.start,
            q_end: q.end,
            r_start: r.start,
            r_end: r.end,
            count: 1,
            score: 0,
            rc: (q.rev != r.rc),
            q_offset: q.offset,
            r_offset: r.offset,
        };
        let new_score = h.score_match(q, r, params);
        h.score = new_score;
        h

    }

    // An empty Hit object.
    pub fn empty() -> Self {
        Hit {
            q_id: String::new(),
            r_id: String::new(),
            q_start: 0,
            q_end: 0,
            r_start: 0,
            r_end: 0,
            count: 0,
            score: 0,
            rc: false,
            q_offset: 0,
            r_offset: 0,
        }
    }

    // Update the Hit start and end locations, count, and score.
    pub fn update(&mut self, q: &Kminmer, r: &Entry, params: &Params) {
        let new_score_add = self.score_extension(q, r, params);
        if self.rc {self.r_start = r.start;}
        else {self.r_end = r.end;}
        self.q_end = q.end;
        self.count += 1;
        self.score += new_score_add;

    }

    // Check if this Hit can be extended by another query Kminmer matching a new reference Entry. 
    pub fn check(&self, q: &Kminmer, r: &Entry, p: &Entry, q_len: usize) -> bool {
        ((r.id == self.r_id) && ((q.rev != r.rc) == self.rc) && 
        (self.rc && (self.r_start > r.start) && (self.r_start - r.start < q_len / 2)) || 
        (!self.rc && (self.r_end < r.end) && (r.end - self.r_end < q_len / 2)))
    }

    // Extend this Hit if it can be extended by the next Kminmer match.
    pub fn extend(&mut self, i: usize, query_mers: &Vec<Kminmer>, index: &Index, p: &Entry, params: &Params, q_len: usize) {
        if i == query_mers.len() - 1 {return;}
        let q = &query_mers[i + 1];
        let (b, r) = index.get_entry(&q);
        //println!("HIT!{}!!QS!{}!QE!{}!QOFF!{}!QRC!{}!RS!{}!RE!{}!ROFF!{}!RRC!{}!", self, q.start, q.end, q.offset, q.rev, r.start, r.end, r.offset, r.rc);
        if b && self.check(q, &r, p, q_len) {
            self.update(q, &r, params);
            self.extend(i + 1, query_mers, index, &r, params, q_len);
        }
    }
    
    // Calculate (approximately) the number of matching bases in this Kminmer match.
    pub fn score_match(&self, q: &Kminmer, r: &Entry, params: &Params) -> usize {
        cmp::min(cmp::min((q.end as i32 - q.start as i32).abs() as usize, (r.start as i32 - r.end as i32).abs() as usize), (params.k * params.l))
    }

    // Calculate (approximately) the number of matching bases if the Kminmer match is extended.
    pub fn score_extension(&self, q: &Kminmer, r: &Entry, params: &Params) -> usize {
        cmp::min(cmp::min((q.end as i32 - q.start as i32).abs() as usize, (r.start as i32 - r.end as i32).abs() as usize),  params.l)
    }

    // Calculate "offset difference" of this Hit:
    
    // An "offset" in the context of a Hit is the index of the first matched k-min-mer in the ordered k-min-mer array (of the query or the reference). An "offset difference" is the absolute difference between the indices of the query and reference k-min-mers in their respective arrays.
    pub fn offset_diff(&self) -> usize {
        (self.q_offset as i32 - self.r_offset as i32).abs() as usize
    }

    // Calculate the number of query bases covered by this Hit.
    pub fn q_span(&self) -> usize {
        (self.q_end as i32 - self.q_start as i32).abs() as usize
    }

    // Calculate the number of reference bases covered by this Hit.
    pub fn r_span(&self) -> usize {
        (self.r_end as i32 - self.r_start as i32).abs() as usize
    }

    pub fn span_diff(&self) -> usize {
        (self.q_span() as i32 - self.r_span() as i32).abs() as usize
    }
}

// Pretty-prints a Hit.
impl fmt::Display for Hit {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let rc = match self.rc {
            true => "-",
            false => "+",
        };
        write!(f, "{}!S!{}!E!{}!SP!{}!OFF!{}!{}!S!{}!E!{}!SP!{}!OFF!{}!CT!{}!SC!{}!OFFD!{}!RC!{}!\n", self.q_id, self.q_start, self.q_end, self.q_span(), self.q_offset, self.r_id, self.r_start, self.r_end, self.r_span(), self.r_offset, self.count, self.score, self.offset_diff(), rc)
    }
}


