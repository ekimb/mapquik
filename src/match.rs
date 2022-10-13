// match.rs
// Contains the "Match" struct, which represents a collection of consecutive k-min-mer matches from query to reference.

use crate::{Entry, Index, ReadOnlyIndex};
use std::{cmp, fmt, iter::Peekable};
use rust_seq2kminmers::{KminmersIterator, KminmerType as Kminmer, Kminmer as KminmerTrait};


#[derive(Clone, Debug, PartialEq)]
pub struct Match {
    pub q_start: usize, // Query start location
    pub q_end: usize, // Query end location
    pub r_start: usize, // Reference start location
    pub r_end: usize, // Reference end location
    pub count: usize, // Number of k-min-mer matches
    pub rc: bool, // Strand direction
}
impl Match {
    
    // A new Match object from a query Kminmer matching a reference Entry.
    pub fn new(q: &Kminmer, r: &Entry) -> Self {
        Match {
            q_start: q.start,
            q_end: q.end,
            r_start: r.start,
            r_end: r.end,
            count: 1,
            rc: (q.rev != r.rc),
        }
    }

    // An empty Match object.
    pub fn empty() -> Self {
        Match {
            q_start: 0,
            q_end: 0,
            r_start: 0,
            r_end: 0,
            count: 0,
            rc: false,
        }
    }

    // Update the Match start and end locations, and count.
    pub fn update(&mut self, q: &Kminmer, r: &Entry) {
        if self.rc {self.r_start = r.start;}
        else {self.r_end = r.end;}
        self.q_end = q.end;
        self.count += 1;

    }

    // Check if this Match can be extended by another query Kminmer matching a new reference Entry. 
    pub fn check(&self, q: &Kminmer, r: &Entry, p: &Entry) -> bool {
        ((r.id == p.id) && ((q.rev != r.rc) == self.rc) && 
        (self.rc && (p.offset as i32 - r.offset as i32 == 1)) || 
        (!self.rc && (r.offset as i32 - p.offset as i32 == 1)))
    }

    // Extend this Match if it can be extended by the next Kminmer match.
    pub fn extend(&mut self, query_it: &mut Peekable<&mut KminmersIterator>, index: &ReadOnlyIndex, p: &Entry) {
        if let Some(q) = query_it.peek() {
            let re = index.get(&q.get_hash());
            if let Some(r) = re {
                if self.check(&q, &r, p) {
                    self.update(&q, &r);
                    query_it.next();
                    self.extend(query_it, index, &r)
                }
            }
            else {query_it.next();}
        }
        else {query_it.next();}
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
impl fmt::Display for Match {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let rc = match self.rc {
            true => "-",
            false => "+",
        };
        write!(f, "S!{}!E!{}!SP!{}!S!{}!E!{}!SP!{}!SC!{}!RC!{}!\n", self.q_start, self.q_end, self.q_span(), self.r_start, self.r_end, self.r_span(), self.count, rc)
    }
}


