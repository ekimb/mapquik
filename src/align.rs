// align.rs
// Contains Smith-Waterman, WFA, WFA2, and wflambda alignment functions (for two byte string slices),
// along with auxiliary functions.


use crate::mers::{AlignCand, Offset};
use dashmap::DashMap;
use bio::alphabets::dna;
use std::borrow::Cow;
use bio::alignment::pairwise::*;

use rust_wfa2::aligner::*;
use wfa::wflambda::wf_align as wflambda_align;
use wfa::wfa::wf_align;
use libwfa::{affine_wavefront::*, bindings::*, mm_allocator::*, penalties::*};



// Obtain reference locations that need to be aligned, and generate pointers to those string slices. These will then be aligned when the queries are parsed for a second time.
pub fn get_slices(ref_id: &str, ref_str: &[u8], aln_coords: &DashMap<String, Vec<AlignCand>>, aln_coords_q:&DashMap<String, Vec<Offset>>, aln_seqs_cow: &DashMap<(String, (usize, usize)), Cow<[u8]>> ) {
    let aln_coord_v = aln_coords.get(ref_id).unwrap();
    for (r_tup, q_id, q_tup, rc) in aln_coord_v.iter() {
        let mut s = ref_str[r_tup.0..r_tup.1].to_vec();
        aln_coords_q.get_mut(q_id).unwrap().push(*q_tup);
        if *rc {
            aln_seqs_cow.insert((q_id.to_string(), *q_tup), Cow::from(dna::revcomp(&s)));
        }
        else {aln_seqs_cow.insert((q_id.to_string(), *q_tup), Cow::from(s));}
        //println!("{}\t{}", score, cigar);
    }
}

pub struct AlignStats {
    pub successful : u64,
    pub failed: u64
}

// Retrieve query locations and pointer to reference string, and run alignment on query and string slice (currently WFA from libwfa).
pub fn align_slices(seq_id: &str, seq_str: &[u8], aln_coords_q: &DashMap<String, Vec<Offset>>, aln_seqs_cow: &DashMap<(String, (usize, usize)), Cow<[u8]>>) -> AlignStats {
    // Not sure if there is a better buffer size to lower max RSS / runtime
    let mut penalties = AffinePenalties {
        match_: 0,
        mismatch: 4,
        gap_opening: 6,
        gap_extension: 2,
    };
    
    let mode = 1; // wfa1
    //let mode = 2; // wfa2
    
    // Failed attempt to use Options to enable either WFA1 or WFA2. Didnt work, due to even having
    // both codes at the same time crash wfa1.
    // --
    // WFA1 and WFA2 don't play nice together. Seems one cannot have both allocators declared for both libs 
    // at the same time.
    //let mut wfa1_alloc : Option<libwfa::mm_allocator::MMAllocator> = None;
    //let mut wfa2_aligner : Option <rust_wfa2::aligner::WFAligner> = None; // just *uncommenting* this line, which supposedly does not much, is enough to crash wfa1 (try it!)
    //if mode == 1
    //{ 
        //wfa1_alloc = Some(libwfa::mm_allocator::MMAllocator::new(BUFFER_SIZE_512M as u64)); 
    //}
    //else
    //{
        //let mut wfa2_aligner = Some(WFAlignerGapAffine::new(4, 6, 2, AlignmentScope::Alignment, MemoryModel::MemoryHigh));
        //wfa2_aligner.set_heuristic(Heuristic::None);
        //wfa2_aligner.as_mut().unwrap().set_heuristic(Heuristic::BandedAdaptive(-10, 10, 1));
        //wfa2_aligner.set_heuristic(Heuristic::BandedAdaptive(-10, 10, 1));
    //}
    
    let wfa1_alloc = libwfa::mm_allocator::MMAllocator::new(BUFFER_SIZE_512M as u64); 
    //let mut wfa2_aligner = WFAlignerGapAffine::new(4, 6, 2, AlignmentScope::Alignment, MemoryModel::MemoryHigh); // just *uncommenting* this line, which supposedly does not much, is enough to crash wfa1 (try it!)
    //wfa2_aligner.set_heuristic(Heuristic::None);
    //wfa2_aligner.set_heuristic(Heuristic::BandedAdaptive(-10, 10, 1));
    
    let mut align_stats = AlignStats { successful: 0, failed: 0};
    
    let aln_coords_v = aln_coords_q.get(seq_id).unwrap();
    for q_tup in aln_coords_v.iter() {
        let mut r = aln_seqs_cow.get(&(seq_id.to_string(), *q_tup)).unwrap();
        //println!("{}\t{}\t{}\t{}", seq_id, q_tup.0, q_tup.1, r.len());
        let q = &seq_str[q_tup.0..q_tup.1];
        let (score, cigar) = if mode == 1 {
                //libwfa(q, &r, wfa1_alloc.as_ref().unwrap(), &mut penalties, &mut align_stats)
                libwfa(q, &r, &wfa1_alloc, &mut penalties, &mut align_stats)
                //(0, String::new())
        }
        else
        {   
                // change this when wfa2 is enabled
                (0, String::new())
                //wfa2(q, &r, &mut wfa2_aligner.as_mut().unwrap(), &mut align_stats)
                //wfa2(q, &r, &mut wfa2_aligner, &mut align_stats)
        };
        println!("{}!{}!{}!\t{}\t{}", seq_id, q_tup.0, q_tup.1, score, cigar);
        // wflambda(q, &r);
        // sw(q, &r);        
    }
    align_stats
}

// WFA implementation from https://github.com/chfi/rs-wfa.
pub fn libwfa(q_str: &[u8], ref_str: &[u8], alloc: &MMAllocator, penalties: &mut AffinePenalties, align_stats : &mut AlignStats) -> (isize, String) {
    let pat_len = q_str.len();
    let text_len = ref_str.len();

    /*
    let mut wavefronts = AffineWavefronts::new_complete(
        pat_len,
        text_len,
        &mut penalties,
        &alloc,
    );
    */
    
    //Also consider adaptive-WFA:
    let mut wavefronts = AffineWavefronts::new_reduced(
        pat_len,
        text_len,
        penalties,
        5,
        25,
        &alloc,
    ); 

    let result = wavefronts.align(q_str, ref_str);

    match result 
    {
        Ok(..) => align_stats.successful += 1,
        _ => align_stats.failed += 1
    }

    let score = wavefronts.edit_cigar_score(penalties);
    let cigar = wavefronts.cigar_bytes();
    let cg_str = std::str::from_utf8(&cigar).unwrap();
    //let score = 0;
    //let cg_str = String::new();
    (score, cg_str.to_string())
}

// WFA2 from https://github.com/tanghaibao/rust-wfa2/
pub fn wfa2(q_str: &[u8], ref_str: &[u8], wfa2_aligner: &mut WFAligner, align_stats : &mut AlignStats) -> (isize, String) {

    let status = wfa2_aligner.align_end_to_end(q_str,ref_str);
    match status
    {
         AlignmentStatus::StatusSuccessful => align_stats.successful += 1,
        _ => align_stats.failed += 1
    }
    //println!("status: {}",status as i32);

    let score = wfa2_aligner.score();
    //let cigar = aligner.cigar();
    let cigar = String::new();
    (score as isize, cigar)
}



// TEST_CONFIG type necessary for @urbanslug's Rust WFA/wflambda implementation.
pub static TEST_CONFIG: wfa::types::Config = wfa::types::Config {
    adapt: false,
    verbosity: 0,
    penalties: wfa::types::Penalties {
        mismatch: 4,
        matches: 0,
        gap_open: 6,
        gap_extend: 2,
    },
};

// wflambda alignment from https://github.com/urbanslug/wfa.

pub fn wflambda(q: &[u8], t: &[u8]) -> (usize, String) {
    let tlen = t.len();
    let qlen = q.len();
    //println!("{:?}, {:?}", q.len(), t.len());
    let mut match_lambda = |v: &mut i32,  h: &mut i32, offset: &mut i32| {
        if *v < 0 || *h < 0 {return false;}
        let v_idx = *v as usize;
        let h_idx = *h as usize;
        let res = h_idx < tlen && v_idx < qlen && t[h_idx] == q[v_idx];
        if res {
            *v += 1;
            *h += 1;
            *offset += 1;
        }
        res
    };
    let mut traceback_lambda = |(q_start, q_stop): (i32, i32), (t_start, t_stop): (i32, i32)| -> bool {
        if q_start < 0
            || q_stop as usize > qlen
            || t_start < 0
            || t_stop as usize > tlen {
            return false;
        }

        let q_start = q_start as usize;
        let q_stop = q_stop as usize;
        let t_start = t_start as usize;
        let t_stop = t_stop as usize;

        q[q_start..q_stop] == t[t_start..t_stop]
    };
    let (score, cigar) = wflambda_align(
        tlen as u32,
        qlen as u32,
        &TEST_CONFIG,
        &mut match_lambda,
        &mut traceback_lambda
    ).unwrap();
    (score, cigar)
}

// Smith-Waterman alignment from rust-bio, mainly for debugging.

pub fn sw(q: &[u8], r: &[u8]) -> (i32, String) {
    let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
    let mut aligner = Aligner::with_capacity(q.len(), r.len(), -5, -1, &score);
    let alignment = aligner.semiglobal(q, r);
    return (alignment.score, alignment.cigar(false));
}






