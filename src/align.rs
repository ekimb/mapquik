// align.rs
// Contains Smith-Waterman, WFA, WFA2, and wflambda alignment functions (for two byte string slices),
// along with auxiliary functions.


use crate::mers::{AlignCand, Offset};
use dashmap::DashMap;
use bio::alphabets::dna;
use std::borrow::Cow;
use bio::alignment::pairwise::*;
use crate::cigar;
use block_aligner::scan_block::*;
use block_aligner::scores::*;
use block_aligner::cigar::Cigar;

/*
use rust_wfa2::aligner::*;
use wfa::wflambda::wf_align as wflambda_align;
use wfa::wfa::wf_align;
use libwfa::{affine_wavefront::*, bindings::*, mm_allocator::*, penalties::*};
*/


// Obtain reference locations that need to be aligned, and generate pointers to those string slices. These will then be aligned when the queries are parsed for a second time.
// Output: coordinate locations that need to be aligned in aln_coords (for the reference) and
//         aln_coords_q (for the query); pointers to the string slices of the reference that need to be aligned,
//         storing them in aln_seqs_cow. 
// Usage: in closures.rs in the closure 'ref_process_read_aux_aln'.
pub fn get_slices(ref_id: usize, ref_str: &[u8], aln_coords: &DashMap<usize, Vec<AlignCand>>, aln_coords_q:&DashMap<String, Vec<Offset>>, aln_seqs_cow: &DashMap<(String, (usize, usize)), (Cow<[u8]>, bool, usize, usize)> ) {
    let aln_coord_v = aln_coords.get(&ref_id).unwrap();
    for (r_tup, q_id, q_tup, rc) in aln_coord_v.iter() {
        let s = ref_str[r_tup.0..r_tup.1].to_vec();
        aln_coords_q.get_mut(q_id).unwrap().push(*q_tup);
        aln_seqs_cow.insert((q_id.to_string(), *q_tup), (Cow::from(s), *rc, ref_id, r_tup.0));
    }
}

pub struct AlignStats {
    pub successful : u64,
    pub failed: u64
}

// Retrieve query locations and pointer to reference string, and run alignment on query and string slice (currently WFA from libwfa).
// Output: an AlignStats struct which  contains the number of successful and failed alignments. 
// Usage: in closures.rs in the closure 'query_process_read_aux_aln'.
pub fn align_slices(seq_id: &str, seq_str: &[u8], aln_coords_q: &DashMap<String, Vec<Offset>>, aln_seqs_cow: &DashMap<(String, (usize, usize)), (Cow<[u8]>, bool, usize, usize)>, ref_map: &DashMap<usize, (String,usize)>) -> (Option<String>, Option<AlignStats>) {
    //let mode = 0; // rust-bio SW
    //let mode = 1; // wfa1
    //let mode = 2; // wfa2
    let mode = 3; // block-aligner 
    
    // 1) Failed attempt to use Options to enable either WFA1 or WFA2. Didnt work, due to even having
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
    
    // 2) So instead I'm just going to comment/uncomment here:
    
    // wfa1
    //let wfa1_alloc = libwfa::mm_allocator::MMAllocator::new(BUFFER_SIZE_512M as u64); 
    /*
    let mut penalties = AffinePenalties {
        match_: 0,
        mismatch: 4,
        gap_opening: 6,
        gap_extension: 2,
    };*/

    // wfa2
    //let mut wfa2_aligner = WFAlignerGapAffine::new(4, 6, 2, AlignmentScope::Alignment, MemoryModel::MemoryHigh); // just *uncommenting* this line, which supposedly does not much, is enough to crash wfa1 (try it!)
    //wfa2_aligner.set_heuristic(Heuristic::None);
    //wfa2_aligner.set_heuristic(Heuristic::BandedAdaptive(-10, 10, 1));
    
    let mut align_stats = AlignStats { successful: 0, failed: 0};
    
    let aln_coords_v = aln_coords_q.get(seq_id).unwrap();
    let mut cigars = Vec::new();
    let mut pos : Option<usize> = None;
    let mut end_pos = 0;
    let mut ref_idx : Option<usize> = None;
    let mut rc : Option<bool> = None;
    let seq_str_rev =  dna::revcomp(seq_str); 
    let seq_len =  seq_str.len(); 
    for q_tup in aln_coords_v.iter() {
        let entry = aln_seqs_cow.get(&(seq_id.to_string(), *q_tup)).unwrap();
        let r = &entry.0;
        let _rc = entry.1;
        if rc.is_none() {
            rc = Some(_rc); 
        }
        let _ref_idx = entry.2;
        let r_pos = entry.3;
        let q :&[u8] = if _rc { &seq_str_rev[(seq_len-1-q_tup.1)..(seq_len-1-q_tup.0)] } else { &seq_str[q_tup.0..q_tup.1] } ;
        if pos.is_none() { pos = Some(r_pos+1); }
        else { pos = Some(std::cmp::min(r_pos+1,pos.unwrap())); }
        if ref_idx.is_none() { ref_idx = Some(_ref_idx); }
        end_pos = std::cmp::max(end_pos,r_pos+r.len());
        let (score, cigar) = match mode {
            0 => {align_stats.successful += 1; sw(q, &r)}, 
            1 => 
                // change this when wfa1 is enabled
                //libwfa(q, &r, wfa1_alloc.as_ref().unwrap(), &mut penalties, &mut align_stats)
                //libwfa(q, &r, &wfa1_alloc, &mut penalties, &mut align_stats),
                (0, String::new()),
            2 => 
                // change this when wfa2 is enabled
                //wfa2(q, &r, &mut wfa2_aligner.as_mut().unwrap(), &mut align_stats)
                //wfa2(q, &r, &mut wfa2_aligner, &mut align_stats)
                (0, String::new()),
            3 => { align_stats.successful += 1; block_aligner(q,&r) }
            _ =>  (0, String::new())
        };
        //println!("{}\t{}\t{}\t{}\t{}\t{}", seq_id, q_tup.0, q_tup.1, r.len(), seq_str.len(), cigar);
        cigars.push(cigar);
        // wflambda(q, &r);
    }
    let mut sam_line : Option<String> = None;
    if cigars.len() > 0 {
        let rc = rc.unwrap();
        let flag = if rc { 16 } else {0};
        let mapq = 60;
        let pos = pos.unwrap();
        let tlen = end_pos-pos;
        let ref_idx = ref_idx.unwrap();
        let ref_id = &ref_map.get(&ref_idx).unwrap().0;
        let cigar = cigar::merge_cigar_strings(cigars);
        // SAM output
        sam_line = Some(format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", seq_id, flag, ref_id, pos, mapq, cigar, "*", "*", tlen, "*", "*"));
    }
    (sam_line, Some(align_stats))
}

// WFA implementation from https://github.com/chfi/rs-wfa
/*
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
*/

// Smith-Waterman alignment from rust-bio, mainly for debugging.

pub fn sw(q: &[u8], r: &[u8]) -> (i32, String) {
    let score = |a: u8, b: u8| if a == b { 2i32 } else { -4i32 };
    let mut aligner = Aligner::with_capacity(q.len(), r.len(), -4, -2, &score);
    let alignment = aligner.semiglobal(q, r);
    return (alignment.score, alignment.cigar(false));
}

// block_aligner
pub fn block_aligner(q: &[u8], r: &[u8]) -> (i32, String) {
    let block_size = 16;
    let run_gaps = Gaps { open: -4, extend: -2 };
    let r_padded = PaddedBytes::from_bytes::<NucMatrix>(r, block_size);
    let q_padded = PaddedBytes::from_bytes::<NucMatrix>(q, block_size);
    let mut block_aligner = Block::<true, false>::new(q.len(), r.len(), block_size);
    block_aligner.align(&q_padded, &r_padded, &NW1, run_gaps, block_size..=block_size, 0);
    let res = block_aligner.res();
    let block_score = res.score as u32;
    let mut cigar = Cigar::new(res.query_idx, res.reference_idx);
    // Compute traceback and resolve =/X (matches/mismatches).
    block_aligner.trace().cigar_eq(&q_padded, &r_padded, res.query_idx, res.reference_idx, &mut cigar);
    (block_score as i32, cigar.to_string())
    //(0,String::new())
}



