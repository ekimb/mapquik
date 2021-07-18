use std::io::stderr;
use std::error::Error;
use std::io::Write;
use std::io::{BufWriter, BufRead, BufReader};
use std::collections::HashMap;
use std::collections::hash_map::Entry;
use itertools::Itertools;
use ::closure;
use std::iter::FromIterator;
use std::thread;
use crate::read::hash_id;
use std::fs::{File,remove_file};
use std::collections::HashSet;
use std::fs;
use std::sync::{Arc, Mutex, MutexGuard};
use seq_io::fasta;
use seq_io::core::BufReader as OtherBufReader;
use seq_io::BaseRecord;
use seq_io::parallel::{read_process_fasta_records, read_process_fastq_records};
use dashmap::DashMap;
use thread_id;
use std::io::Result;
use super::strobes;
use super::strobes::{NAM as StrobeNAM};
use super::kstrobes;
use super::kstrobes::{NAM as KStrobeNAM};
use super::kminmers;
use super::kminmers::{NAM as KminmerNAM};
use super::ksyncmers;
use super::ksyncmers::{NAM as KsyncmerNAM};
use super::{MersVector, MersVectorReduced, MersVectorRead, PosIndex, KmerLookup, KmerLookupMod, PosIndexMod};
use std::path::PathBuf;
use dashmap::DashSet;
use super::Params;
use crate::get_reader;

pub fn run_kstrobes(filename: &PathBuf, ref_filename: &PathBuf, params: &Params, threads: usize, queue_len: usize, fasta_reads: bool, ref_fasta_reads: bool, output_prefix: &PathBuf) {
    let mut id_hashes : Arc<DashMap<u64, String>> = Arc::new(DashMap::new());
    let mut id_lengths : Arc<DashMap<u64, usize>> = Arc::new(DashMap::new());
    let mut hash_counts : Arc<DashMap<u64, usize>> = Arc::new(DashMap::new());
    let mut tmp_index : Arc<PosIndex> = Arc::new(PosIndex::new());
    let mut mers_index : Arc<KmerLookupMod> = Arc::new(KmerLookupMod::new());
    let mut filter_cutoff : Arc<DashMap<usize, usize>> = Arc::new(DashMap::new());
    let mut all_nams = Vec::<(Vec<KStrobeNAM>, u64)>::new();
    let path = format!("{}{}", output_prefix.to_str().unwrap(), ".paf");
    let mut paf_file = match File::create(&path) {
        Err(why) => panic!("Couldn't create {}: {}", path, why.description()),
        Ok(paf_file) => paf_file,
    };
    let mut index_kstrobes = |seq_id_hash: u64, seq_id: &str, seq: &[u8], params: &Params, read: bool| -> MersVectorRead {
        let mut kstrobes = kstrobes::index(seq, seq_id_hash, params, read);
        //println!("Seq: {:?}, read_seq: {:?}", seq, read_seq);
        id_hashes.insert(seq_id_hash, seq_id.to_string());
        id_lengths.insert(seq_id_hash, seq.len());
        tmp_index.insert(seq_id_hash, kstrobes.to_vec()); 
        return kstrobes;   
    };
    let ref_process_read_aux_kstrobe = |ref_str: &[u8], ref_id: &str| -> Option<(MersVectorRead, u64)> {
        let thread_id :usize =  thread_id::get();
        //let seq = std::str::from_utf8(ref_str).unwrap().replace("\n", "").replace("\r", "");
        //let seq = seq_for_ref; 
        let ref_id_hash = hash_id(&ref_id);
        let mut kstrobes = index_kstrobes(ref_id_hash, ref_id, ref_str, params, false);
        return Some((kstrobes, ref_id_hash))
    
    };
    let ref_process_read_fasta_kstrobe = |record: seq_io::fasta::RefRecord, found: &mut Option<(MersVectorRead, u64)>| {
        let ref_str = record.seq(); 
        let ref_id = record.id().unwrap().to_string();
        *found = ref_process_read_aux_kstrobe(&ref_str, &ref_id);
    };
    let ref_process_read_fastq_kstrobe = |record: seq_io::fastq::RefRecord, found: &mut Option<(MersVectorRead, u64)>| {
        let ref_str = record.seq(); 
        let ref_id = record.id().unwrap().to_string();
        *found = ref_process_read_aux_kstrobe(&ref_str, &ref_id);
    };
    let mut ref_main_thread_kstrobe = |found: &mut Option<(MersVectorRead, u64)>| { // runs in main thread
        let (kstrobes, seq_id_hash) = found.as_ref().unwrap();
        None::<()>
    };
    let query_process_read_aux_kstrobe = |seq_str: &[u8], seq_id: &str| -> Option<(Vec<KStrobeNAM>, u64)> {
        let thread_id :usize =  thread_id::get();
        let mut output : Vec<KStrobeNAM> = Vec::new();
        //let seq = std::str::from_utf8(seq_str).unwrap(); // see https://docs.rs/seq_io/0.4.0-alpha.0/seq_io/fasta/index.html
        let seq_id_hash = hash_id(&seq_id);
        let mut kstrobes = index_kstrobes(seq_id_hash, seq_id, seq_str, params, true);
        let fc = 0; //*filter_cutoff.get(&0).unwrap();
        output = kstrobes::find_nams(&kstrobes, &tmp_index, &mers_index, &hash_counts, params.l, seq_str, fc, params.k, &id_hashes, &id_lengths);
        if output.len() > 0 {
            return Some((output, seq_id_hash))
        }
        else {return None}
    };
    
    let query_process_read_fasta_kstrobe = |record: seq_io::fasta::RefRecord, found: &mut Option<(Vec<KStrobeNAM>, u64)>| {
        let seq_str = record.seq(); 
        let seq_id = record.id().unwrap().to_string();
        *found = query_process_read_aux_kstrobe(&seq_str, &seq_id);
    
    };
    let query_process_read_fastq_kstrobe = |record: seq_io::fastq::RefRecord, found: &mut Option<(Vec<KStrobeNAM>, u64)>| {
        let seq_str = record.seq(); 
        let seq_id = record.id().unwrap().to_string();
        *found = query_process_read_aux_kstrobe(&seq_str, &seq_id);
    };
    let mut main_thread_kstrobe = |found: &Option<(Vec<KStrobeNAM>, u64)>| { // runs in main thread
        if found.is_some() {
            let (nams, seq_id_hash) = found.as_ref().unwrap();
            all_nams.push((nams.to_vec(), *seq_id_hash));
        }
        // Some(value) will stop the reader, and the value will be returned.
        // In the case of never stopping, we need to give the compiler a hint about the
        // type parameter, thus the special 'turbofish' notation is needed.
        None::<()>
    };
    let mut construct_flat_vector = || -> (MersVectorRead, u64) {
        let mut flat_vector = MersVectorRead::new();
        for entry in tmp_index.iter_mut() {
            for t in entry.iter() {flat_vector.push(*t);}
        }
        flat_vector.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        let mut prev_k = flat_vector[0].0;
        let mut unique_elements = 1;
        for entry in flat_vector.iter() {
            let mut curr_k = entry.0;
            if curr_k != prev_k {unique_elements += 1;}
            prev_k = curr_k;
        }
        return (flat_vector, unique_elements);
    };

    let mut index_vector = |tmp_index: &PosIndex| {
        for entry in tmp_index.into_iter() {
            let ref_id = entry.key();
            let flat_vector = entry.value();
            println!("Flat vector size: {}", flat_vector.len());
            let mut offset = 0;
            for entry in flat_vector.iter() {
                let curr_k = entry.0;
                let mut s = (*ref_id, offset);
                mers_index.insert(curr_k, s);
                offset += 1;      
            }
        }
    };

    let mut count_vector = |flat_vector: &MersVectorRead, f: f64| -> usize {
        println!("Flat vector size: {}", flat_vector.len());
        let mut offset = 0;
        let mut prev_offset = 0;
        let mut count = 0;
        let mut tot_occur_once = 0;
        let mut tot_high_ab = 0;
        let mut tot_mid_ab = 0;
        let mut kstrobes_counts = Vec::<usize>::new();
        let mut t = flat_vector[0];
        let mut prev_k = t.0;
        let mut curr_k = 0;
        for entry in flat_vector.iter() {
            curr_k = entry.0;
            if curr_k == prev_k {count += 1;}
            else {
                if count == 1 {tot_occur_once += 1;}
                else if count > 100 {tot_high_ab += 1;}
                else {tot_mid_ab += 1;}
                kstrobes_counts.push(count);
                hash_counts.insert(prev_k, count);
                count = 1;
                prev_k = curr_k;
                prev_offset = offset;
            }
            offset += 1;      
        }
        if count == 1 {tot_occur_once += 1;}
        else if count > 100 {tot_high_ab += 1;}
        else {tot_mid_ab += 1;}
        kstrobes_counts.push(count);
        hash_counts.insert(curr_k, count);
        println!("Total k-strobemers count:\t{}", offset);
        println!("Total k-strobemers unique:\t{}", tot_occur_once);
        println!("Total k-strobemers high-ab:\t{}", tot_high_ab);
        println!("Total k-strobemers mid-ab:\t{}", tot_mid_ab);
        if params.f != 0.0 {
            kstrobes_counts.sort_unstable();
            kstrobes_counts.reverse();
            let mut index_cutoff = (kstrobes_counts.len() as f64 * f) as usize;
            println!("Filtered cutoff index:\t{}", index_cutoff);
            let mut filter_cutoff = kstrobes_counts[index_cutoff];
            println!("Filtered cutoff count:\t{}", filter_cutoff);
            return filter_cutoff;
        }
        else {return 0;}
    };
    let buf = get_reader(&ref_filename);
    if ref_fasta_reads {
        let reader = seq_io::fasta::Reader::new(buf);
        read_process_fasta_records(reader, threads as u32, queue_len, ref_process_read_fasta_kstrobe, |record, found| {ref_main_thread_kstrobe(found)});
        let (mut flat_vector, unique_elements) = construct_flat_vector();
        println!("Unique k-strobemers:\t{}", unique_elements);
        let f = params.f;
        index_vector(&tmp_index);
        filter_cutoff.insert(0, count_vector(&flat_vector, f));

    }
    else {
        let reader = seq_io::fastq::Reader::new(buf);
        read_process_fastq_records(reader, threads as u32, queue_len, ref_process_read_fastq_kstrobe, |record, found| {ref_main_thread_kstrobe(found)});
        let (mut flat_vector, unique_elements) = construct_flat_vector();
        println!("Unique k-strobemers:\t{}", unique_elements);
        let f = params.f;
        index_vector(&tmp_index);
        filter_cutoff.insert(0, count_vector(&flat_vector, f));
    }
    let buf = get_reader(&filename);
    if fasta_reads {
        let reader = seq_io::fasta::Reader::new(buf);
        read_process_fasta_records(reader, threads as u32, queue_len, query_process_read_fasta_kstrobe, |record, found| {main_thread_kstrobe(found)});
        let id_hashes_view = Arc::try_unwrap(id_hashes).unwrap().into_read_only();
        let id_lengths_view = Arc::try_unwrap(id_lengths).unwrap().into_read_only();
        kstrobes::output_paf(&mut all_nams, params.l, id_hashes_view, id_lengths_view, &mut paf_file, params.k);

    }
    else {
        let reader = seq_io::fastq::Reader::new(buf);
        read_process_fastq_records(reader, threads as u32, queue_len, query_process_read_fastq_kstrobe, |record, found| {main_thread_kstrobe(found)});
        let id_hashes_view = Arc::try_unwrap(id_hashes).unwrap().into_read_only();
        let id_lengths_view = Arc::try_unwrap(id_lengths).unwrap().into_read_only();
        kstrobes::output_paf(&mut all_nams, params.l, id_hashes_view, id_lengths_view, &mut paf_file, params.k);
    }
}



pub fn run_randstrobes(filename: &PathBuf, ref_filename: &PathBuf, params: &Params, threads: usize, queue_len: usize, fasta_reads: bool, ref_fasta_reads: bool, output_prefix: &PathBuf) {
    let mut id_hashes : Arc<DashMap<u64, String>> = Arc::new(DashMap::new());
    let mut id_lengths : Arc<DashMap<u64, usize>> = Arc::new(DashMap::new());
    let mut tmp_index : Arc<PosIndex> = Arc::new(PosIndex::new());
    let mut mers_index : Arc<KmerLookup> = Arc::new(KmerLookup::new());
    let mut filter_cutoff : Arc<DashMap<usize, usize>> = Arc::new(DashMap::new());
    let mut all_nams = Vec::<(Vec<StrobeNAM>, u64)>::new();
    let mut all_mers_vector : Arc<DashMap<usize, MersVectorReduced>> = Arc::new(DashMap::new());
    let path = format!("{}{}", output_prefix.to_str().unwrap(), ".paf");
    let mut paf_file = match File::create(&path) {
        Err(why) => panic!("Couldn't create {}: {}", path, why.description()),
        Ok(paf_file) => paf_file,
    };
    let mut index_randstrobes = |seq_id_hash: u64, seq_id: &str, seq: &[u8], params: &Params, read: bool| -> MersVectorRead {
        let mut randstrobes = strobes::index(seq, seq_id_hash, params, read);
        //println!("Seq: {:?}, read_seq: {:?}", seq, read_seq);
        id_hashes.insert(seq_id_hash, seq_id.to_string());
        id_lengths.insert(seq_id_hash, seq.len());
        tmp_index.insert(seq_id_hash, randstrobes.to_vec()); 
        return randstrobes;   
    };
    let ref_process_read_aux_strobe = |ref_str: &[u8], ref_id: &str| -> Option<(MersVectorRead, u64)> {
        let thread_id :usize =  thread_id::get();
        //let seq = std::str::from_utf8(ref_str).unwrap().replace("\n", "").replace("\r", "");
        //let seq = seq_for_ref; 
        let ref_id_hash = hash_id(&ref_id);
        let mut randstrobes = index_randstrobes(ref_id_hash, ref_id, ref_str, params, false);
        return Some((randstrobes, ref_id_hash))
    
    };
    let ref_process_read_fasta_strobe = |record: seq_io::fasta::RefRecord, found: &mut Option<(MersVectorRead, u64)>| {
        let ref_str = record.seq(); 
        let ref_id = record.id().unwrap().to_string();
        *found = ref_process_read_aux_strobe(&ref_str, &ref_id);
    };
    let ref_process_read_fastq_strobe = |record: seq_io::fastq::RefRecord, found: &mut Option<(MersVectorRead, u64)>| {
        let ref_str = record.seq(); 
        let ref_id = record.id().unwrap().to_string();
        *found = ref_process_read_aux_strobe(&ref_str, &ref_id);
    };
    let mut ref_main_thread_strobe = |found: &mut Option<(MersVectorRead, u64)>| { // runs in main thread
        let (randstrobes, seq_id_hash) = found.as_ref().unwrap();
        None::<()>
    };
    let query_process_read_aux_strobe = |seq_str: &[u8], seq_id: &str| -> Option<(Vec<StrobeNAM>, u64)> {
        let thread_id :usize =  thread_id::get();
        let mut output : Vec<StrobeNAM> = Vec::new();
        //let seq = std::str::from_utf8(seq_str).unwrap(); // see https://docs.rs/seq_io/0.4.0-alpha.0/seq_io/fasta/index.html
        let seq_id_hash = hash_id(&seq_id);
        let mut randstrobes = index_randstrobes(seq_id_hash, seq_id, seq_str, params, true);
        let fc = *filter_cutoff.get(&0).unwrap();
        let all_mers = all_mers_vector.get(&1).unwrap();
        let hit_upper_window_lim = (params.l-params.s+1)*params.wmax;
        output = strobes::find_nams(&randstrobes, &all_mers, &mers_index, params.l, seq_str, hit_upper_window_lim, fc);
        if output.len() > 0 {
            return Some((output, seq_id_hash))
        }
        else {return None}
    };
    
    let query_process_read_fasta_strobe = |record: seq_io::fasta::RefRecord, found: &mut Option<(Vec<StrobeNAM>, u64)>| {
        let seq_str = record.seq(); 
        let seq_id = record.id().unwrap().to_string();
        *found = query_process_read_aux_strobe(&seq_str, &seq_id);
    
    };
    let query_process_read_fastq_strobe = |record: seq_io::fastq::RefRecord, found: &mut Option<(Vec<StrobeNAM>, u64)>| {
        let seq_str = record.seq(); 
        let seq_id = record.id().unwrap().to_string();
        *found = query_process_read_aux_strobe(&seq_str, &seq_id);
    };
    let mut main_thread_strobe = |found: &Option<(Vec<StrobeNAM>, u64)>| { // runs in main thread
        if found.is_some() {
            let (nams, seq_id_hash) = found.as_ref().unwrap();
            all_nams.push((nams.to_vec(), *seq_id_hash));
        }
        // Some(value) will stop the reader, and the value will be returned.
        // In the case of never stopping, we need to give the compiler a hint about the
        // type parameter, thus the special 'turbofish' notation is needed.
        None::<()>
    };
    let mut construct_flat_vector = || -> (MersVectorRead, u64) {
        let mut flat_vector = MersVectorRead::new();
        for entry in tmp_index.iter_mut() {
            for t in entry.iter() {flat_vector.push(*t);}
        }
        flat_vector.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        let mut prev_k = flat_vector[0].0;
        let mut unique_elements = 1;
        for entry in flat_vector.iter() {
            let mut curr_k = entry.0;
            if curr_k != prev_k {unique_elements += 1;}
            prev_k = curr_k;
        }
        return (flat_vector, unique_elements);
    };

    let mut index_vector = |flat_vector: &MersVectorRead, f: f64| -> usize {
        println!("Flat vector size: {}", flat_vector.len());
        let mut offset = 0;
        let mut prev_offset = 0;
        let mut count = 0;
        let mut tot_occur_once = 0;
        let mut tot_high_ab = 0;
        let mut tot_mid_ab = 0;
        let mut strobemer_counts = Vec::<usize>::new();
        let mut t = flat_vector[0];
        let mut prev_k = t.0;
        let mut curr_k = 0;
        for entry in flat_vector.iter() {
            curr_k = entry.0;
            if curr_k == prev_k {count += 1;}
            else {
                if count == 1 {tot_occur_once += 1;}
                else if count > 100 {tot_high_ab += 1;}
                else {tot_mid_ab += 1;}
                strobemer_counts.push(count);
                let mut s = (prev_offset, count);
                mers_index.insert(prev_k, s);
                count = 1;
                prev_k = curr_k;
                prev_offset = offset;
            }
            offset += 1;      
        }
        let mut s = (prev_offset, count);
        mers_index.insert(curr_k, s);
        if count == 1 {tot_occur_once += 1;}
        else if count > 100 {tot_high_ab += 1;}
        else {tot_mid_ab += 1;}
        strobemer_counts.push(count);
        println!("Total strobemers count:\t{}", offset);
        println!("Total strobemers unique:\t{}", tot_occur_once);
        println!("Total strobemers high-ab:\t{}", tot_high_ab);
        println!("Total strobemers mid-ab:\t{}", tot_mid_ab);
        println!("Total distinct strobemers:\t{}", mers_index.len());
        if tot_high_ab >= 1 {println!("Ratio distinct/hi-ab:\t{}", mers_index.len()/tot_high_ab);}
        if tot_mid_ab >= 1 {println!("Ratio distinct/non-distinct:\t{}", mers_index.len()/(tot_high_ab+tot_mid_ab));}
        if params.f != 0.0 {
            let mut index_cutoff = (strobemer_counts.len() as f64 * f) as usize;
            println!("Filtered cutoff index:\t{}", index_cutoff);
            let mut filter_cutoff = strobemer_counts[index_cutoff];
            println!("Filtered cutoff count:\t{}", filter_cutoff);
            return filter_cutoff;
        }
        else {return 0;}
    };
    let mut remove_kmer_hash_from_flat_vector = |flat_vector: &MersVectorRead| -> MersVectorReduced {
        let mut mers_vector_reduced = MersVectorReduced::new();
        for entry in flat_vector.iter() {
            let s = (entry.1, entry.2, entry.3);
            mers_vector_reduced.push(s);
        }
        return mers_vector_reduced;
    };
    let buf = get_reader(&ref_filename);
    if ref_fasta_reads {
        let reader = seq_io::fasta::Reader::new(buf);
        read_process_fasta_records(reader, threads as u32, queue_len, ref_process_read_fasta_strobe, |record, found| {ref_main_thread_strobe(found)});
        let (mut all_mers_vector_tmp, unique_elements) = construct_flat_vector();
        println!("Unique strobemers:\t{}", unique_elements);
        let f = params.f;
        filter_cutoff.insert(0, index_vector(&all_mers_vector_tmp, f));
        tmp_index.clear();
        all_mers_vector.insert(1, remove_kmer_hash_from_flat_vector(&all_mers_vector_tmp));

    }
    else {
        let reader = seq_io::fastq::Reader::new(buf);
        read_process_fastq_records(reader, threads as u32, queue_len, ref_process_read_fastq_strobe, |record, found| {ref_main_thread_strobe(found)});
        let (mut all_mers_vector_tmp, unique_elements) = construct_flat_vector();
        println!("Unique strobemers:\t{}", unique_elements);
        let f = params.f;
        filter_cutoff.insert(0, index_vector(&all_mers_vector_tmp, f));
        tmp_index.clear();
        all_mers_vector.insert(1, remove_kmer_hash_from_flat_vector(&all_mers_vector_tmp));
    }
    let buf = get_reader(&filename);
    if fasta_reads {
        let reader = seq_io::fasta::Reader::new(buf);
        read_process_fasta_records(reader, threads as u32, queue_len, query_process_read_fasta_strobe, |record, found| {main_thread_strobe(found)});
        let id_hashes_view = Arc::try_unwrap(id_hashes).unwrap().into_read_only();
        let id_lengths_view = Arc::try_unwrap(id_lengths).unwrap().into_read_only();
        strobes::output_paf(&mut all_nams, params.l, id_hashes_view, id_lengths_view, &mut paf_file);

    }
    else {
        let reader = seq_io::fastq::Reader::new(buf);
        read_process_fastq_records(reader, threads as u32, queue_len, query_process_read_fastq_strobe, |record, found| {main_thread_strobe(found)});
        let id_hashes_view = Arc::try_unwrap(id_hashes).unwrap().into_read_only();
        let id_lengths_view = Arc::try_unwrap(id_lengths).unwrap().into_read_only();
        strobes::output_paf(&mut all_nams, params.l, id_hashes_view, id_lengths_view, &mut paf_file);
    }
}

pub fn run_kminmers(filename: &PathBuf, ref_filename: &PathBuf, params: &Params, threads: usize, queue_len: usize, fasta_reads: bool, ref_fasta_reads: bool, output_prefix: &PathBuf) {
    let mut id_hashes : Arc<DashMap<u64, String>> = Arc::new(DashMap::new());
    let mut id_lengths : Arc<DashMap<u64, usize>> = Arc::new(DashMap::new());
    let mut tmp_index : Arc<PosIndex> = Arc::new(PosIndex::new());
    let mut mers_index : Arc<KmerLookup> = Arc::new(KmerLookup::new());
    let mut filter_cutoff : Arc<DashMap<usize, usize>> = Arc::new(DashMap::new());
    let mut all_nams = Vec::<(Vec<KminmerNAM>, u64)>::new();
    let mut all_mers_vector : Arc<DashMap<usize, MersVectorReduced>> = Arc::new(DashMap::new());
    let path = format!("{}{}", output_prefix.to_str().unwrap(), ".paf");
    let mut paf_file = match File::create(&path) {
        Err(why) => panic!("Couldn't create {}: {}", path, why.description()),
        Ok(paf_file) => paf_file,
    };
    let mut index_kminmers = |seq_id_hash: u64, seq_id: &str, seq: &[u8], params: &Params, read: bool| -> MersVectorRead {
        let mut kminmers = kminmers::index(seq, seq_id_hash, params, read);
        //println!("Seq: {:?}, read_seq: {:?}", seq, read_seq);
        id_hashes.insert(seq_id_hash, seq_id.to_string());
        id_lengths.insert(seq_id_hash, seq.len());
        tmp_index.insert(seq_id_hash, kminmers.to_vec()); 
        return kminmers;   
    };
    let ref_process_read_aux_kminmer = |ref_str: &[u8], ref_id: &str| -> Option<(MersVectorRead, u64)> {
        let thread_id :usize =  thread_id::get();
        //let seq = std::str::from_utf8(ref_str).unwrap().replace("\n", "").replace("\r", "");
        //let seq = seq_for_ref; 
        let ref_id_hash = hash_id(&ref_id);
        let mut kminmers = index_kminmers(ref_id_hash, ref_id, ref_str, params, false);
        return Some((kminmers, ref_id_hash))
    
    };
    let ref_process_read_fasta_kminmer = |record: seq_io::fasta::RefRecord, found: &mut Option<(MersVectorRead, u64)>| {
        let ref_str = record.seq(); 
        let ref_id = record.id().unwrap().to_string();
        *found = ref_process_read_aux_kminmer(&ref_str, &ref_id);
    };
    let ref_process_read_fastq_kminmer = |record: seq_io::fastq::RefRecord, found: &mut Option<(MersVectorRead, u64)>| {
        let ref_str = record.seq(); 
        let ref_id = record.id().unwrap().to_string();
        *found = ref_process_read_aux_kminmer(&ref_str, &ref_id);
    };
    let mut ref_main_thread_kminmer = |found: &mut Option<(MersVectorRead, u64)>| { // runs in main thread
        let (kminmers, seq_id_hash) = found.as_ref().unwrap();
        None::<()>
    };
    let query_process_read_aux_kminmer = |seq_str: &[u8], seq_id: &str| -> Option<(Vec<KminmerNAM>, u64)> {
        let thread_id :usize =  thread_id::get();
        let mut output : Vec<KminmerNAM> = Vec::new();
        //let seq = std::str::from_utf8(seq_str).unwrap(); // see https://docs.rs/seq_io/0.4.0-alpha.0/seq_io/fasta/index.html
        let seq_id_hash = hash_id(&seq_id);
        let mut kminmers = index_kminmers(seq_id_hash, seq_id, seq_str, params, true);
        let fc = *filter_cutoff.get(&0).unwrap();
        let all_mers = all_mers_vector.get(&1).unwrap();
        output = kminmers::find_nams(&kminmers, &all_mers, &mers_index, params.l, seq_str, fc);
        if output.len() > 0 {
            return Some((output, seq_id_hash))
        }
        else {return None}
    };
    
    let query_process_read_fasta_kminmer = |record: seq_io::fasta::RefRecord, found: &mut Option<(Vec<KminmerNAM>, u64)>| {
        let seq_str = record.seq(); 
        let seq_id = record.id().unwrap().to_string();
        *found = query_process_read_aux_kminmer(&seq_str, &seq_id);
    
    };
    let query_process_read_fastq_kminmer = |record: seq_io::fastq::RefRecord, found: &mut Option<(Vec<KminmerNAM>, u64)>| {
        let seq_str = record.seq(); 
        let seq_id = record.id().unwrap().to_string();
        *found = query_process_read_aux_kminmer(&seq_str, &seq_id);
    };
    let mut main_thread_kminmer = |found: &Option<(Vec<KminmerNAM>, u64)>| { // runs in main thread
        if found.is_some() {
            let (nams, seq_id_hash) = found.as_ref().unwrap();
            all_nams.push((nams.to_vec(), *seq_id_hash));
        }
        // Some(value) will stop the reader, and the value will be returned.
        // In the case of never stopping, we need to give the compiler a hint about the
        // type parameter, thus the special 'turbofish' notation is needed.
        None::<()>
    };
    let mut construct_flat_vector = || -> (MersVectorRead, u64) {
        let mut flat_vector = MersVectorRead::new();
        for entry in tmp_index.iter_mut() {
            for t in entry.iter() {flat_vector.push(*t);}
        }
        flat_vector.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        let mut prev_k = flat_vector[0].0;
        let mut unique_elements = 1;
        for entry in flat_vector.iter() {
            let mut curr_k = entry.0;
            if curr_k != prev_k {unique_elements += 1;}
            prev_k = curr_k;
        }
        return (flat_vector, unique_elements);
    };

    let mut index_vector = |flat_vector: &MersVectorRead, f: f64| -> usize {
        println!("Flat vector size: {}", flat_vector.len());
        let mut offset = 0;
        let mut prev_offset = 0;
        let mut count = 0;
        let mut tot_occur_once = 0;
        let mut tot_high_ab = 0;
        let mut tot_mid_ab = 0;
        let mut kminmer_counts = Vec::<usize>::new();
        let mut t = flat_vector[0];
        let mut prev_k = t.0;
        let mut curr_k = 0;
        for entry in flat_vector.iter() {
            curr_k = entry.0;
            if curr_k == prev_k {count += 1;}
            else {
                if count == 1 {tot_occur_once += 1;}
                else if count > 100 {tot_high_ab += 1;}
                else {tot_mid_ab += 1;}
                kminmer_counts.push(count);
                let mut s = (prev_offset, count);
                mers_index.insert(prev_k, s);
                count = 1;
                prev_k = curr_k;
                prev_offset = offset;
            }
            offset += 1;      
        }
        let mut s = (prev_offset, count);
        mers_index.insert(curr_k, s);
        mers_index.insert(prev_k, s);
        if count == 1 {tot_occur_once += 1;}
        else if count > 100 {tot_high_ab += 1;}
        else {tot_mid_ab += 1;}
        kminmer_counts.push(count);
        println!("Total k-min-mers count:\t{}", offset);
        println!("Total k-min-mers unique:\t{}", tot_occur_once);
        println!("Total k-min-mers high-ab:\t{}", tot_high_ab);
        println!("Total k-min-mers mid-ab:\t{}", tot_mid_ab);
        println!("Total distinct k-min-mers:\t{}", mers_index.len());
        if tot_high_ab >= 1 {println!("Ratio distinct/hi-ab:\t{}", mers_index.len()/tot_high_ab);}
        if tot_mid_ab >= 1 {println!("Ratio distinct/non-distinct:\t{}", mers_index.len()/(tot_high_ab+tot_mid_ab));}
        if params.f != 0.0 {
            let mut index_cutoff = (kminmer_counts.len() as f64 * f) as usize;
            println!("Filtered cutoff index:\t{}", index_cutoff);
            let mut filter_cutoff = kminmer_counts[index_cutoff];
            println!("Filtered cutoff count:\t{}", filter_cutoff);
            return filter_cutoff;
        }
        else {return 0;}
    };
    let mut remove_kmer_hash_from_flat_vector = |flat_vector: &MersVectorRead| -> MersVectorReduced {
        let mut mers_vector_reduced = MersVectorReduced::new();
        for entry in flat_vector.iter() {
            let s = (entry.1, entry.2, entry.3);
            mers_vector_reduced.push(s);
        }
        return mers_vector_reduced;
    };
    let buf = get_reader(&ref_filename);
    if ref_fasta_reads {
        let reader = seq_io::fasta::Reader::new(buf);
        read_process_fasta_records(reader, threads as u32, queue_len, ref_process_read_fasta_kminmer, |record, found| {ref_main_thread_kminmer(found)});
        let (mut all_mers_vector_tmp, unique_elements) = construct_flat_vector();
        println!("Unique k-min-mers:\t{}", unique_elements);
        let f = params.f;
        filter_cutoff.insert(0, index_vector(&all_mers_vector_tmp, f));
        tmp_index.clear();
        all_mers_vector.insert(1, remove_kmer_hash_from_flat_vector(&all_mers_vector_tmp));

    }
    else {
        let reader = seq_io::fastq::Reader::new(buf);
        read_process_fastq_records(reader, threads as u32, queue_len, ref_process_read_fastq_kminmer, |record, found| {ref_main_thread_kminmer(found)});
        let (mut all_mers_vector_tmp, unique_elements) = construct_flat_vector();
        println!("Unique k-min-mers:\t{}", unique_elements);
        let f = params.f;
        filter_cutoff.insert(0, index_vector(&all_mers_vector_tmp, f));
        tmp_index.clear();
        all_mers_vector.insert(1, remove_kmer_hash_from_flat_vector(&all_mers_vector_tmp));
    }
    let buf = get_reader(&filename);
    if fasta_reads {
        let reader = seq_io::fasta::Reader::new(buf);
        read_process_fasta_records(reader, threads as u32, queue_len, query_process_read_fasta_kminmer, |record, found| {main_thread_kminmer(found)});
        let id_hashes_view = Arc::try_unwrap(id_hashes).unwrap().into_read_only();
        let id_lengths_view = Arc::try_unwrap(id_lengths).unwrap().into_read_only();
        kminmers::output_paf(&mut all_nams, params.l, id_hashes_view, id_lengths_view, &mut paf_file);

    }
    else {
        let reader = seq_io::fastq::Reader::new(buf);
        read_process_fastq_records(reader, threads as u32, queue_len, query_process_read_fastq_kminmer, |record, found| {main_thread_kminmer(found)});
        let id_hashes_view = Arc::try_unwrap(id_hashes).unwrap().into_read_only();
        let id_lengths_view = Arc::try_unwrap(id_lengths).unwrap().into_read_only();
        kminmers::output_paf(&mut all_nams, params.l, id_hashes_view, id_lengths_view, &mut paf_file);
    }
}

pub fn run_ksyncmers(filename: &PathBuf, ref_filename: &PathBuf, params: &Params, threads: usize, queue_len: usize, fasta_reads: bool, ref_fasta_reads: bool, output_prefix: &PathBuf) {
    let mut id_hashes : Arc<DashMap<u64, String>> = Arc::new(DashMap::new());
    let mut id_lengths : Arc<DashMap<u64, usize>> = Arc::new(DashMap::new());
    let mut hash_counts : Arc<DashMap<u64, usize>> = Arc::new(DashMap::new());
    let mut tmp_index : Arc<PosIndex> = Arc::new(PosIndex::new());
    let mut id_syncmers_pos : Arc<DashMap<u64, (Vec<u64>, Vec<usize>)>> = Arc::new(DashMap::new());
    let mut mers_index : Arc<KmerLookupMod> = Arc::new(KmerLookupMod::new());
    let mut filter_cutoff : Arc<DashMap<usize, usize>> = Arc::new(DashMap::new());
    let mut all_nams = Vec::<(Vec<KsyncmerNAM>, u64)>::new();
    let mut all_mers_vector : Arc<DashMap<usize, MersVectorRead>> = Arc::new(DashMap::new());
    let path = format!("{}{}", output_prefix.to_str().unwrap(), ".paf");
    let mut paf_file = match File::create(&path) {
        Err(why) => panic!("Couldn't create {}: {}", path, why.description()),
        Ok(paf_file) => paf_file,
    };
    let mut index_ksyncmers = |seq_id_hash: u64, seq_id: &str, seq: &[u8], params: &Params, read: bool| -> MersVectorRead {
        let (mut string_hashes, mut pos_to_seq_coord, mut ksyncmers) = ksyncmers::index(seq, seq_id_hash, params, read);
        //println!("Seq: {:?}, read_seq: {:?}", seq, read_seq);
        id_hashes.insert(seq_id_hash, seq_id.to_string());
        id_lengths.insert(seq_id_hash, seq.len());
        tmp_index.insert(seq_id_hash, ksyncmers.to_vec());

        return ksyncmers;   
    };
    let ref_process_read_aux_ksyncmer = |ref_str: &[u8], ref_id: &str| -> Option<(MersVectorRead, u64)> {
        let thread_id :usize =  thread_id::get();
        //let seq = std::str::from_utf8(ref_str).unwrap().replace("\n", "").replace("\r", "");
        //let seq = seq_for_ref; 
        let ref_id_hash = hash_id(&ref_id);
        let mut ksyncmers = index_ksyncmers(ref_id_hash, ref_id, ref_str, params, false);
        return Some((ksyncmers, ref_id_hash))
    
    };
    let ref_process_read_fasta_ksyncmer = |record: seq_io::fasta::RefRecord, found: &mut Option<(MersVectorRead, u64)>| {
        let ref_str = record.seq(); 
        let ref_id = record.id().unwrap().to_string();
        *found = ref_process_read_aux_ksyncmer(&ref_str, &ref_id);
    };
    let ref_process_read_fastq_ksyncmer = |record: seq_io::fastq::RefRecord, found: &mut Option<(MersVectorRead, u64)>| {
        let ref_str = record.seq(); 
        let ref_id = record.id().unwrap().to_string();
        *found = ref_process_read_aux_ksyncmer(&ref_str, &ref_id);
    };
    let mut ref_main_thread_ksyncmer = |found: &mut Option<(MersVectorRead, u64)>| { // runs in main thread
        let (ksyncmers, seq_id_hash) = found.as_ref().unwrap();
        None::<()>
    };
    let query_process_read_aux_ksyncmer = |seq_str: &[u8], seq_id: &str| -> Option<(Vec<KsyncmerNAM>, u64)> {
        let thread_id :usize =  thread_id::get();
        let mut output : Vec<KsyncmerNAM> = Vec::new();
        //let seq = std::str::from_utf8(seq_str).unwrap(); // see https://docs.rs/seq_io/0.4.0-alpha.0/seq_io/fasta/index.html
        let seq_id_hash = hash_id(&seq_id);
        let mut ksyncmers = index_ksyncmers(seq_id_hash, seq_id, seq_str, params, true);
        let fc = *filter_cutoff.get(&0).unwrap();
        output = ksyncmers::find_nams(&ksyncmers, &tmp_index, &mers_index, &hash_counts, params.l, seq_str, fc, params.k, &id_hashes, &id_lengths);
        if output.len() > 0 {
            return Some((output, seq_id_hash))
        }
        else {return None}
    };
    
    let query_process_read_fasta_ksyncmer = |record: seq_io::fasta::RefRecord, found: &mut Option<(Vec<KsyncmerNAM>, u64)>| {
        let seq_str = record.seq(); 
        let seq_id = record.id().unwrap().to_string();
        *found = query_process_read_aux_ksyncmer(&seq_str, &seq_id);
    
    };
    let query_process_read_fastq_ksyncmer = |record: seq_io::fastq::RefRecord, found: &mut Option<(Vec<KsyncmerNAM>, u64)>| {
        let seq_str = record.seq(); 
        let seq_id = record.id().unwrap().to_string();
        *found = query_process_read_aux_ksyncmer(&seq_str, &seq_id);
    };
    let mut main_thread_ksyncmer = |found: &Option<(Vec<KsyncmerNAM>, u64)>| { // runs in main thread
        if found.is_some() {
            let (nams, seq_id_hash) = found.as_ref().unwrap();
            all_nams.push((nams.to_vec(), *seq_id_hash));
        }
        // Some(value) will stop the reader, and the value will be returned.
        // In the case of never stopping, we need to give the compiler a hint about the
        // type parameter, thus the special 'turbofish' notation is needed.
        None::<()>
    };
    let mut construct_flat_vector = || -> (MersVectorRead, u64) {
        let mut flat_vector = MersVectorRead::new();
        for entry in tmp_index.iter_mut() {
            for t in entry.iter() {flat_vector.push(*t);}
        }
        flat_vector.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        let mut prev_k = flat_vector[0].0;
        let mut unique_elements = 1;
        for entry in flat_vector.iter() {
            let mut curr_k = entry.0;
            if curr_k != prev_k {unique_elements += 1;}
            prev_k = curr_k;
        }
        return (flat_vector, unique_elements);
    };

    let mut index_vector = |tmp_index: &PosIndex| {
        for entry in tmp_index.into_iter() {
            let ref_id = entry.key();
            let flat_vector = entry.value();
            println!("Flat vector size: {}", flat_vector.len());
            let mut offset = 0;
            for entry in flat_vector.iter() {
                let curr_k = entry.0;
                let mut s = (*ref_id, offset);
                mers_index.insert(curr_k, s);
                offset += 1;      
            }
        }
    };

    let mut count_vector = |flat_vector: &MersVectorRead, f: f64| -> usize {
        println!("Flat vector size: {}", flat_vector.len());
        let mut offset = 0;
        let mut prev_offset = 0;
        let mut count = 0;
        let mut tot_occur_once = 0;
        let mut tot_high_ab = 0;
        let mut tot_mid_ab = 0;
        let mut kminmer_counts = Vec::<usize>::new();
        let mut t = flat_vector[0];
        let mut prev_k = t.0;
        let mut curr_k = 0;
        for entry in flat_vector.iter() {
            curr_k = entry.0;
            if curr_k == prev_k {count += 1;}
            else {
                if count == 1 {tot_occur_once += 1;}
                else if count > 100 {tot_high_ab += 1;}
                else {tot_mid_ab += 1;}
                kminmer_counts.push(count);
                hash_counts.insert(prev_k, count);
                count = 1;
                prev_k = curr_k;
                prev_offset = offset;
            }
            offset += 1;      
        }
        if count == 1 {tot_occur_once += 1;}
        else if count > 100 {tot_high_ab += 1;}
        else {tot_mid_ab += 1;}
        kminmer_counts.push(count);
        hash_counts.insert(curr_k, count);
        println!("Total k-syncmers count:\t{}", offset);
        println!("Total k-syncmers unique:\t{}", tot_occur_once);
        println!("Total k-syncmers high-ab:\t{}", tot_high_ab);
        println!("Total k-syncmers mid-ab:\t{}", tot_mid_ab);
        if params.f != 0.0 {
            kminmer_counts.sort_unstable();
            kminmer_counts.reverse();
            let mut index_cutoff = (kminmer_counts.len() as f64 * f) as usize;
            println!("Filtered cutoff index:\t{}", index_cutoff);
            let mut filter_cutoff = kminmer_counts[index_cutoff];
            println!("Filtered cutoff count:\t{}", filter_cutoff);
            return filter_cutoff;
        }
        else {return 0;}
    };
    let mut remove_kmer_hash_from_flat_vector = |flat_vector: &MersVectorRead| -> MersVectorReduced {
        let mut mers_vector_reduced = MersVectorReduced::new();
        for entry in flat_vector.iter() {
            let s = (entry.1, entry.2, entry.3);
            mers_vector_reduced.push(s);
        }
        return mers_vector_reduced;
    };
    let buf = get_reader(&ref_filename);
    if ref_fasta_reads {
        let reader = seq_io::fasta::Reader::new(buf);
        read_process_fasta_records(reader, threads as u32, queue_len, ref_process_read_fasta_ksyncmer, |record, found| {ref_main_thread_ksyncmer(found)});
        let (mut flat_vector, unique_elements) = construct_flat_vector();
        println!("Unique k-syncmers:\t{}", unique_elements);
        let f = params.f;
        index_vector(&tmp_index);
        filter_cutoff.insert(0, count_vector(&flat_vector, f));

    }
    else {
        let reader = seq_io::fastq::Reader::new(buf);
        read_process_fastq_records(reader, threads as u32, queue_len, ref_process_read_fastq_ksyncmer, |record, found| {ref_main_thread_ksyncmer(found)});
        let (mut flat_vector, unique_elements) = construct_flat_vector();
        println!("Unique k-syncmers:\t{}", unique_elements);
        let f = params.f;
        index_vector(&tmp_index);
        filter_cutoff.insert(0, count_vector(&flat_vector, f));
    }
    let buf = get_reader(&filename);
    if fasta_reads {
        let reader = seq_io::fasta::Reader::new(buf);
        read_process_fasta_records(reader, threads as u32, queue_len, query_process_read_fasta_ksyncmer, |record, found| {main_thread_ksyncmer(found)});
        let id_hashes_view = Arc::try_unwrap(id_hashes).unwrap().into_read_only();
        let id_lengths_view = Arc::try_unwrap(id_lengths).unwrap().into_read_only();
        ksyncmers::output_paf(&mut all_nams, params.l, id_hashes_view, id_lengths_view, &mut paf_file, params.k);

    }
    else {
        let reader = seq_io::fastq::Reader::new(buf);
        read_process_fastq_records(reader, threads as u32, queue_len, query_process_read_fastq_ksyncmer, |record, found| {main_thread_ksyncmer(found)});
        let id_hashes_view = Arc::try_unwrap(id_hashes).unwrap().into_read_only();
        let id_lengths_view = Arc::try_unwrap(id_lengths).unwrap().into_read_only();
        ksyncmers::output_paf(&mut all_nams, params.l, id_hashes_view, id_lengths_view, &mut paf_file, params.k);
    }
}






















/*let ref_process_read_fasta = |record: seq_io::fasta::RefRecord, found: &mut Option<(Vec<(Kmer, String, bool, String, (usize, usize))>, Option<Read>, Option<ReadSync>)>| {
    let ref_str = record.seq(); 
    let ref_id = record.id().unwrap().to_string();
    *found = ref_process_read_aux(&ref_str, &ref_id);
};
let ref_process_read_fastq = |record: seq_io::fastq::RefRecord, found: &mut Option<(Vec<(Kmer, String, bool, String, (usize, usize))>,Option<Read>, Option<ReadSync>)>| {
    let ref_str = record.seq(); 
    let ref_id = record.id().unwrap().to_string();
    *found = ref_process_read_aux(&ref_str, &ref_id);
};
let mut add_ref_kminmer = |ref_kminmer: Kmer, seq: Option<&str>, seq_reversed: &bool, origin: &str, shift: &(usize, usize), sequences_file: &mut SeqFileType, thread_id: usize, read_seq: Option<&str>, read_offsets: Option<(usize, usize, usize)>| {
    //println!("Seq: {:?}, read_seq: {:?}", seq, read_seq);
    let seq = if *seq_reversed {utils::revcomp(&seq.unwrap())} else {seq.unwrap().to_string()};
    let seqlen = read_seq.unwrap().len() as u32;
    let seq_line = format!("{}\t{}\t{}\t{}\t{:?}", ref_kminmer.print_as_string(), seq, "*", origin, shift);
    seq_write(sequences_file, format!("{}\n", seq_line));
    let mut contains_key = kmer_table.contains_key(&ref_kminmer);
    if contains_key {
        let mut entry_vec = kmer_table.get_mut(&ref_kminmer).unwrap();
        entry_vec.push(HashEntry{origin: origin.to_string(), seq: seq, seqlen: seqlen, shift: *shift, seq_rev: *seq_reversed}); 
    }
    else {
        let mut entry_vec = Vec::<HashEntry>::new();
        entry_vec.push(HashEntry{origin: origin.to_string(), seq: seq, seqlen: seqlen, shift: *shift, seq_rev: *seq_reversed});
        kmer_table.insert(ref_kminmer.clone(), entry_vec);
    }
    let mut contains_abund = ref_abundance_table.contains_key(&ref_kminmer);
    if contains_abund {*ref_abundance_table.get_mut(&ref_kminmer).unwrap() += 1;}
    else {ref_abundance_table.insert(ref_kminmer.clone(), 1);}
};
// worker thread
let ref_process_read_aux = |ref_str: &[u8], ref_id: &str| -> Option<(Vec<(Kmer, String, bool, String, (usize, usize))>, Option<Read>, Option<ReadSync>)> {
    let thread_id :usize =  thread_id::get();
    let mut output : Vec<(Kmer, String, bool, String, (usize, usize))> = Vec::new();
    let seq = std::str::from_utf8(ref_str).unwrap(); // see https://docs.rs/seq_io/0.4.0-alpha.0/seq_io/fasta/index.html
        // those two next lines do a string copy of the read sequence, because in the case of a
        // reference, we'll need to remove newlines. also, that sequence will be moved to the
        // Read object. could be optimized later
        let seq_for_ref = seq.replace("\n", "").replace("\r", "");
        let seq = seq_for_ref; 
    let mut ref_obj = Read::extract(&ref_id, seq, &params, &minimizer_to_int, &int_to_minimizer, &uhs_bloom, &lcp_bloom);
    if sequences_files.get(&thread_id).is_none() {
        sequences_files.insert(thread_id, create_sequences_file(thread_id));
    }
    let mut sequences_file = sequences_files.get_mut(&thread_id).unwrap();
    if ref_obj.transformed.len() > k {
        for i in 0..(ref_obj.transformed.len() - k + 1) {
            let mut ref_kminmer : Kmer = Kmer::make_from(&ref_obj.transformed[i..i+k]);
            // uncomment the line below to track where the kmer is coming from (but not needed in production)
            //let origin = "*".to_string(); 
            let mut seq_reversed = false;
            if revcomp_aware { 
                let (ref_kminmer_norm, reversed) = ref_kminmer.normalize(); 
                ref_kminmer = ref_kminmer_norm;
                seq_reversed = reversed;
            } 
            
            let minimizers_pos = &ref_obj.minimizers_pos;
            let position_of_first_minimizer = match seq_reversed {
                true => minimizers_pos[i+k-1] + params.l - 1,
                false => minimizers_pos[i]
            };
            let position_of_last_minimizer = match seq_reversed {
                true => minimizers_pos[i] + params.l - 1,
                false => minimizers_pos[i+k-1]
            };
            let mut shift = (minimizers_pos[i], minimizers_pos[i+k-1] + params.l - 1);
            let ref_offsets = (ref_obj.minimizers_pos[i] as usize, (ref_obj.minimizers_pos[i+k-1] as usize + l - 1), (ref_obj.minimizers_pos[i+k-1] + l - ref_obj.minimizers_pos[i]));
            //println!("Shift: {:?}", shift);
            let mut span_seq = &ref_obj.seq[shift.0..shift.1 + 1];
        // println!("Span_seq: {:?}, len: {:?}", span_seq, span_seq.len());
            //start of kminmer : pos[i]
            //end of kminmer : pos [i+k-1]+l
            add_ref_kminmer(ref_kminmer, Some(span_seq), &seq_reversed, &ref_id, &shift, &mut sequences_file, thread_id, Some(&ref_obj.seq), Some(ref_offsets));
        }
    }
    Some((output, Some(ref_obj), None))
};
let mut query_kminmer_lookup =|inp_kminmer: Kmer, seq: Option<&str>, seq_reversed: &bool, origin: &str, shift: &(usize, usize), sequences_file: &mut SeqFileType, thread_id: usize, read_seq: Option<&str>, read_offsets: Option<(usize, usize, usize)>| -> Option<QueryMatch>
{
    let mut contains_abund = query_abundance_table.contains_key(&inp_kminmer);
    if contains_abund {*query_abundance_table.get_mut(&inp_kminmer).unwrap() += 1;}
    else {query_abundance_table.insert(inp_kminmer.clone(), 1);}
    let seq = if *seq_reversed {utils::revcomp(&seq.unwrap())} else {seq.unwrap().to_string()};
    let seqlen = read_seq.unwrap().len() as u32;
    let mut contains_key = kmer_table.contains_key(&inp_kminmer);
    if contains_key {
        let entry_vec = kmer_table.get_mut(&inp_kminmer).unwrap();
        let seq_line = format!("{}\t{}\t{}\t{}\t{:?}", inp_kminmer.print_as_string(), seq, "*", origin, shift);
        seq_write(sequences_file, format!("{}\n", seq_line));
        return Some(QueryMatch{kminmer: inp_kminmer, target: entry_vec.to_vec(), query: HashEntry{origin: origin.to_string(), seq: seq, seqlen: seqlen, shift: *shift, seq_rev: *seq_reversed}});
    }
    else {return None};
};

let query_process_read_aux = |seq_str: &[u8], seq_id: &str| -> Option<(Vec<QueryMatch>, Option<Read>, Option<ReadSync>)> {
    let thread_id :usize =  thread_id::get();
    let mut output : Vec<(Kmer, String, bool, String, (usize, usize))> = Vec::new();
    let seq = std::str::from_utf8(seq_str).unwrap(); // see https://docs.rs/seq_io/0.4.0-alpha.0/seq_io/fasta/index.html
    // those two next lines do a string copy of the read sequence, because in the case of a
    // reference, we'll need to remove newlines. also, that sequence will be moved to the
    // Read object. could be optimized later
    let mut read_obj = Read::extract(&seq_id, seq.to_string(), &params, &minimizer_to_int, &int_to_minimizer, &uhs_bloom, &lcp_bloom);
    let mut lookups = Vec::<QueryMatch>::new();
    if sequences_files.get(&thread_id).is_none() {
        sequences_files.insert(thread_id, create_sequences_file(thread_id));
    }
    let mut sequences_file = sequences_files.get_mut(&thread_id).unwrap();
    if read_obj.transformed.len() > k {
        for i in 0..(read_obj.transformed.len() - k + 1) {
            let mut kminmer : Kmer = Kmer::make_from(&read_obj.transformed[i..i+k]);
            let mut seq_reversed = false;
            if revcomp_aware { 
                let (kminmer_norm, reversed) = kminmer.normalize(); 
                kminmer = kminmer_norm;
                seq_reversed = reversed;
            } 
            let minimizers_pos = &read_obj.minimizers_pos;
            let position_of_first_minimizer = match seq_reversed {
                true => minimizers_pos[i+k-1] + params.l - 1,
                false => minimizers_pos[i]
            };
            let position_of_last_minimizer = match seq_reversed {
                true => minimizers_pos[i] + params.l - 1,
                false => minimizers_pos[i+k-1]
            };
            let mut shift = (minimizers_pos[i], minimizers_pos[i+k-1] + params.l - 1);
            let read_offsets = (read_obj.minimizers_pos[i] as usize, (read_obj.minimizers_pos[i+k-1] as usize + l - 1), (read_obj.minimizers_pos[i+k-1] + l - read_obj.minimizers_pos[i]));
            //println!("Shift: {:?}", shift);
            let mut span_seq = &read_obj.seq[shift.0..shift.1 + 1];
            //println!("Span_seq: {:?}, len: {:?}", span_seq, span_seq.len());
            //start of kminmer : pos[i]
            //end of kminmer : pos [i+k-1]+l
            let lookup = query_kminmer_lookup(kminmer, Some(span_seq), &seq_reversed, &seq_id, &shift, &mut sequences_file, thread_id, Some(&read_obj.seq), None);
            match lookup {
                Some(x) => lookups.push(x),
                None => continue
            }
        }
    }
    if lookups.len() > 0 {
        return Some((lookups, Some(read_obj), None))
    }
    else {return None}
};

// parallel fasta parsing, with a main thread that writes to disk and populates hash tables
let mut ref_main_thread = |found: &Option<(Vec<(Kmer, String, bool, String, (usize, usize))>, Option<Read>, Option<ReadSync>)>| { // runs in main thread
    //println!("Received read in main thread, nb kmers: {}", vec.len());
    let debug_only_display_read_and_minimizers = false;
    if debug_only_display_read_and_minimizers {
        // debug: just displays the read id and the list of minimizers
        let (vec, ref_obj, ref_obj_sync) = found.as_ref().unwrap();
        if ref_obj.is_some(){println!("{} {}", ref_obj.as_ref().unwrap().id.to_string(), ref_obj.as_ref().unwrap().transformed.to_vec().iter().join(" "));}
        else {println!("{} {}", ref_obj_sync.as_ref().unwrap().id.to_string(), ref_obj_sync.as_ref().unwrap().transformed.to_vec().iter().join(" "));}
    }
    else {
        let (vec, ref_obj, ref_obj_sync) = found.as_ref().unwrap();
    }
    // Some(value) will stop the reader, and the value will be returned.
    // In the case of never stopping, we need to give the compiler a hint about the
    // type parameter, thus the special 'turbofish' notation is needed.
    None::<()>
};
let query_process_read_fasta = |record: seq_io::fasta::RefRecord, found: &mut Option<(Vec<QueryMatch>, Option<Read>, Option<ReadSync>)>| {
    let seq_str = record.seq(); 
    let seq_id = record.id().unwrap().to_string();
    *found = query_process_read_aux(&seq_str, &seq_id);
};
let query_process_read_fastq = |record: seq_io::fastq::RefRecord, found: &mut Option<(Vec<QueryMatch>, Option<Read>, Option<ReadSync>)>| {
    let seq_str = record.seq(); 
    let seq_id = record.id().unwrap().to_string();
    *found = query_process_read_aux(&seq_str, &seq_id);
};


// parallel fasta parsing, with a main thread that writes to disk and populates hash tables
let mut main_thread = |found: &Option<(Vec<QueryMatch>, Option<Read>, Option<ReadSync>)>| { // runs in main thread
    if found.is_some() {
        //println!("Received read in main thread, nb kmers: {}", vec.len());
        let debug_only_display_matches = false;
        let (query_lookups, read_obj, read_sync_obj) = found.as_ref().unwrap();
        if debug_only_display_matches {
            for qmatch in query_lookups.iter() {
                println!("K-min-mer {:?}, target(s) {:?}, query {:?}", qmatch.kminmer.print_as_string(), qmatch.target, qmatch.query);
            } // kminmer: Kmer, target: Vec<HashEntry>, query: HashEntry, unique: bool};
        }
        if read_obj.is_some() {query_matches.insert(read_obj.as_ref().unwrap().id.to_string(), query_lookups.to_vec());}
        else if read_sync_obj.is_some(){query_matches.insert(read_sync_obj.as_ref().unwrap().id.to_string(), query_lookups.to_vec());}
    }
    // Some(value) will stop the reader, and the value will be returned.
    // In the case of never stopping, we need to give the compiler a hint about the
    // type parameter, thus the special 'turbofish' notation is needed.
    None::<()>
};*/