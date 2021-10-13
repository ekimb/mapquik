use std::io::{self};
use std::error::Error;
use std::io::{BufRead, BufReader};
use std::path::Path;
use crate::BufReadDecompressor;
use histo::Histogram;
use std::fs::{File};
use std::sync::{Arc};
use seq_io::BaseRecord;
use seq_io::parallel::{read_process_fasta_records, read_process_fastq_records};
use dashmap::DashMap;
use thread_id;
use super::mers;
use crate::mers::{Match};
use std::path::PathBuf;
use super::Params;
use crate::{get_reader, hash_id};
use indicatif::ProgressBar;
use std::collections::HashSet;
use std::hash::Hash;
use std::sync::atomic::AtomicUsize;
use crate::graph::{DbgEntry};
use super::graph;
use crate::graph::DbgIndex;
use dashmap::DashSet;
use crate::graph::Kmer;
use crate::mers::seq_to_kmers;

fn has_unique_elements<T>(iter: T) -> bool
where
    T: IntoIterator,
    T::Item: Eq + Hash,
{
    let mut uniq = HashSet::new();
    iter.into_iter().all(move |x| uniq.insert(x))
}

pub fn run_mers(filename: &PathBuf, ref_filename: &PathBuf, params: &Params, threads: usize, queue_len: usize, fasta_reads: bool, ref_fasta_reads: bool, output_prefix: &PathBuf) {
    let NODE_INDEX : AtomicUsize = AtomicUsize::new(0);
    let dbg_nodes : Arc<DashMap<Vec<u64>, DbgEntry>> =  Arc::new(DashMap::new()); 
    let dbg_edges : Arc<DashMap<DbgIndex, Vec<Kmer>>> =  Arc::new(DashMap::new());
    let f = params.f;
    let mut fmax = params.fmax;
    let mut fmin = params.fmin;
    let query_mers_index : Arc<DashMap<u64, usize>> =  Arc::new(DashMap::new());
    let mut all_matches = Vec::<(Match, String)>::new();
    let path = format!("{}{}", output_prefix.to_str().unwrap(), ".paf");
    let mut paf_file = match File::create(&path) {
        Err(why) => panic!("Couldn't create {}: {}", path, why.description()),
        Ok(paf_file) => paf_file,
    };
    let lens : Arc<DashMap<String, usize>> = Arc::new(DashMap::new());
    let count_lmers = |seq_id: &str, seq: &[u8], params: &Params| {
        let (mut string_hashes, mut pos_to_seq_coord) = mers::extract_mers(seq, params);
        for i in 0..string_hashes.len() {
            let mut entry = query_mers_index.get_mut(&string_hashes[i]);
            if entry.is_some() {*entry.unwrap() += 1;}
            else {query_mers_index.insert(string_hashes[i], 1);}
        }
    };
    let query_count_read_aux_mer = |seq_str: &[u8], seq_id: &str| -> Option<usize> {
        let mut output = Vec::<Match>::new();
        count_lmers(seq_id, seq_str, params);
        return Some(1)
    }; 
    let query_count_read_fasta_mer = |record: seq_io::fasta::RefRecord, found: &mut Option<usize>| {
        let seq_str = record.seq(); 
        let seq_id = record.id().unwrap().to_string();
        *found = query_count_read_aux_mer(&seq_str, &seq_id);
    
    };
    let query_count_read_fastq_mer = |record: seq_io::fastq::RefRecord, found: &mut Option<usize>| {
        let seq_str = record.seq(); 
        let seq_id = record.id().unwrap().to_string();
        *found = query_count_read_aux_mer(&seq_str, &seq_id);
    
    };
    let mut count_thread_mer = |found: &mut Option<usize>| { // runs in main thread
        None::<()>
    };
    let buf = get_reader(&filename);
    if fasta_reads {
        let reader = seq_io::fasta::Reader::new(buf); 
        read_process_fasta_records(reader, threads as u32, queue_len, query_count_read_fasta_mer, |record, found| {count_thread_mer(found)});
    }
    else {
        let reader = seq_io::fastq::Reader::new(buf); 
        read_process_fastq_records(reader, threads as u32, queue_len, query_count_read_fastq_mer, |record, found| {count_thread_mer(found)});
    }
    if params.f == 0.0 {
        if fmax == usize::max_value() {println!("Retaining l-mers with count >= {}.", fmin);}
        else {println!("Retaining l-mers with count >= {} and <= {}.", fmin, fmax);}
    }
    else {
        let mut frac_index = (f * query_mers_index.len() as f64) as usize;
        let mut counts : Vec<usize> = query_mers_index.iter().map(|(ent)| *ent.value()).collect();
        counts.sort_by(|a, b| a.cmp(&b));
        let mut histogram = Histogram::with_buckets(200);
        for count in counts.iter() {histogram.add(*count as u64);}
        println!("{}", histogram);
        counts.reverse();
        fmax = counts.iter().nth(frac_index).unwrap().clone();
        println!("Removed top {}% fraction of {} l-mers (count >= {}).", params.f*100.0, query_mers_index.len(), fmax);
    }

    let index_mers = |seq_id: &str, seq: &[u8], params: &Params, read: bool| -> (usize, Vec<Kmer>, Vec<Kmer>) {
        let (seq_len, mut mers, mut mers_rev) = mers::seq_to_kmers(seq, seq_id, params, read, &query_mers_index, (fmin, fmax));
        if read == true {lens.insert(seq_id.to_string(), seq_len);}
        else {
            
        }
        return (seq_len, mers, mers_rev);
    };
    let ref_process_read_aux_mer = |ref_str: &[u8], ref_id: &str| -> Option<u64> {
        let (ref_len, mers, mers_rev) = seq_to_kmers(ref_str, ref_id, params, false, &query_mers_index, (fmin, fmax));
        graph::add_kminmers(ref_id, &mers, &dbg_nodes, &NODE_INDEX, &dbg_edges);
        lens.insert(ref_id.to_string(), ref_str.len());          

        println!("Indexed {} unique reference k-mers.", mers.len());
        return Some(1)
    
    };
    let ref_process_read_fasta_mer = |record: seq_io::fasta::RefRecord, found: &mut Option<u64>| {
        let ref_str = record.seq(); 
        let ref_id = record.id().unwrap().to_string();
        *found = ref_process_read_aux_mer(&ref_str, &ref_id);
    };
    let ref_process_read_fastq_mer = |record: seq_io::fastq::RefRecord, found: &mut Option<u64>| {
        let ref_str = record.seq(); 
        let ref_id = record.id().unwrap().to_string();
        *found = ref_process_read_aux_mer(&ref_str, &ref_id);
    };
    let ref_main_thread_mer = |found: &mut Option<u64>| { // runs in main thread
        let seq_id_hash = found.as_ref().unwrap();
        None::<()>
    };
    let query_process_read_aux_mer = |seq_str: &[u8], seq_id: &str| -> Option<(Match, String)> {
       // println!("{:?}", seq_id);
        let output = mers::new_query_graph(&seq_id, seq_str, &params, &dbg_nodes, &dbg_edges, &lens, false, &query_mers_index, (fmin, fmax));
        let output_rev = mers::new_query_graph(&seq_id, seq_str, &params, &dbg_nodes, &dbg_edges, &lens, true, &query_mers_index, (fmin, fmax));
        if output.8 > output_rev.8 {return Some((output, seq_id.to_string()));}
        else if output.8 < output_rev.8 {return Some((output_rev, seq_id.to_string()));} 
        else {
           if (output.5 as i32 - output.4 as i32) > (output_rev.5 as i32 - output_rev.4 as i32) {
               if output.5 != 0 {
                return Some((output, seq_id.to_string()));
               }
               else {return None;}
            }
            else {
                if output_rev.5 != 0 {
                    return Some((output_rev, seq_id.to_string()));
                }
                else {return None;}
            }
        }

    };
    let query_process_read_fasta_mer = |record: seq_io::fasta::RefRecord, found: &mut Option<(Match, String)>| {
        let seq_str = record.seq(); 
        let seq_id = record.id().unwrap().to_string();
        *found = query_process_read_aux_mer(&seq_str, &seq_id);
    
    };
    
    let query_process_read_fastq_mer = |record: seq_io::fastq::RefRecord, found: &mut Option<(Match, String)>| {
        let seq_str = record.seq(); 
        let seq_id = record.id().unwrap().to_string();
        *found = query_process_read_aux_mer(&seq_str, &seq_id);
    };
    let mut main_thread_mer = |found: &mut Option<(Match, String)>| { // runs in main thread
        if found.is_some() {
            let (m, seq_id) = found.as_ref().unwrap();
            all_matches.push((m.clone(), seq_id.to_string()));
        }
        None::<()>
    };
    
    let buf = get_reader(&ref_filename);
    if ref_fasta_reads {
        let reader = seq_io::fasta::Reader::new(buf);
        read_process_fasta_records(reader, threads as u32, queue_len, ref_process_read_fasta_mer, |record, found| {ref_main_thread_mer(found)});
    }
    else {
        let reader = seq_io::fastq::Reader::new(buf);
        read_process_fastq_records(reader, threads as u32, queue_len, ref_process_read_fastq_mer, |record, found| {ref_main_thread_mer(found)});
    }
    let buf = get_reader(&filename);    
    if fasta_reads {
        let reader = seq_io::fasta::Reader::new(buf);
        read_process_fasta_records(reader, threads as u32, queue_len, query_process_read_fasta_mer, |record, found| {main_thread_mer(found)});
        mers::output_paf(&mut all_matches, &mut paf_file);
    }
    else {
        let reader = seq_io::fastq::Reader::new(buf);
        read_process_fastq_records(reader, threads as u32, queue_len, query_process_read_fastq_mer, |record, found| {main_thread_mer(found)});
        mers::output_paf(&mut all_matches, &mut paf_file);
    }
}