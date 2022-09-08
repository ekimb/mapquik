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
use crate::kminmer::Kminmer;
use crate::mers::{Match};
use std::path::PathBuf;
use super::Params;
use crate::get_reader;
use indicatif::ProgressBar;
pub type Entry = (usize, usize, bool, usize);

pub fn run_mers(filename: &PathBuf, ref_filename: &PathBuf, params: &Params, threads: usize, queue_len: usize, fasta_reads: bool, ref_fasta_reads: bool, output_prefix: &PathBuf) {
    //let mers_index : Arc<DashMap<String, DashMap<u64, (Mer, usize)>>> =  Arc::new(DashMap::new());
    let mers_index : Arc<DashMap<u64, Vec<(String, Entry)>>> =  Arc::new(DashMap::new());
    let query_mers_index : Arc<DashMap<u64, usize>> =  Arc::new(DashMap::new());
    let mut all_matches = Vec::<(Vec<Match>, String)>::new();
    let path = format!("{}{}", output_prefix.to_str().unwrap(), ".paf");
    let mut paf_file = match File::create(&path) {
        Err(why) => panic!("Couldn't create {}: {}", path, why.description()),
        Ok(paf_file) => paf_file,
    };
    let unmap_path = format!("{}{}", output_prefix.to_str().unwrap(), ".unmapped.out");
    let mut unmap_file = match File::create(&unmap_path) {
        Err(why) => panic!("Couldn't create {}: {}", unmap_path, why.description()),
        Ok(unmap_file) => unmap_file,
    };
    let lens : DashMap<String, usize> = DashMap::new();
    let buf = get_reader(&filename);
    let index_mers = |seq_id: &str, seq: &[u8], params: &Params| {
        let (mut sk, mut pos) = mers::extract(seq, params);
        mers::ref_kminmers(seq_id, &sk, &pos, params, &mers_index);
        lens.insert(seq_id.to_string(), seq.len());
    };
    let ref_process_read_aux_mer = |ref_str: &[u8], ref_id: &str| -> Option<u64> {
        index_mers(ref_id, ref_str, params);
        println!("Indexed reference k-mers for {}.", ref_id);
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
    let query_process_read_aux_mer = |seq_str: &[u8], seq_id: &str| -> Option<(Vec<Match>, String)> {
        let mut output = Vec::<Match>::new();
        output = mers::find_hits(&seq_id, seq_str.len(), &seq_str, &lens, &mers_index, params);
        return Some((output, seq_id.to_string()))
    };
    let query_process_read_fasta_mer = |record: seq_io::fasta::RefRecord, found: &mut Option<(Vec<Match>, String)>| {
        let seq_str = record.seq(); 
        let seq_id = record.id().unwrap().to_string();
        *found = query_process_read_aux_mer(&seq_str, &seq_id);
    
    };
    
    let query_process_read_fastq_mer = |record: seq_io::fastq::RefRecord, found: &mut Option<(Vec<Match>, String)>| {
        let seq_str = record.seq(); 
        let seq_id = record.id().unwrap().to_string();
        *found = query_process_read_aux_mer(&seq_str, &seq_id);
    };
    let mut main_thread_mer = |found: &mut Option<(Vec<Match>, String)>| { // runs in main thread
        if found.is_some() {
            let (matches, seq_id) = found.as_ref().unwrap();
            all_matches.push((matches.to_vec(), seq_id.to_string()));
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
        mers::output_paf(&mut all_matches, &mut paf_file, &mut unmap_file, params);
    }
    else {
        let reader = seq_io::fastq::Reader::new(buf);
        read_process_fastq_records(reader, threads as u32, queue_len, query_process_read_fastq_mer, |record, found| {main_thread_mer(found)});
        mers::output_paf(&mut all_matches, &mut paf_file, &mut unmap_file, params);
    }
}