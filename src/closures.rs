use std::io::{self};
use std::error::Error;
use std::io::{BufRead, BufReader};
use std::path::Path;
use crate::BufReadDecompressor;
use std::fs::{File};
use std::sync::{Arc};
use seq_io::BaseRecord;
use seq_io::parallel::{read_process_fasta_records, read_process_fastq_records};
use dashmap::DashMap;
use thread_id;
use super::mers;
use crate::mers::{Match, Mer};
use super::{MersVector, MersVectorReduced, MersVectorRead, PosIndex, KmerLookup, KmerLookupMod, PosIndexMod};
use std::path::PathBuf;
use super::Params;
use crate::{get_reader, hash_id};

pub fn run_mers(filename: &PathBuf, ref_filename: &PathBuf, params: &Params, threads: usize, queue_len: usize, fasta_reads: bool, ref_fasta_reads: bool, output_prefix: &PathBuf) {
    let mers_index : Arc<DashMap<u64, Vec<(String, usize)>>> =  Arc::new(DashMap::new());
    let mut all_matches = Vec::<(Vec<Match>, String)>::new();
    let path = format!("{}{}", output_prefix.to_str().unwrap(), ".paf");
    let mut paf_file = match File::create(&path) {
        Err(why) => panic!("Couldn't create {}: {}", path, why.description()),
        Ok(paf_file) => paf_file,
    };
    let ref_mers : DashMap<String, (usize, Vec<Mer>)> = DashMap::new();
    let query_mers : DashMap<String, (usize, Vec<Mer>)> = DashMap::new();
    let index_mers = |seq_id: &str, seq: &[u8], params: &Params, read: bool| -> (usize, Vec<Mer>) {
        let mut unique_mers = Vec::<Mer>::new();
        let (seq_len, mers) = mers::seq_to_kmers(seq, seq_id, params, read);
        
        if read == true {
            println!("Seq ID\t{}\tmers len\t{}", seq_id, mers.len());
            let mut i_p = 0;
            for i in 0..mers.len() {
                let mut mer = mers[i];
                let entry = mers_index.get_mut(&mer.0);
                if entry.is_some() {
                    if entry.unwrap().len() == 1 {
                        mer.3 = i_p;
                        unique_mers.push(mer);
                        i_p += 1;
                    }
                }
            }
            query_mers.insert(seq_id.to_string(), (seq_len, unique_mers.to_vec()));
        }
        else {
            let mut i_p = 0;
            for i in 0..mers.len() {
                let mut mer = mers[i];
                let entry = mers_index.get_mut(&mer.0);
                if entry.is_some() {
                    if entry.unwrap().len() >= 1 {mers_index.remove(&mer.0);}
                }
                else {
                    mer.3 = i_p;
                    mers_index.insert(mer.0, vec![(seq_id.to_string(), i_p)]);
                    unique_mers.push(mer);
                    i_p += 1;
                }
            }  
            ref_mers.insert(seq_id.to_string(), (seq_len, unique_mers.to_vec()));   
        }
        println!("Seq ID\t{}\tunique mers len\t{}", seq_id, unique_mers.len());
        return (seq_len, unique_mers);
    };
    let ref_process_read_aux_mer = |ref_str: &[u8], ref_id: &str| -> Option<u64> {
        let (seq_len, unique_mers) = index_mers(ref_id, ref_str, params, false);
        println!("Indexed {} unique reference k-mers.", unique_mers.len());
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
        let (seq_len, mers) = index_mers(seq_id, seq_str, params, true);
        output = mers::find_hits(&seq_id, seq_len, &mers, &ref_mers, &mers_index, params.l, params.k);
        if output.len() > 0 {
            return Some((output, seq_id.to_string()))
        }
        else {return None}
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
        mers::output_paf(&mut all_matches, &mut paf_file);
    }
    else {
        let reader = seq_io::fastq::Reader::new(buf);
        read_process_fastq_records(reader, threads as u32, queue_len, query_process_read_fastq_mer, |record, found| {main_thread_mer(found)});
        mers::output_paf(&mut all_matches, &mut paf_file);
    }
}