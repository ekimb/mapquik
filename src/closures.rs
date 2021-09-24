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
use crate::mers::{Match, Mer};
use super::{MersVector, MersVectorReduced, MersVectorRead, PosIndex, KmerLookup, KmerLookupMod, PosIndexMod};
use std::path::PathBuf;
use super::Params;
use crate::{get_reader, hash_id};
use indicatif::ProgressBar;

pub fn run_mers(filename: &PathBuf, ref_filename: &PathBuf, params: &Params, threads: usize, queue_len: usize, fasta_reads: bool, ref_fasta_reads: bool, output_prefix: &PathBuf) {
    //let mers_index : Arc<DashMap<String, DashMap<u64, (Mer, usize)>>> =  Arc::new(DashMap::new());
    let mers_index : Arc<DashMap<u64, Vec<(Mer, String)>>> =  Arc::new(DashMap::new());
    let query_mers_index : Arc<DashMap<u64, usize>> =  Arc::new(DashMap::new());
    let mut all_matches = Vec::<(Vec<Match>, String)>::new();
    let path = format!("{}{}", output_prefix.to_str().unwrap(), ".paf");
    let mut threshold_max = 0;
    let mut threshold_min = u64::max_value();
    let mut paf_file = match File::create(&path) {
        Err(why) => panic!("Couldn't create {}: {}", path, why.description()),
        Ok(paf_file) => paf_file,
    };
    let lens : DashMap<String, usize> = DashMap::new();

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
    if params.f != 0.0 {
        if fasta_reads {
            let reader = seq_io::fasta::Reader::new(buf); 
            read_process_fasta_records(reader, threads as u32, queue_len, query_count_read_fasta_mer, |record, found| {count_thread_mer(found)});
        }
        else {
            let reader = seq_io::fastq::Reader::new(buf); 
            read_process_fastq_records(reader, threads as u32, queue_len, query_count_read_fastq_mer, |record, found| {count_thread_mer(found)});
        }
        let mut frac_index = (params.f * query_mers_index.len() as f64) as usize;
        let mut counts : Vec<usize> = query_mers_index.iter().map(|(ent)| *ent.value()).collect();
        counts.sort_by(|a, b| a.cmp(&b));
        let mut histogram = Histogram::with_buckets(200);
        for count in counts.iter() {
            histogram.add(*count as u64);
        }
        println!("{}", histogram);
        counts.reverse();
        threshold_max = counts.iter().nth(frac_index).unwrap().clone();
        let mut prev_count = u64::max_value();
        for bucket in histogram.buckets() {
            if  bucket.count() < prev_count {
                threshold_min = bucket.end();
                prev_count = bucket.count();
            }
            else {break;}
        }
    }
    println!("Removed top {}% fraction of {} l-mers (count > {}).", params.f*100.0, query_mers_index.len(), threshold_max);
    println!("Removed likely-erroneous l-mers based on histogram (count < {}).", threshold_min);

    let index_mers = |seq_id: &str, seq: &[u8], params: &Params, read: bool| -> (usize, Vec<Mer>) {
        let (seq_len, mut mers) = mers::seq_to_kmers(seq, seq_id, params, read, &query_mers_index, (threshold_min as usize, threshold_max as usize));
        if read == true {
            lens.insert(seq_id.to_string(), seq_len);
        }
        else {
            for mer in mers.iter() {
                let mut entry =  mers_index.get_mut(&mer.0);
                if entry.is_some() {
                    let mut e = entry.unwrap();
                    let mut match_kmer = e.value().iter().find(|x| x.1 == seq_id);
                    if match_kmer.is_some() {continue;}
                    else {
                        e.push((*mer, seq_id.to_string()));
                    }
                }
                else {mers_index.insert(mer.0, vec![(*mer,  seq_id.to_string())]);}
            }
            lens.insert(seq_id.to_string(), seq_len);
        }
        /*else {
            mers_index.insert(seq_id.to_string(), DashMap::new());
            mers.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
            let mut prev_mer = mers[0];
            let mut entry_map = mers_index.get_mut(&seq_id.to_string()).unwrap();
            entry_map.insert(prev_mer.0, (prev_mer, 1));
            for i in 1..mers.len() {
                let mut curr_mer = mers[i];
                if curr_mer.0 == prev_mer.0 {
                    let entry = entry_map.get_mut(&curr_mer.0);
                    entry.unwrap().1 += 1;
                }
                else {
                    entry_map.insert(curr_mer.0, (curr_mer, 1));
                }
                prev_mer = curr_mer;
            }  
            lens.insert(seq_id.to_string(), seq_len);
        }*/
        return (seq_len, mers);
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
        output = mers::find_hits(&seq_id, seq_len, &mers, &lens, &query_mers_index, &mers_index, params.l, params.k);
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