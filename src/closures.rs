// closures.rs
// Functions for FASTA parsing and invoking all necessary functions for mapping and alignment.

use std::error::Error;
use std::io::{Write, BufWriter};
use std::fs::{File};
use seq_io::BaseRecord;
use seq_io::parallel::{read_process_fasta_records, read_process_fastq_records, read_process_fastx_records};
use dashmap::DashMap;
use super::mers;
use std::path::{Path, PathBuf};
use super::Params;
use crate::get_reader;
use std::time::Instant;
use crate::index::{Index, ReadOnlyIndex};
//use crate::align::{get_slices, align_slices, AlignStats};
use std::sync::atomic::{AtomicUsize, Ordering};
use rust_parallelfastx::parallel_fastx;
use std::sync::mpsc;

// Main function for all FASTA parsing + mapping / alignment functions.
pub fn run_mers(filename: &PathBuf, ref_filename: &PathBuf, params: &Params, ref_threads: usize, threads: usize, ref_queue_len: usize, queue_len: usize, fasta_reads: bool, ref_fasta_reads: bool, output_prefix: &Path) {

    let mers_index = Index::new(); // Index of reference k-min-mer entries
    //let mut aln_coords : Arc<DashMap<String, Vec<AlignCand>>> =  Arc::new(DashMap::new()); // Index of AlignCand objects (see mers.rs for a definition) per reference
    //let mut aln_coords_q : Arc<DashMap<String, Vec<Offset>>> =  Arc::new(DashMap::new()); // Index of intervals that need to be aligned per query
    //let mut aln_seqs_cow : Arc<DashMap<(String, Offset), Cow<[u8]>>> =  Arc::new(DashMap::new()); // Index of pointers to string slices that need to be aligned per reference
    let ref_i = AtomicUsize::new(0);
    let ref_map : DashMap<usize, (String, usize)> = DashMap::new(); // Sequence lengths per reference

    // PAF file generation
    let paf_filename = format!("{}{}", output_prefix.to_str().unwrap(), ".paf");
    let mut paf_file = match File::create(&paf_filename) {
        Err(why) => panic!("Couldn't create {}: {}", paf_filename, why.description()),
        Ok(paf_file) => BufWriter::new(paf_file),
    };

    // Unmapped read file generation
    /*let unmap_path = format!("{}{}", output_prefix.to_str().unwrap(), ".unmapped.out");
    let unmap_file = match File::create(&unmap_path) {
        Err(why) => panic!("Couldn't create {}: {}", unmap_path, why.to_string()),
        Ok(unmap_file) => unmap_file,
    };*/

    // Closure for indexing reference k-min-mers
    let index_mers = |seq_id: &str, seq: &[u8], params: &Params| -> usize {
        let ref_idx = ref_i.fetch_add(1, Ordering::Relaxed);
        let nb_mers = mers::ref_extract(ref_idx, seq, params, &mers_index);
        ref_map.insert(ref_idx, (seq_id.to_string(), seq.len()));
        nb_mers
    };

    // Closures for obtaining k-min-mers from references

    let ref_process_read_aux_mer = |ref_str: &[u8], ref_id: &str| -> Option<u64> {
        let nb_mers = index_mers(ref_id, ref_str, params);
        //if params.a {aln_coords.insert(ref_id.to_string(), Vec::new());}
        println!("Indexed reference {}: {} k-min-mers.", ref_id, nb_mers);
        Some(1)
    };

    let ref_process_read_fasta_mer = |record: seq_io::fasta::RefRecord, found: &mut Option<u64>| {
        let ref_str = record.seq().to_ascii_uppercase(); 
        let ref_id = record.id().unwrap();
        *found = ref_process_read_aux_mer(&ref_str, ref_id);

    };
    let ref_process_read_fastq_mer = |record: seq_io::fastq::RefRecord, found: &mut Option<u64>| {
        let ref_str = record.seq().to_ascii_uppercase(); 
        let ref_id = record.id().unwrap();
        *found = ref_process_read_aux_mer(&ref_str, ref_id);
    };
    let ref_main_thread_mer = |_found: &mut Option<u64>| { // runs in main thread
        None::<()>
    };

    //
    
    // Start processing references

    let start = Instant::now();
    let (buf,_dontcare) = get_reader(ref_filename);
    if ref_fasta_reads {
        let reader = seq_io::fasta::Reader::with_capacity(buf, 64*1024*params.b);
        read_process_fasta_records(reader, ref_threads as u32, ref_queue_len, ref_process_read_fasta_mer, |_record, found| {ref_main_thread_mer(found)}).ok();
    }
    else {
        let reader = seq_io::fastq::Reader::with_capacity(buf, 64*1024*params.b);
        read_process_fastq_records(reader, ref_threads as u32, ref_queue_len, ref_process_read_fastq_mer, |_record, found| {ref_main_thread_mer(found)}).ok();
    }
    let duration = start.elapsed();
    println!("Indexed {} unique k-min-mers in {:?}.", mers_index.get_count(), duration);

    let mers_index = ReadOnlyIndex::new(mers_index.index);

    // Done, start processing queries
    
    // Closures for mapping queries to references

    let query_process_read_aux_mer = |seq_str: &[u8], seq_id: &str| -> (String, Option<String>) {
        //if params.a {aln_coords_q.insert(seq_id.to_string(), vec![]);}
        let match_opt = mers::find_matches(seq_id, seq_str.len(), seq_str, &ref_map, &mers_index, params); //&aln_coords);
        (seq_id.to_string(), match_opt)
    };
    let query_process_read_fasta_mer = |record: seq_io::fasta::RefRecord, found: &mut (String, Option<String>)| {
        let seq_str = record.seq().to_ascii_uppercase(); 
        let seq_id = record.id().unwrap();
        *found = query_process_read_aux_mer(&seq_str, seq_id);

    };
    let query_process_read_fastq_mer = |record: seq_io::fastq::RefRecord, found: &mut (String, Option<String>)| {
        let seq_str = record.seq().to_ascii_uppercase(); 
        let seq_id = record.id().unwrap();
        *found = query_process_read_aux_mer(&seq_str, seq_id);
    };
    // main thread for writing to PAF
    let mut main_thread_mer = |found: &mut (String, Option<String>)| { // runs in main thread
        let (_, match_opt) = found;
        if let Some(l) = match_opt {
            write!(paf_file, "{}\n", l).expect("Error writing line.");
        }
        None::<()>
    };

    //

    // Closures for base-level alignment (obtaining reference sequence slices)
    /*
    let ref_process_read_aux_aln = |ref_str: &[u8], ref_id: &str| -> Option<u64> {
        get_slices(ref_id, ref_str, &aln_coords, &aln_coords_q, &aln_seqs_cow);
        return Some(1)
    };
    let ref_process_read_fasta_aln = |record: seq_io::fasta::RefRecord, found: &mut Option<u64>| {
        let ref_str = record.seq(); 
        let ref_id = record.id().unwrap().to_string();
        *found = ref_process_read_aux_aln(&ref_str, &ref_id);
    };
    let ref_process_read_fastq_aln = |record: seq_io::fastq::RefRecord, found: &mut Option<u64>| {
        let ref_str = record.seq(); 
        let ref_id = record.id().unwrap().to_string();
        *found = ref_process_read_aux_aln(&ref_str, &ref_id);
    };

    let ref_main_thread_aln = |found: &mut Option<u64>| { // runs in main thread
        None::<()>
    };

    //

    //  Closures for base-level alignment (obtaining query sequence slices and alignment)
    
    let query_process_read_aux_aln = |seq_str: &[u8], seq_id: &str| -> Option<AlignStats> {
        let align_stats = align_slices(seq_id, seq_str, &aln_coords_q, &aln_seqs_cow);
        return Some(align_stats)
    };
    let query_process_read_fasta_aln = |record: seq_io::fasta::RefRecord, found: &mut Option<AlignStats>| {
        let seq_str = record.seq(); 
        let seq_id = record.id().unwrap().to_string();
        *found = query_process_read_aux_aln(&seq_str, &seq_id);
        //if params.a {eprintln!("{}", seq_id);}
    
    };
    let query_process_read_fastq_aln = |record: seq_io::fastq::RefRecord, found: &mut Option<AlignStats>| {
        let seq_str = record.seq(); 
        let seq_id = record.id().unwrap().to_string();
        *found = query_process_read_aux_aln(&seq_str, &seq_id);
    };
    let mut global_align_stats = AlignStats{ successful: 0, failed: 0 };
    let mut main_thread_aln = |align_stats: &mut Option<AlignStats>| { // runs in main thread
        global_align_stats.successful += align_stats.as_ref().unwrap().successful;
        global_align_stats.failed     += align_stats.as_ref().unwrap().failed;
        None::<()>
    };

    */

    let query_start = Instant::now();
    let (buf, are_reads_compressed) = get_reader(filename);
    if are_reads_compressed || (!params.use_pfx) {  // fall-back to seq_io parallel
        // spawn read processing threads
        if fasta_reads {
            let reader = seq_io::fasta::Reader::with_capacity(buf, 64*1024*params.b);
            let _ = read_process_fasta_records(reader, threads as u32, queue_len, query_process_read_fasta_mer, |record, found| {main_thread_mer(found)});
        }
        else {
            let reader = seq_io::fastq::Reader::with_capacity(buf, 64*1024*params.b);
            let _ = read_process_fastq_records(reader, threads as u32, queue_len, query_process_read_fastq_mer, |record, found| {main_thread_mer(found)});
        }
    } else { // rust-parallelfastx is a more efficient fastx parser than seq_io when reading from disk is fast and file is uncompressed
        // the only downside is that it will display a large RSS footprint as the reads will be
        // loaded in memory (though, that memory isn't needed by hifimap, it will just use as much as possible)
        println!("Warning: using experimental rust-parallelfastx (exciting!)");
        let (paf_mpsc_send, paf_mpsc_recv) = mpsc::sync_channel(1000);
        let task = |seq_str: &[u8], seq_id: &str|  {
            let (_osef, match_opt) = query_process_read_aux_mer(seq_str, seq_id);
            if let Some(l) = match_opt {
                paf_mpsc_send.send(Some(l));
            }
        };
        let writer = std::thread::spawn(move || {
            while let Some(paf_line) = paf_mpsc_recv.recv().unwrap() {
                write!(paf_file, "{}\n", paf_line).expect("Error writing line.");
            }
        });
        parallel_fastx(&filename.to_string_lossy(), threads, task);
        paf_mpsc_send.send(None); // signal we're done
        writer.join().unwrap();
    }

    let query_duration = query_start.elapsed();
    println!("Mapped query sequences in {:?}.", query_duration);



    // Done, start alignment (optional)
    /*
    if params.a {
        let start = Instant::now();
        let buf_aln = get_reader(&ref_filename);

        // Obtain reference sequences

        if ref_fasta_reads {
            let reader = seq_io::fasta::Reader::new(buf_aln);
            read_process_fasta_records(reader, ref_threads as u32, 4, ref_process_read_fasta_aln, |record, found| {ref_main_thread_aln(found)});
        }
        else {
            let reader = seq_io::fastq::Reader::new(buf_aln);
            read_process_fastq_records(reader, ref_threads as u32, 4, ref_process_read_fastq_aln, |record, found| {ref_main_thread_aln(found)});
        }
        let duration = start.elapsed();
        println!("Obtained references for alignment in {:?}.", duration);

        // Done, obtain query sequences and align

        let query_start = Instant::now();
        let buf_aln = get_reader(&filename);
        if fasta_reads {
            let reader = seq_io::fasta::Reader::new(buf_aln);
            read_process_fasta_records(reader, threads as u32, threads, query_process_read_fasta_aln, |record, found| {main_thread_aln(found)});
            let query_duration = query_start.elapsed();
            println!("Aligned query sequences in {:?}.", query_duration);
        }
        else {
            let reader = seq_io::fastq::Reader::new(buf_aln);
            read_process_fastq_records(reader, threads as u32, threads, query_process_read_fastq_aln, |record, found| {main_thread_aln(found)});
            let query_duration = query_start.elapsed();
            println!("Aligned query sequences in {:?}.", query_duration);
        }
        if global_align_stats.failed > 0
        {
            println!("[warning] Alignment stats: {} successful, {} failed", global_align_stats.successful, global_align_stats.failed);
        }
    }
    */
    //println!("current time before exiting closures {:?}",Utc::now());
}
