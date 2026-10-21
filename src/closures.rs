// closures.rs
// Functions for FASTA parsing and invoking all necessary functions for mapping and alignment.

use std::borrow::Cow;
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
use crate::report::{self, CoverageAccumulator, PafRecord};
use std::time::Instant;
use crate::index::{Index, ReadOnlyIndex};
//use crate::align::{get_slices, align_slices, AlignStats};
use std::sync::atomic::{AtomicUsize, Ordering};
use rust_parallelfastx::parallel_fastx;
use std::sync::mpsc;

/// Returns `seq` with all line breaks (`\n`, `\r`) removed and ASCII bases uppercased.
///
/// `seq_io`'s `BaseRecord::seq()` returns the raw record bytes for FASTA records, i.e.
/// including embedded line breaks when the record spans multiple lines (and
/// `rust-parallelfastx` does the same), so sequences coming from those paths must be
/// normalised before k-min-mer extraction / matching. No allocation happens when the
/// input is already clean and uppercase.
pub fn normalize_seq(seq: &[u8]) -> Cow<[u8]> {
    if !seq.iter().any(|&b| b == b'\n' || b == b'\r' || b.is_ascii_lowercase()) {
        return Cow::Borrowed(seq);
    }
    let mut out = Vec::with_capacity(seq.len());
    out.extend(seq.iter().copied().filter(|&b| b != b'\n' && b != b'\r'));
    out.make_ascii_uppercase();
    Cow::Owned(out)
}

// Main function for all FASTA parsing + mapping / alignment functions.
pub fn run_mers(filename: &PathBuf, ref_filename: &PathBuf, params: &Params, ref_threads: usize, threads: usize, ref_queue_len: usize, queue_len: usize, fasta_reads: bool, ref_fasta_reads: bool, output_prefix: &Path, json_out: Option<&Path>) {

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

    // Coverage report accumulation (only when --json was given). Contigs are
    // read once more, sequentially, to know their order in the reference file:
    // reference ids are handed out from an atomic counter in worker threads.
    let mut report_accumulator: Option<CoverageAccumulator> = match json_out {
        Some(_) => {
            let contigs = report::contigs_in_file_order(ref_filename, ref_fasta_reads)
                .unwrap_or_else(|e| panic!("Error reading reference {} for the JSON report: {}", ref_filename.to_str().unwrap(), e));
            Some(CoverageAccumulator::new(contigs))
        }
        None => None,
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
        let full_seq = record.full_seq(); // seq() is not joined for multi-line FASTA, full_seq() is
        let ref_str = normalize_seq(&full_seq);
        let ref_id = record.id().unwrap();
        *found = ref_process_read_aux_mer(ref_str.as_ref(), ref_id);

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

    let query_process_read_aux_mer = |seq_str: &[u8], seq_id: &str| -> (String, usize, Option<PafRecord>) {
        //if params.a {aln_coords_q.insert(seq_id.to_string(), vec![]);}
        let match_opt = mers::find_matches(seq_id, seq_str.len(), seq_str, &ref_map, &mers_index, params); //&aln_coords);
        (seq_id.to_string(), seq_str.len(), match_opt)
    };
    let query_process_read_fasta_mer = |record: seq_io::fasta::RefRecord, found: &mut (String, usize, Option<PafRecord>)| {
        let full_seq = record.full_seq(); // seq() is not joined for multi-line FASTA, full_seq() is
        let seq_str = normalize_seq(&full_seq);
        let seq_id = record.id().unwrap();
        *found = query_process_read_aux_mer(seq_str.as_ref(), seq_id);

    };
    let query_process_read_fastq_mer = |record: seq_io::fastq::RefRecord, found: &mut (String, usize, Option<PafRecord>)| {
        let seq_str = record.seq().to_ascii_uppercase(); 
        let seq_id = record.id().unwrap();
        *found = query_process_read_aux_mer(&seq_str, seq_id);
    };
    // main thread for writing to PAF (and collecting the JSON report data)
    let mut main_thread_mer = |found: &mut (String, usize, Option<PafRecord>)| { // runs in main thread
        let (_, q_len, match_opt) = found;
        if let Some(acc) = report_accumulator.as_mut() {
            acc.add_input_read(*q_len as u64);
        }
        if let Some(record) = match_opt {
            if let Some(acc) = report_accumulator.as_mut() {
                acc.add_paf_record(record);
            }
            writeln!(paf_file, "{}", record.to_line()).expect("Error writing line.");
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
        let ref_str = record.full_seq(); // multi-line FASTA safe
        let ref_id = record.id().unwrap().to_string();
        *found = ref_process_read_aux_aln(ref_str.as_ref(), &ref_id);
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
        let seq_str = record.full_seq(); // multi-line FASTA safe
        let seq_id = record.id().unwrap().to_string();
        *found = query_process_read_aux_aln(seq_str.as_ref(), &seq_id);
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
        let (paf_mpsc_send, paf_mpsc_recv) = mpsc::sync_channel::<Option<PafRecord>>(1000);
        let input_reads = AtomicUsize::new(0);
        let input_bases = AtomicUsize::new(0);
        let task = |seq_str: &[u8], seq_id: &str|  {
            // rust-parallelfastx passes the raw record bytes, which still contain line
            // breaks for multi-line records, so normalise before mapping.
            let seq_str = normalize_seq(seq_str);
            input_reads.fetch_add(1, Ordering::Relaxed);
            input_bases.fetch_add(seq_str.len(), Ordering::Relaxed);
            let (_osef, _q_len, match_opt) = query_process_read_aux_mer(seq_str.as_ref(), seq_id);
            if let Some(record) = match_opt {
                let _ = paf_mpsc_send.send(Some(record));
            }
        };
        let report_accumulator_for_writer = report_accumulator.take();
        let writer = std::thread::spawn(move || {
            let mut accumulator = report_accumulator_for_writer;
            while let Some(record) = paf_mpsc_recv.recv().unwrap() {
                if let Some(acc) = accumulator.as_mut() {
                    acc.add_paf_record(&record);
                }
                writeln!(paf_file, "{}", record.to_line()).expect("Error writing line.");
            }
            accumulator
        });
        parallel_fastx(&filename.to_string_lossy(), threads, task);
        let _ = paf_mpsc_send.send(None); // signal we're done
        report_accumulator = writer.join().unwrap();
        if let Some(acc) = report_accumulator.as_mut() {
            acc.set_input_totals(
                input_reads.load(Ordering::Relaxed) as u64,
                input_bases.load(Ordering::Relaxed) as u64,
            );
        }
    }

    let query_duration = query_start.elapsed();
    println!("Mapped query sequences in {:?}.", query_duration);

    // Write the JSON coverage report (only if --json was given).
    if let (Some(json_path), Some(accumulator)) = (json_out, report_accumulator.as_ref()) {
        let meta = report::ReportMeta {
            version: env!("CARGO_PKG_VERSION").to_string(),
            reads: filename.display().to_string(),
            reference: ref_filename.display().to_string(),
            parameters: report::Parameters {
                k: params.k,
                l: params.l,
                density: params.density,
                chain: params.c,
                seed: params.s,
                gap_diff: params.g,
                b: params.b,
                q: params.q,
                threads,
                low_memory: params.low_memory,
                hpc: params.use_hpc,
                simd: params.use_simd,
                parallelfastx: params.use_pfx,
                debug: params.debug,
            },
        };
        if let Err(e) = report::write_report(json_path, accumulator, &meta) {
            panic!("Error writing JSON report to {}: {}", json_path.display(), e);
        }
        println!("Wrote JSON coverage report to {}", json_path.display());
    }



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

#[cfg(test)]
mod tests {
    use super::normalize_seq;
    use seq_io::BaseRecord;

    // Two records: the first spans three lines (last one shorter, mixed case) and
    // contains a description, the second uses CRLF line endings.
    const MULTILINE_FASTA: &[u8] = b">seq1 some description\nacgtACGT\ntttt\nAa\n>seq2\r\nGGGG\r\ncccc\r\n";

    // Same parsing path as `query_process_read_fasta_mer` / `ref_process_read_fasta_mer`.
    fn read_fasta(bytes: &[u8]) -> Vec<(String, Vec<u8>)> {
        let mut reader = seq_io::fasta::Reader::new(bytes);
        let mut records = Vec::new();
        while let Some(record) = reader.next() {
            let record = record.expect("malformed FASTA");
            let full_seq = record.full_seq();
            let seq = normalize_seq(&full_seq);
            records.push((record.id().unwrap().to_string(), seq.into_owned()));
        }
        records
    }

    #[test]
    fn multiline_fasta_sequences_are_joined() {
        let records = read_fasta(MULTILINE_FASTA);
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].0, "seq1");
        assert_eq!(records[0].1, b"ACGTACGTTTTTAA");
        assert_eq!(records[1].0, "seq2");
        assert_eq!(records[1].1, b"GGGGCCCC");
    }

    #[test]
    fn raw_seq_still_contains_line_breaks_for_multiline_records() {
        // Documents why `full_seq()` + normalization is needed: `seq()` is not joined.
        let mut reader = seq_io::fasta::Reader::new(MULTILINE_FASTA);
        let record = reader.next().unwrap().unwrap();
        assert_eq!(record.seq(), b"acgtACGT\ntttt\nAa");
        assert_eq!(record.full_seq().as_ref(), b"acgtACGTttttAa");
    }

    #[test]
    fn normalize_seq_strips_line_breaks_and_uppercases() {
        assert_eq!(normalize_seq(b"acgtACGT\ntttt\nAa").as_ref(), b"ACGTACGTTTTTAA");
        assert_eq!(normalize_seq(b"acgt\r\nACGT\r\n").as_ref(), b"ACGTACGT");
        assert_eq!(normalize_seq(b"ACGT").as_ref(), b"ACGT");
        assert_eq!(normalize_seq(b"").as_ref(), b"");
        // Already clean input is borrowed, not copied.
        assert!(matches!(normalize_seq(b"ACGT"), std::borrow::Cow::Borrowed(_)));
    }
}
