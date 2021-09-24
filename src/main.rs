// hifimap v0.1.0
// Copyright 2020-2021 Baris Ekim, , Kristoffer Sahlin, Rayan Chikhi.
// Licensed under the MIT license (http://opensource.org/licenses/MIT).
// This file may not be copied, modified, or distributed except according to those terms.

#![allow(unused_variables)]
#![allow(non_upper_case_globals)]
#![allow(warnings)]
use indicatif::ProgressBar;
use std::io::stderr;
use std::error::Error;
use std::io::Write;
use std::io::{BufRead, BufReader};
use std::collections::HashMap;
use std::fs::{File};
use std::collections::HashSet;
extern crate array_tool;
use std::fs;
use structopt::StructOpt;
use std::sync::{Arc, Mutex};
use std::path::PathBuf;
use std::time::{Instant};
use std::mem::{MaybeUninit};
use seq_io::BaseRecord;
use lzzzz::lz4f::{WriteCompressor, BufReadDecompressor, Preferences};
use xx_bloomfilter::Bloom;
use flate2::read::GzDecoder;
use dashmap::DashMap;
use std::cell::UnsafeCell;
use std::io::Result;
use core::hash::{Hash, Hasher};
use std::collections::hash_map::DefaultHasher;
mod utils;
mod paf_output;
mod seq_output;
mod kmer_vec;
mod closures;
mod mers;
pub type MersVectorRead = Vec<(u64, u64, usize, usize, bool)>;
pub type MersVectorReadMod = Vec<(Vec<u64>, u64, usize, usize, bool)>;
pub type MersVector = Vec<(u64, u64, usize, usize)>;
pub type MersVectorReduced = Vec<(u64, usize, usize)>;
pub type KmerLookup = DashMap<u64, (usize, usize)>;
pub type KmerLookupMod = DashMap<u64, (u64, usize)>;
pub type KmerLookupModMod = DashMap<Vec<u64>, (u64, usize)>;
pub type PosIndex = DashMap<u64, MersVectorRead>;
pub type PosIndexMod = DashMap<u64, MersVectorReadMod>;
const revcomp_aware : bool = true; // shouldn't be set to false except for strand-directed data or for debugging
type Kmer = kmer_vec::KmerVec;
type Overlap = kmer_vec::KmerVec;
#[derive(Clone, Debug)] // seems necessary to move out of the Arc into dbg_nodes_view
pub struct HashEntry {origin: String, seq: String, seqlen: u32, shift: (usize, usize), seq_rev: bool}
#[derive(Clone, Debug)] // seems necessary to move out of the Arc into dbg_nodes_view
pub struct QueryMatch {kminmer: Kmer, target: Vec<HashEntry>, query: HashEntry}
pub struct SeqFileType(WriteCompressor<File>);
unsafe impl Sync for SeqFileType {} // same trick as below. we won't share files among threads but Rust can't know that.
impl SeqFileType {
    fn new(v: File) -> SeqFileType {
        SeqFileType(WriteCompressor::new(v, Preferences::default()).unwrap())
    }
}
impl Write for SeqFileType {
    fn write(&mut self, buf: &[u8]) -> Result<usize> {
        self.0.write(buf)
    }

    fn flush(&mut self) -> Result<()> {
        self.0.flush()
    }
}
type ThreadIdType = usize;
pub struct Params {
    l: usize,
    k: usize,
    density: f64,
    s: usize,
    use_hpc: bool,
    f: f64,
    wmin: usize,
    wmax: usize,
    use_strobe: bool,
}

pub fn hash_id<T: Hash>(t: &T) -> u64 {
    let mut s = DefaultHasher::new();
    t.hash(&mut s);
    s.finish()
}

fn debug_output_read_minimizers(seq_str: &str, read_minimizers: &[String], read_minimizers_pos: &[u32]) {
    println!("\nseq: {}", seq_str);
    print!("min: ");
    let mut current_minimizer :String = "".to_string();
    for i in 0..seq_str.len() {
        if read_minimizers_pos.contains(&(i as u32)) {
            let index = read_minimizers_pos.iter().position(|&r| r == i as u32).unwrap();
            current_minimizer = read_minimizers[index].clone();
            let c = current_minimizer.remove(0);
            if c == seq_str.chars().nth(i).unwrap() {print!("X");}
            else {print!("x");}
            continue;
        }
        if !current_minimizer.is_empty() {
            let c = current_minimizer.remove(0);
            print!("{}", c);
        }
        else {print!(".");}
    }
    println!("");
}

/// Try to get memory usage (resident set size) in bytes using the `getrusage()` function from libc.
// from https://github.com/digama0/mm0/blob/bebd670c5a77a1400913ebddec2c6248e76f90fe/mm0-rs/src/util.rs
fn get_memory_rusage() -> usize {
  let usage = unsafe {
    let mut usage = MaybeUninit::uninit();
    assert_eq!(libc::getrusage(libc::RUSAGE_SELF, usage.as_mut_ptr()), 0);
    usage.assume_init()
  };
  usage.ru_maxrss as usize * 1024
}

// thread helpers
fn thread_update_hashmap<U, V>(hashmap_all: &Arc<Mutex<HashMap<usize, HashMap<U, V>>>>, hashmap: HashMap<U, V>, thread_num: usize) {
    let mut hashmap_all = hashmap_all.lock().unwrap();
    let entry = hashmap_all.entry(thread_num).or_insert_with(HashMap::new);
    *entry = hashmap; // I believe hashmap is moved in this function as per https://stackoverflow.com/a/29490907 
}

pub fn thread_update_vec<U>(vec_all: &Arc<Mutex<HashMap<usize, Vec<U>>>>, vec: Vec<U>, thread_num: usize) {
    let mut vec_all = vec_all.lock().unwrap();
    let entry = vec_all.entry(thread_num).or_insert_with(Vec::new);
    *entry = vec;
}

fn get_reader(path: &PathBuf) -> Box<dyn BufRead + Send> {
    let mut filetype = "unzip";
    let filename_str = path.to_str().unwrap();
    let file = match File::open(path) {
            Ok(file) => file,
            Err(error) => panic!("Error opening compressed file: {:?}.", error),
        };
    if filename_str.ends_with(".gz")  {filetype = "zip";}
    if filename_str.ends_with(".lz4") {filetype = "lz4";}
    let reader :Box<dyn BufRead + Send> = match filetype { 
        "zip" => Box::new(BufReader::new(GzDecoder::new(file))), 
        "lz4" => Box::new(BufReadDecompressor::new(BufReader::new(file)).unwrap()),
        _ =>     Box::new(BufReader::new(file)), 
    }; 
    reader
}

/*fn read_first_n_reads(filename: &PathBuf, fasta_reads: bool, max_reads: usize) -> (usize, usize) {
    let mut mean_length = 0;
    let mut max_length = 0;
    let mut nb_reads = 0;
    let buf = get_reader(&filename);
    // todo should factorize
    if fasta_reads {
        let mut reader = seq_io::fasta::Reader::new(buf);
        while let Some(record) = reader.next() {
            let record = record.unwrap();
            let seq_str = String::from_utf8_lossy(record.seq()).to_string(); // might induce a copy? can probably be optimized (see https://docs.rs/seq_io/0.4.0-alpha.0/seq_io/fasta/index.html)
            let l = seq_str.len();
            mean_length += l;
            max_length = std::cmp::max(l, max_length);
            nb_reads += 1;
            if nb_reads == max_reads {break;}
        }
    }
    else {
        let mut reader = seq_io::fastq::Reader::new(buf);
        while let Some(record) = reader.next() {
            let record = record.unwrap();
            let seq_str = String::from_utf8_lossy(record.seq()).to_string(); // might induce a copy? can probably be optimized (see https://docs.rs/seq_io/0.4.0-alpha.0/seq_io/fasta/index.html)
            let l = seq_str.len();
            mean_length += l;
            max_length = std::cmp::max(l, max_length);
            nb_reads += 1;
            if nb_reads == max_reads {break;}
        }
    }
    mean_length /= nb_reads;
    (mean_length, max_length) 
}*/

/*fn autodetect_k_l_d(filename: &PathBuf, fasta_reads: bool) -> (usize, usize, f64) {
    println!("Parsing input sequences to estimate mean read length...");
    let (mean_length, max_length) = read_first_n_reads(&filename, fasta_reads, 100);
    println!("Detected mean read length of {} bp.",mean_length);
    // a bit crude, but let's try
    let d = 0.003;
    let coeff : f64 = 3.0/4.0;
    let slightly_below_readlen : f64 = mean_length as f64;
    let k = (d * slightly_below_readlen) as usize;
    let l = 12;
    println!("Setting k = {} l = {} density = {}.", k, l, d);
    (k, l, d)
}*/

#[derive(Debug, StructOpt)]
#[structopt(name = "hifimap")]
/// Original implementation of hifimap, a fast HiFi read mapper.
struct Opt {
    /// Activate debug mode
    ///
    /// Debug mode shows the base-space sequence and 
    /// the minimizer-space representation obtained from each read.
    #[structopt(long)]
    debug: bool,
    /// Input file (raw or gzip-/lz4-compressed FASTX)
    ///
    /// Input file can be FASTA/FASTQ, as well as gzip-compressed (.gz) or
    /// lz4-compressed (.lz4). Lowercase bases are currently not supported;
    /// see documentation for formatting.
    #[structopt(parse(from_os_str))]
    reads: Option<PathBuf>,
    #[structopt(long)]
    #[structopt(long)]
    /// Output prefix for GFA and .sequences files
    ///
    /// All files generated by rust-mdbg (including the GFA
    /// and .sequences output files) will have this file name.
    #[structopt(parse(from_os_str), short, long)]
    prefix: Option<PathBuf>,
    /// k-min-mer length
    ///
    /// The length of each node of the mdBG. If
    /// fewer l-mers than this value are obtained
    /// from a read, they will be ignored.
    #[structopt(short, long)]
    k: Option<usize>,
    /// l-mer (minimizer) length
    ///
    /// The length of each minimizer selected using
    /// the minimizer scheme from base-space sequences.
    #[structopt(short, long)]
    l: Option<usize>,
    /// Density threshold for density-based selection scheme
    /// 
    /// The density threshold is analogous to the
    /// fraction of l-mers that will be selected as
    /// minimizers from a read.
    #[structopt(short, long)]
    density: Option<f64>,
    /// Syncmer length for syncmer-based selection scheme.
    #[structopt(short, long)]
    s: Option<usize>,
    /// Filter cutoff value
    ///
    /// Removes this fraction of repetitive syncmers.
    #[structopt(short, long)]
    f: Option<f64>,
    /// Syncmer-based window minimum value.
    #[structopt(short, long)]
    wmin: Option<usize>,
    /// Syncmer-based window maximum value.
    #[structopt(short, long)]
    wmax: Option<usize>,
    /// Reference genome input
    ///
    /// Reference to be indexed and mapped to. 
    /// Allows multi-line FASTA and
    /// doesn't filter any kminmers.
    #[structopt(parse(from_os_str), long)]
    reference: Option<PathBuf>,
    /// Enable strobemers.
    #[structopt(long)]
    strobe: bool,
    /// Homopolymer-compressed (HPC) input
    ///
    /// Both raw and homopolymer-compressed (HPC) reads can
    /// be provided as input. If the reads are not compressed,
    /// rust-mdbg manually performs HPC, but uses the raw sequences
    /// for transformation into base-space.
    #[structopt(long)]
    hpc: bool,
    /// Number of threads
    ///
    /// rust-mdbg is highly parallelized to decrease running
    /// time, but can be run on a single core as well.
    #[structopt(long)]
    threads: Option<usize>,

}

fn main() {
    let start = Instant::now();
    let opt = Opt::from_args();      
    let mut filename = PathBuf::new();
    let mut ref_filename = PathBuf::new();
    let mut output_prefix;
    let mut k : usize = 10;
    let mut l : usize = 12;
    let mut wmax : usize = 0;
    let mut wmin : usize = 0;
    let mut f : f64 = 0.0;
    let mut s : usize = 0;
    let mut density : f64 = 0.10;
    let reference : bool = false;
    let mut use_strobe : bool = false;
    let mut use_hpc : bool = false;
    let mut threads : usize = 8;
    if opt.reads.is_some() {filename = opt.reads.unwrap();} 
    if opt.reference.is_some() {ref_filename = opt.reference.unwrap();} 
    if filename.as_os_str().is_empty() {panic!("Please specify an input file.");}
    if ref_filename.as_os_str().is_empty() {panic!("Please specify a reference file.");}
    let mut fasta_reads : bool = false;
    let mut ref_fasta_reads : bool = false;
    let filename_str = filename.to_str().unwrap();
    if filename_str.contains(".fasta.") || filename_str.contains(".fa.") || filename_str.ends_with(".fa") || filename_str.ends_with(".fasta") { // not so robust but will have to do for now
        fasta_reads = true;
        println!("Input file: {}", filename_str);
        println!("Format: FASTA");
    }
    let ref_filename_str = ref_filename.to_str().unwrap();
    if ref_filename_str.contains(".fasta.") || ref_filename_str.contains(".fa.") || ref_filename_str.ends_with(".fa") || ref_filename_str.ends_with(".fasta") { // not so robust but will have to do for now
        ref_fasta_reads = true;
        println!("Reference file: {}", ref_filename_str);
        println!("Format: FASTA");
    }
    if opt.k.is_some() {k = opt.k.unwrap()} else {println!("Warning: Using default k value ({}).", k);} 
    if opt.l.is_some() {l = opt.l.unwrap()} else {println!("Warning: Using default l value ({}).", l);}
    if opt.density.is_some() {density = opt.density.unwrap()} else if opt.s.is_none() {println!("Warning: Using default density value ({}%).", density * 100.0);}
    if opt.threads.is_some() {threads = opt.threads.unwrap();} else {println!("Warning: Using default number of threads (8).");}
    if opt.strobe {use_strobe = true;}
    if opt.hpc {use_hpc = true;}
    output_prefix = PathBuf::from(format!("hifimap-k{}-d{}-l{}", k, density, l));
    if opt.f.is_some() {f = opt.f.unwrap()} else {println!("Warning: Not using filtering.");} 
    if opt.s.is_some() { 
        s = opt.s.unwrap(); 
        println!("Syncmer parameter s: {}", s);
        if opt.wmin.is_some() {wmin = opt.wmin.unwrap()} else {println!("Warning: Using default wmin value ({}).", l/(l-s+1)+2);} 
        if opt.wmax.is_some() {wmax = opt.wmax.unwrap()} else {println!("Warning: Using default wmax value ({}).", l/(l-s+1)+10);} 

    } 
    if opt.prefix.is_some() {output_prefix = opt.prefix.unwrap();} else {println!("Warning: Using default output prefix ({}).", output_prefix.to_str().unwrap());}
    let params = Params { 
        l,
        k,
        density,
        s,
        use_hpc,
        f,
        wmin,
        wmax,
        use_strobe,
    };
    // init some useful objects
    // get file size for progress bar
    let metadata = fs::metadata(&filename).expect("Error opening input file.");
    let ref_metadata = fs::metadata(&ref_filename).expect("Error opening reference file.");
    let file_size = metadata.len();
    let queue_len = threads; // https://doc.rust-lang.org/std/sync/mpsc/fn.sync_channel.html
                             // also: controls how many reads objects are buffered during fasta/fastq
                             // parsing

    //let mut bloom : RacyBloom = RacyBloom::new(Bloom::new_with_rate(if use_bf {100_000_000} else {1}, 1e-7)); // a bf to avoid putting stuff into kmer_table too early
    let paf_path = format!("{}{}", output_prefix.to_str().unwrap(), ".paf");
    let paf_file = match File::create(&paf_path) {
        Err(why) => panic!("Couldn't create {}: {}", paf_path, why.description()),
        Ok(paf_file) => paf_file,
    };
    closures::run_mers(&filename, &ref_filename, &params, threads, queue_len, fasta_reads, ref_fasta_reads, &output_prefix);
    let duration = start.elapsed();
    println!("Total execution time: {:?}", duration);
    println!("Maximum RSS: {:?}GB", (get_memory_rusage() as f32) / 1024.0 / 1024.0 / 1024.0 / 1024.0);
}



