// hifimap v0.1.0
// Copyright 2020-2021 Baris Ekim, , Kristoffer Sahlin, Rayan Chikhi.
// Licensed under the MIT license (http://opensource.org/licenses/MIT).
// This file may not be copied, modified, or distributed except according to those terms.

#![allow(unused_variables)]
#![allow(non_upper_case_globals)]
#![allow(warnings)]
use pbr::ProgressBar;
use std::io::stderr;
use std::error::Error;
use std::io::Write;
use std::io::{BufWriter, BufRead, BufReader};
use std::collections::HashMap;
use std::collections::hash_map::Entry;
use itertools::Itertools;
use ::closure;
use std::iter::FromIterator;
use crate::kmer_vec::get;
use crate::read::Read;
use crate::read::ReadSync;
use std::thread;
use crate::read::hash_id;
use std::fs::{File,remove_file};
use std::collections::HashSet;
extern crate array_tool;
use std::fs;
use structopt::StructOpt;
use std::sync::{Arc, Mutex, MutexGuard};
use std::path::PathBuf;
use std::time::{Duration, Instant};
use std::mem::{self, MaybeUninit};
use editdistancewf as wf;
use seq_io::fasta;
use seq_io::core::BufReader as OtherBufReader;
use seq_io::BaseRecord;
use seq_io::parallel::{read_process_fasta_records, read_process_fastq_records};
use lzzzz::lz4f::{WriteCompressor, BufReadDecompressor, Preferences, PreferencesBuilder, CLEVEL_HIGH};
use xx_bloomfilter::Bloom;
use flate2::read::GzDecoder;
use std::sync::atomic::{AtomicUsize, Ordering};
use glob::glob;
use dashmap::DashMap;
use thread_id;
use std::cell::UnsafeCell;
use std::io::Result;
use std::fmt::Arguments;
mod utils;
mod paf_output;
mod seq_output;
mod minimizers;
//mod ec_reads;
mod kmer_vec;
//mod poa;
mod read;
mod strobes;
mod closures;
mod kminmers;
mod ksyncmers;
mod kstrobes;
//mod pairwise;
//mod presimp;
use std::env;
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
pub struct RacyBloom(UnsafeCell<Bloom>); // intentionnally allowing data races as a tradeoff for bloom speed
unsafe impl Sync for RacyBloom {}
// follows https://sodocumentation.net/rust/topic/6018/unsafe-guidelines
// don't try this at home
impl RacyBloom {
    fn new(v: Bloom) -> RacyBloom {
        RacyBloom(UnsafeCell::new(v))
    }
    fn get(&self) -> &mut Bloom {
        // UnsafeCell::get() returns a raw pointer to the value it contains
        // Dereferencing a raw pointer is also "unsafe"
        unsafe {&mut *self.0.get()}
    }
}
type ThreadIdType = usize;
pub struct Params {
    l: usize,
    k: usize,
    density: f64,
    size_miniverse: u32, //delete
    average_lmer_count: f64,
    lmer_counts_min: u32,
    lmer_counts_max: u32,
    uhs: bool,
    lcp: bool,
    s: usize,
    has_lmer_counts: bool,
    use_bf: bool,
    use_hpc: bool,
    debug: bool,
    f: f64,
    wmin: usize,
    wmax: usize
}

fn debug_output_read_minimizers(seq_str: &String, read_minimizers: &Vec<String>, read_minimizers_pos: &Vec<u32>) {
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
        if current_minimizer.len() > 0 {
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
    let mut entry = hashmap_all.entry(thread_num).or_insert(HashMap::new());
    *entry = hashmap; // I believe hashmap is moved in this function as per https://stackoverflow.com/a/29490907 
}

pub fn thread_update_vec<U>(vec_all: &Arc<Mutex<HashMap<usize, Vec<U>>>>, vec: Vec<U>, thread_num: usize) {
    let mut vec_all = vec_all.lock().unwrap();
    let mut entry = vec_all.entry(thread_num).or_insert(Vec::new());
    *entry = vec;
}

fn get_reader(path: &PathBuf) -> Box<BufRead + Send> {
    let mut filetype = "unzip";
    let filename_str = path.to_str().unwrap();
    let file = match File::open(path) {
            Ok(file) => file,
            Err(error) => panic!("Error opening compressed file: {:?}.", error),
        };
    if filename_str.ends_with(".gz")  {filetype = "zip";}
    if filename_str.ends_with(".lz4") {filetype = "lz4";}
    let reader :Box<BufRead + Send> = match filetype { 
        "zip" => Box::new(BufReader::new(GzDecoder::new(file))), 
        "lz4" => Box::new(BufReadDecompressor::new(BufReader::new(file)).unwrap()),
        _ =>     Box::new(BufReader::new(file)), 
    }; 
    reader
}

fn read_first_n_reads(filename: &PathBuf, fasta_reads: bool, max_reads: usize) -> (usize, usize) {
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
}

fn autodetect_k_l_d(filename: &PathBuf, fasta_reads: bool) -> (usize, usize, f64) {
    println!("Parsing input sequences to estimate mean read length...");
    let (mean_length, max_length) = read_first_n_reads(&filename, fasta_reads, 100);
    println!("Detected mean read length of {} bp.",mean_length);
    // a bit crude, but let's try
    let d = 0.003;
    let coeff : f64 = 3.0/4.0;
    let slightly_below_readlen : f64 = (mean_length as f64);
    let k = (d * slightly_below_readlen) as usize;
    let l = 12;
    println!("Setting k = {} l = {} density = {}.", k, l, d);
    (k, l, d)
}

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
    /// Universal k-mer file (enables universal hitting sets (UHS))
    ///
    /// Universal k-mers need to be provided as a single file,
    /// one universal k-mer per line. The minimizers will be selected
    /// if they are a universal k-mer and if they satisfy the density
    /// bound.
    uhs: Option<String>,
    #[structopt(long)]
    /// Core substring file (enables locally consistent parsing (LCP))
    ///
    /// Core substrings need to be provided as a single file,
    /// one core substring per line. The minimizers will be selected
    /// if they are a core substring and if they satisfy the density
    /// bound.
    lcp: Option<String>,
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
    /// Enable Bloom filters
    ///
    /// Bloom filters can be used to reduce memory usage,
    /// but results in slightly less contiguous assemblies.
    #[structopt(long)]
    bf: bool,
    /// Homopolymer-compressed (HPC) input
    ///
    /// Both raw and homopolymer-compressed (HPC) reads can
    /// be provided as input. If the reads are not compressed,
    /// rust-mdbg manually performs HPC, but uses the raw sequences
    /// for transformation into base-space.
    #[structopt(long)]
    hpc: bool,
    /// l-mer counts (enables downweighting of frequent l-mers)
    ///
    /// Frequencies of l-mers in the reads (obtained using k-mer counters)
    /// can be provided in order to downweight frequently-occurring l-mers 
    /// and increase contiguity.
    #[structopt(parse(from_os_str), long)]
    lmer_counts: Option<PathBuf>,
    /// Minimum l-mer count threshold
    ///
    /// l-mers with frequencies below this threshold will be
    /// downweighted.
    #[structopt(long)]
    lmer_counts_min: Option<u32>,
    /// Maximum l-mer count threshold
    ///
    /// l-mers with frequencies above this threshold will be
    /// downweighted.
    #[structopt(long)]
    lmer_counts_max: Option<u32>,
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
    let mut uhs : bool = false;
    let mut lcp : bool = false;
    let mut filename = PathBuf::new();
    let mut ref_filename = PathBuf::new();
    let mut lmer_counts_filename = PathBuf::new();
    let mut uhs_filename = String::new();
    let mut lcp_filename = String::new();
    let mut output_prefix;
    let mut k : usize = 10;
    let mut l : usize = 12;
    let mut wmax : usize = 0;
    let mut wmin : usize = 0;
    let mut f : f64 = 0.0;
    let mut s : usize = 0;
    let mut density : f64 = 0.10;
    let mut reference : bool = false;
    let mut windowed : bool = false;
    let mut has_lmer_counts : bool = false;
    let mut lmer_counts_min : u32 = 2;
    let mut lmer_counts_max : u32 = 100000;
    let mut use_bf : bool = false;
    let mut use_strobe : bool = false;
    let mut use_hpc : bool = false;
    let mut threads : usize = 8;
    if !opt.reads.is_none() {filename = opt.reads.unwrap().clone();} 
    if !opt.reference.is_none() {ref_filename = opt.reference.unwrap().clone();} 
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
    if opt.k.is_none() && opt.l.is_none() && opt.density.is_none() {
        println!("Autodetecting values for k, l, and density.");
        let (ak, al, ad) = autodetect_k_l_d(&ref_filename, ref_fasta_reads);
        k = ak; l = al; density = ad;
    }
    else {
        if !opt.k.is_none() {k = opt.k.unwrap()} else {println!("Warning: Using default k value ({}).", k);} 
        if !opt.l.is_none() {l = opt.l.unwrap()} else {println!("Warning: Using default l value ({}).", l);}
        if !opt.density.is_none() {density = opt.density.unwrap()} else if opt.s.is_none() {println!("Warning: Using default density value ({}%).", density * 100.0);}
    }
    if !opt.threads.is_none() {threads = opt.threads.unwrap();} else {println!("Warning: Using default number of threads (8).");}
    if opt.bf {use_bf = true;}
    if opt.strobe {use_strobe = true;}
    if opt.hpc {use_hpc = true;}
    output_prefix = PathBuf::from(format!("hifimap-k{}-d{}-l{}", k, density, l));
    if !opt.lmer_counts.is_none() { 
        has_lmer_counts = true;
        lmer_counts_filename = opt.lmer_counts.unwrap().clone(); 
        if !opt.lmer_counts_min.is_none() {lmer_counts_min = opt.lmer_counts_min.unwrap();} else {println!("Warning: Using default l-mer minimum count ({}).", lmer_counts_min);}
        if !opt.lmer_counts_max.is_none() {lmer_counts_max = opt.lmer_counts_max.unwrap();} else {println!("Warning: Using default l-mer maximum count ({}).", lmer_counts_max);}
    } 
    if !opt.uhs.is_none() { 
        uhs = true;
        uhs_filename = opt.uhs.unwrap(); 
    }
    if !opt.lcp.is_none() { 
        lcp = true;
        lcp_filename = opt.lcp.unwrap(); 
    } 
    if !opt.s.is_none() { 
        s = opt.s.unwrap(); 
        println!("Syncmer parameter s: {}", s);
        if !opt.f.is_none() {f = opt.f.unwrap()} else {println!("Warning: Using default filter cutoff value ({}).", f);} 
        if !opt.wmin.is_none() {wmin = opt.wmin.unwrap()} else {println!("Warning: Using default wmin value ({}).", l/(l-s+1)+2);} 
        if !opt.wmax.is_none() {wmax = opt.wmax.unwrap()} else {println!("Warning: Using default wmax value ({}).", l/(l-s+1)+10);} 

    } 
    if !opt.prefix.is_none() {output_prefix = opt.prefix.unwrap();} else {println!("Warning: Using default output prefix ({}).", output_prefix.to_str().unwrap());}
    let debug = opt.debug;
    let size_miniverse = match revcomp_aware {
        false => 4f32.powf(l as f32) as u32,
        true => 4f32.powf(l as f32) as u32 / 2
    };
    let mut params = Params { 
        l,
        k,
        density,
        size_miniverse,
        average_lmer_count: 0.0,
        lmer_counts_min,
        lmer_counts_max,
        uhs,
        lcp,
        s,
        has_lmer_counts,
        use_bf,
        use_hpc,
        debug,
        f,
        wmin,
        wmax,
    };
    // init some useful objects
    let mut nb_minimizers_per_read : f64 = 0.0;
    let mut nb_reads : u64 = 0;
    // get file size for progress bar
    let metadata = fs::metadata(&filename).expect("Error opening input file.");
    let ref_metadata = fs::metadata(&ref_filename).expect("Error opening reference file.");
    let file_size = metadata.len();
    let mut pb = ProgressBar::on(stderr(),file_size);
    let mut queue_len = 200; // https://doc.rust-lang.org/std/sync/mpsc/fn.sync_channel.html
                             // also: controls how many reads objects are buffered during fasta/fastq
                             // parsing

    // dbg_nodes is a hash table containing (kmers -> (index,count))
    // it will keep only those with count > 1
   

    //let mut bloom : RacyBloom = RacyBloom::new(Bloom::new_with_rate(if use_bf {100_000_000} else {1}, 1e-7)); // a bf to avoid putting stuff into kmer_table too early
    let paf_path = format!("{}{}", output_prefix.to_str().unwrap(), ".paf");
    let mut paf_file = match File::create(&paf_path) {
        Err(why) => panic!("Couldn't create {}: {}", paf_path, why.description()),
        Ok(paf_file) => paf_file,
    };
    if params.s != 0 && opt.k.is_some() && use_strobe {
        closures::run_kstrobes(&filename, &ref_filename, &params, threads, queue_len, fasta_reads, ref_fasta_reads, &output_prefix);
    }
    else if params.s != 0 && opt.k.is_some() {
        closures::run_ksyncmers(&filename, &ref_filename, &params, threads, queue_len, fasta_reads, ref_fasta_reads, &output_prefix);
    }
    else if params.s != 0 && use_strobe {
        closures::run_randstrobes(&filename, &ref_filename, &params, threads, queue_len, fasta_reads, ref_fasta_reads, &output_prefix);
    }
    else {
        closures::run_kminmers(&filename, &ref_filename, &params, threads, queue_len, fasta_reads, ref_fasta_reads, &output_prefix);
    }
    let duration = start.elapsed();
    println!("Total execution time: {:?}", duration);
    println!("Maximum RSS: {:?}GB", (get_memory_rusage() as f32) / 1024.0 / 1024.0 / 1024.0);
}



