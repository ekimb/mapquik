// mapquik v0.1.0
// Copyright 2020-2021 Baris Ekim, Rayan Chikhi.
// Licensed under the MIT license (http://opensource.org/licenses/MIT).
// This file may not be copied, modified, or distributed except according to those terms.

#![allow(unused_variables)]
#![allow(non_upper_case_globals)]
//#![allow(warnings)]
#![feature(iter_advance_by)]
#[cfg(not(target_env = "msvc"))]
#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;
use crate::index::{Entry, Index, ReadOnlyIndex};
use crate::stats::Stats;
use std::fs::{self, File};
use std::io::{BufRead, BufReader};
use std::mem::{MaybeUninit};
use std::path::PathBuf;
use std::time::{Instant};
use flate2::read::GzDecoder;
use lzzzz::lz4f::{BufReadDecompressor};
use rust_seq2kminmers::{FH, KH};
use structopt::StructOpt;
mod chain;
mod closures;
mod index;
mod r#match;
mod mers;
mod stats;

pub type PseudoChainCoords = (bool, usize, usize, usize, usize, usize, usize);
pub type PseudoChainCoordsTuple<'a> = (usize, PseudoChainCoords);
pub struct Params {
    k: usize,
    l: usize,
    density: FH,
    use_hpc: bool,
    use_simd: bool,
    use_pfx: bool,
    debug: bool,
    a: bool,
    c: usize, // minimum chain length
    s: usize, // minimum match score (# of matching seeds)
    g: usize, // maximum gap difference
    b: usize, // buffer increase
    q: usize, // queue length
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

fn get_reader(path: &PathBuf) -> (Box<dyn BufRead + Send>, bool) {
    let mut filetype = "unzip";
    let filename_str = path.to_str().unwrap();
    let file = match File::open(path) {
            Ok(file) => file,
            Err(error) => panic!("Error opening compressed file: {:?}.", error),
        };
    if filename_str.ends_with(".gz")  {filetype = "zip";}
    if filename_str.ends_with(".lz4") {filetype = "lz4";}
    let reader :(Box<dyn BufRead + Send>,bool) = match filetype { 
        "zip" => (Box::new(BufReader::new(GzDecoder::new(file))),true), 
        "lz4" => (Box::new(BufReadDecompressor::new(BufReader::new(file)).unwrap()),true),
        _ =>     (Box::new(BufReader::new(file)),false), 
    }; 
    reader
}

#[derive(Debug, StructOpt)]
#[structopt(name = "mapquik")]
/// Original implementation of mapquik, a fast HiFi read mapper.
struct Opt {
    /// Activate debug mode
    ///
    #[structopt(long)]
    debug: bool,
    /// Input file (raw or gzip-/lz4-compressed FASTX)
    ///
    /// Input file can be FASTA/FASTQ, as well as gzip-compressed (.gz) or
    /// lz4-compressed (.lz4). Lowercase bases are currently not supported;
    /// see documentation for formatting.
    #[structopt(parse(from_os_str))]
    reads: Option<PathBuf>,
    /// Output prefix for PAF file
    /// 
    #[structopt(parse(from_os_str), short, long)]
    prefix: Option<PathBuf>,
    /// k-min-mer length
    ///
    /// The length of each k-min-mer. If
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
    density: Option<FH>,
    /// Minimum chain length
    ///
    /// Only outputs chains of length >= c.
    #[structopt(short, long)]
    chain: Option<usize>,
    /// Minimum number of matching seeds
    ///
    /// Only outputs chains of >= s matching seeds.
    #[structopt(short, long)]
    seed: Option<usize>,
    /// Maximum nucleotide gap length difference
    ///
    /// Allows chaining of hits with a gap difference of < g.
    #[structopt(short, long)]
    gap_diff: Option<usize>,
    /// Reference genome input
    ///
    /// Reference to be indexed and mapped to. 
    /// Allows multi-line FASTA and
    /// doesn't filter any kminmers.
    #[structopt(parse(from_os_str), long)]
    reference: Option<PathBuf>,
    /// Number of threads
    /// 
    #[structopt(long)]
    threads: Option<usize>,
    /// Enable base-level alignment (WIP)
    /// 
    //#[structopt(short, long)]
   // align: bool,
    /// Enable low-memory reference FASTA parsing
    /// 
    #[structopt(long)]
    low_memory: bool,
    /// Deactivate SIMD (AVX2,AVX512) functions (for old processors)
    #[structopt(long)]
    nosimd: bool,
    /// Deactivate HomoPolymer Compression
    #[structopt(long)]
    nohpc: bool,
    /// Use parallelfastx (faster uncompressed reads parsing)
    #[structopt(long)]
    parallelfastx: bool,
    /// buffer size multiplier
    #[structopt(short, long)]
    b: Option<usize>,
    /// queue length
    #[structopt(short, long)]
    q: Option<usize>,

}

fn main() {
    let start = Instant::now();
    let opt = Opt::from_args();      
    let mut filename = PathBuf::new();
    let mut ref_filename = PathBuf::new();
    let mut output_prefix;
    let mut k : usize = 5;
    let mut l : usize = 31;
    let mut c = 4;
    let mut s = 11;
    let mut g = 2000;
    let mut b = 1;
    let mut q = 200;
    let low_memory = opt.low_memory;
    let a = false;
    let mut density : FH = 0.01;
    let reference : bool = false;
    let mut use_hpc : bool = true; 
    let mut use_simd : bool = true; 
    let mut use_pfx : bool = false; 
    let mut threads : usize = 8;
    if opt.reads.is_some() {filename = opt.reads.unwrap();} 
    if opt.reference.is_some() {ref_filename = opt.reference.unwrap();} 
    if filename.as_os_str().is_empty() {panic!("Please specify an input file.");}
    if ref_filename.as_os_str().is_empty() {panic!("Please specify a reference file.");}
    let mut reads_are_fasta : bool = false;
    let mut ref_is_fasta    : bool = false;
    let filename_str = filename.to_str().unwrap();
    if filename_str.contains(".fasta.") || filename_str.ends_with(".fna") || filename_str.contains(".fna.") || filename_str.contains(".fa.") || filename_str.ends_with(".fa") || filename_str.ends_with(".fasta") { // not so robust but will have to do for now
        reads_are_fasta = true;
        println!("Input file: {}", filename_str);
        println!("Format: FASTA");
    }
    let ref_filename_str = ref_filename.to_str().unwrap();
    if ref_filename_str.contains(".fasta.") || ref_filename_str.ends_with(".fna") || ref_filename_str.contains(".fna.") || ref_filename_str.contains(".fa.") || ref_filename_str.ends_with(".fa") || ref_filename_str.ends_with(".fasta") { // not so robust but will have to do for now
        ref_is_fasta = true;
        println!("Reference file: {}", ref_filename_str);
        println!("Format: FASTA");
    }
    if opt.k.is_some() {k = opt.k.unwrap()} else {println!("Warning: Using default k value ({}).", k);} 
    if opt.l.is_some() {l = opt.l.unwrap()} else {println!("Warning: Using default l value ({}).", l);}
    if opt.b.is_some() {b = opt.b.unwrap()} else {println!("Warning: Using default buffer size ({}X).", b);}
    if opt.q.is_some() {q = opt.q.unwrap()} else {println!("Warning: Using default queue length ({}).", q);}
    if opt.density.is_some() {density = opt.density.unwrap()} else {println!("Warning: Using default density value ({}%).", density * 100.0);}
    if opt.threads.is_some() {threads = opt.threads.unwrap();} else {println!("Warning: Using default number of threads (8).");}
    if opt.chain.is_some() {c = opt.chain.unwrap()} else {println!("Warning: Using default minimum chain length ({}).", c);}
    if opt.seed.is_some() {s = opt.seed.unwrap()} else {println!("Warning: Using default minimum number of matching seeds ({}).", s);}
    if opt.gap_diff.is_some() {g = opt.gap_diff.unwrap()} else {println!("Warning: Using default maximum seed gap difference ({}).", g);}
    output_prefix = PathBuf::from(format!("mapquik-k{}-d{}-l{}", k, density, l));
    if opt.prefix.is_some() {output_prefix = opt.prefix.unwrap();} else {println!("Warning: Using default output prefix ({}).", output_prefix.to_str().unwrap());}
    let debug = opt.debug;
    if opt.nohpc  { use_hpc = false; }
    if opt.nosimd { use_simd = false; }
    if opt.parallelfastx { use_pfx = true; }
    if ! std::is_x86_feature_detected!("avx512f") { 
        println!("Warning: No AVX-512 CPU found, falling back to scalar implementation");
        use_simd = false; 
    }
    if use_hpc {
        if use_simd {
            println!("Using HPC ntHash, with SIMD");
        }
            else {
            println!("Using HPC ntHash, scalar");
        }
    } else {
        if use_simd {
            println!("Using regular ntHash (not HPC), with SIMD");
        }
            else {
            println!("Using regular ntHash (not HPC), scalar");
        }
    }
    let params = Params { 
        k,
        l,
        density,
        use_hpc,
        use_simd,
        use_pfx,
        debug,
        a,
        c,
        s,
        g,
        b,
        q,
    };
    // init some useful objects
    // get file size for progress bar
    let _metadata = fs::metadata(&filename).expect("Error opening input file.");
    let _ref_metadata = fs::metadata(&ref_filename).expect("Error opening reference file.");
    let ref_threads = threads;
    let mut ref_queue_len = threads;
    if low_memory {ref_queue_len = 1;}
    let queue_len = params.q; // https://doc.rust-lang.org/std/sync/mpsc/fn.sync_channel.html
                             // also: controls how many reads objects are buffered during fasta/fastq
                             // parsing
    Stats::init(threads, output_prefix.to_str().unwrap());
    closures::run_mers(&filename, &ref_filename, &params, ref_threads, threads, ref_queue_len, queue_len, reads_are_fasta, ref_is_fasta, &output_prefix);
    //println!("current time after exiting closures {:?}",Utc::now());
    let duration = start.elapsed();
    println!("Total execution time: {:?}", duration);
    println!("Maximum RSS: {:?}GB", (get_memory_rusage() as f32) / 1024.0 / 1024.0 / 1024.0);
}



