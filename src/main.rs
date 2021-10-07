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
use std::fs;
use structopt::StructOpt;
use std::sync::{Arc, Mutex};
use std::path::PathBuf;
use std::time::{Instant};
use std::mem::{MaybeUninit};
use seq_io::BaseRecord;
use lzzzz::lz4f::{WriteCompressor, BufReadDecompressor, Preferences};
use flate2::read::GzDecoder;
use dashmap::DashMap;
use std::io::Result;
use core::hash::{Hash, Hasher};
use std::collections::hash_map::DefaultHasher;
mod utils;
mod closures;
mod mers;
mod graph;
mod kmer_vec;

pub struct Params {
    k: usize,
    l: usize,
    d: f64,
    s: usize,
    wmax: usize,
    wmin: usize,
    f: f64,
    fmax: usize,
    fmin: usize,
    hpc: bool,
    strobe: bool,
    z: f64,
    nam: bool,
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

#[derive(Debug, StructOpt)]
#[structopt(name = "hifimap")]
/// Original implementation of hifimap, a fast HiFi read mapper.
struct Opt {
    /// Activate debug mode
    #[structopt(long)]
    debug: bool,
    /// Input file (raw or gzip-/lz4-compressed FASTX)
    #[structopt(parse(from_os_str))]
    reads: Option<PathBuf>,
    /// Output prefix for PAF file
    #[structopt(parse(from_os_str), short, long)]
    prefix: Option<PathBuf>,
    /// k-min-mer length
    #[structopt(short, long)]
    k: Option<usize>,
    /// l-mer (minimizer) length
    #[structopt(short, long)]
    l: Option<usize>,
    /// Density threshold for density-based selection scheme
    #[structopt(short, long)]
    density: Option<f64>,
    /// Syncmer length for syncmer-based selection scheme
    #[structopt(short, long)]
    s: Option<usize>,
    /// Filter cutoff value
    #[structopt(short, long)]
    f: Option<f64>,
    /// Minimum l-mer count
    #[structopt(long)]
    fmin: Option<usize>,
    /// Maximum l-mer count
    #[structopt(long)]
    fmax: Option<usize>,
    /// Syncmer-based window minimum value
    #[structopt(long)]
    wmin: Option<usize>,
    /// Syncmer-based window maximum value
    #[structopt(long)]
    wmax: Option<usize>,
    /// Reference genome input
    ///
    /// Reference to be indexed and mapped to
    /// Allows multi-line FASTA and
    /// doesn't filter any kminmers
    #[structopt(parse(from_os_str), long)]
    reference: Option<PathBuf>,
    /// Enable strobemers
    #[structopt(long)]
    strobe: bool,
    /// Homopolymer-compressed (HPC) input
    #[structopt(long)]
    hpc: bool,
    /// Use z scores
    #[structopt(short, long)]
    z: Option<f64>,
    /// Number of threads
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
    let mut d : f64 = 0.10;
    let mut s : usize = 0;
    let mut wmax : usize = 0;
    let mut wmin : usize = 0;
    let mut f : f64 = 0.0;
    let mut fmax : usize = usize::max_value();
    let mut fmin : usize = 2; 
    let mut z = 0.0;
    let mut nam = true;
    let mut strobe : bool = false;
    let mut hpc : bool = false;
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
    if opt.density.is_some() {d = opt.density.unwrap()} else if opt.s.is_none() {println!("Warning: Using default density value ({}%).", d * 100.0);}
    output_prefix = PathBuf::from(format!("hifimap-k{}-d{}-l{}", k, d, l));
    if opt.threads.is_some() {threads = opt.threads.unwrap();} else {println!("Warning: Using default number of threads (8).");}
    if opt.strobe {strobe = true;}
    if opt.hpc {hpc = true;}
    if opt.z.is_some() {z = opt.z.unwrap(); nam = false;}
    if opt.f.is_some() {f = opt.f.unwrap()} else {println!("Warning: Not using relative frequency filtering.");} 
    if opt.s.is_some() { 
        s = opt.s.unwrap(); 
        println!("Syncmer parameter s: {}", s);
        if opt.wmin.is_some() {wmin = opt.wmin.unwrap()} else {println!("Warning: Using default wmin value ({}).", l/(l-s+1)+2);} 
        if opt.wmax.is_some() {wmax = opt.wmax.unwrap()} else {println!("Warning: Using default wmax value ({}).", l/(l-s+1)+10);} 
    } 
    if opt.fmax.is_some() {fmax = opt.fmax.unwrap();}
    if opt.fmin.is_some() {fmin = opt.fmin.unwrap();}
    if opt.prefix.is_some() {output_prefix = opt.prefix.unwrap();} else {println!("Warning: Using default output prefix ({}).", output_prefix.to_str().unwrap());}

    let params = Params {
        k,
        l,
        d,
        s,
        wmax,
        wmin,
        f,
        fmax,
        fmin,
        hpc,
        strobe,
        z,
        nam,
    };
    let metadata = fs::metadata(&filename).expect("Error opening input file.");
    let ref_metadata = fs::metadata(&ref_filename).expect("Error opening reference file.");
    let file_size = metadata.len();
    let queue_len = threads;
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



