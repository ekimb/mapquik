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
use closure::closure;
use std::iter::FromIterator;
use crate::kmer_vec::get;
use crate::read::Read;
use crate::read::ReadSync;
use std::fs::{File,remove_file};
use std::collections::HashSet;
extern crate array_tool;
use std::fs;
use crossbeam_utils::{thread};
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
mod bit;
//mod pairwise;
//mod presimp;
use std::env;

const revcomp_aware : bool = true; // shouldn't be set to false except for strand-directed data or for debugging
type Kmer = kmer_vec::KmerVec;
type Overlap = kmer_vec::KmerVec;
#[derive(Clone, Debug)] // seems necessary to move out of the Arc into dbg_nodes_view
pub struct HashEntry {origin: String, seq: String, seqlen: u32, shift: (usize, usize), seq_rev: bool}
#[derive(Clone, Debug)] // seems necessary to move out of the Arc into dbg_nodes_view
pub struct QueryMatch {kminmer: Kmer, target: Vec<HashEntry>, query: HashEntry}
struct SeqFileType(WriteCompressor<File>);
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
    syncmer: usize,
    has_lmer_counts: bool,
    use_bf: bool,
    use_hpc: bool,
    debug: bool,
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
    /// Syncmer-based selection scheme
    #[structopt(short, long)]
    syncmer: Option<usize>,
    /// Reference genome input
    ///
    /// Reference to be indexed and mapped to. 
    /// Allows multi-line FASTA and
    /// doesn't filter any kminmers.
    #[structopt(parse(from_os_str), long)]
    reference: Option<PathBuf>,
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
    let mut w : usize = 0;
    let mut syncmer : usize = 0;
    let mut density : f64 = 0.10;
    let mut reference : bool = false;
    let mut windowed : bool = false;
    let mut has_lmer_counts : bool = false;
    let mut lmer_counts_min : u32 = 2;
    let mut lmer_counts_max : u32 = 100000;
    let mut use_bf : bool = false;
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
        if !opt.density.is_none() {density = opt.density.unwrap()} else if opt.syncmer.is_none() {println!("Warning: Using default density value ({}%).", density * 100.0);}
    }
    if !opt.threads.is_none() {threads = opt.threads.unwrap();} else {println!("Warning: Using default number of threads (8).");}
    if opt.bf {use_bf = true;}
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
    if !opt.syncmer.is_none() { 
        syncmer = opt.syncmer.unwrap(); 
        println!("Syncmer parameter s: {}", syncmer);
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
        syncmer,
        has_lmer_counts,
        use_bf,
        use_hpc,
        debug,
    };
    // init some useful objects
    let mut nb_minimizers_per_read : f64 = 0.0;
    let mut nb_reads : u64 = 0;
    // get file size for progress bar
    let metadata = fs::metadata(&filename).expect("Error opening input file.");
    let ref_metadata = fs::metadata(&ref_filename).expect("Error opening reference file.");
    let file_size = metadata.len();
    let mut pb = ProgressBar::on(stderr(),file_size);
    let mut lmer_counts : HashMap<String, u32> = HashMap::new();

    if has_lmer_counts {
        let lmer_counts_file = match File::open(lmer_counts_filename) {
            Err(why) => panic!("Couldn't load l-mer counts: {}.", why.description()),
            Ok(lmer_counts_file) => lmer_counts_file,
        }; 
        let mut br = BufReader::new(lmer_counts_file);
        loop {
            let mut line = String::new();
            let new_line = |line: &mut String, br: &mut BufReader<File>| {line.clear(); br.read_line(line).ok();};
            if let Err(e) = br.read_line(&mut line) {break;}
            if line.len() == 0                      {break;}
            let trimmed = line.trim().to_string();   
            let vec : Vec<String> = trimmed.split(" ").map(String::from).collect();
            let lmer = vec[0].to_string();
            let lmer_rev = utils::revcomp(&lmer);
            let lmer = if lmer > lmer_rev {lmer} else {lmer_rev}; //don't trust the kmer counter to normalize like we do
            let count = vec[1].parse::<u32>().unwrap();
            lmer_counts.insert(lmer, count);               
            new_line(&mut line, &mut br);
        }
    }
    let mut minimizer_to_int : HashMap<String,u64> = HashMap::new();
    let mut int_to_minimizer : HashMap<u64,String> = HashMap::new();
    // only need to initialize the minimizer_to_int / int_to_minimizer array if we do POA or use robust minimizers
    // they can be costly for k=14
    if has_lmer_counts {
        let res = minimizers::minimizers_preparation(&mut params, &filename, file_size, &lmer_counts);
        minimizer_to_int = res.0;
        int_to_minimizer = res.1;
    }
    let (mean_length, max_length) = read_first_n_reads(&filename, fasta_reads, 10);
    let mut queue_len = 200; // https://doc.rust-lang.org/std/sync/mpsc/fn.sync_channel.html
                             // also: controls how many reads objects are buffered during fasta/fastq
                             // parsing

    let mut uhs_bloom : RacyBloom = RacyBloom::new(Bloom::new(if use_bf {500_000_000} else {1}, 1_000_000_000_000_000));
    let mut lcp_bloom : RacyBloom = RacyBloom::new(Bloom::new(if use_bf {500_000_000} else {1}, 1_000_000_000_000_000));

    if params.uhs {uhs_bloom = minimizers::uhs_preparation(&mut params, &uhs_filename)}
    if params.lcp {lcp_bloom = minimizers::lcp_preparation(&mut params, &lcp_filename)}

    // dbg_nodes is a hash table containing (kmers -> (index,count))
    // it will keep only those with count > 1
    let mut kmer_table     : Arc<DashMap<Kmer, Vec<HashEntry>>> = Arc::new(DashMap::new()); // it's a Counter
    let mut ref_abundance_table     : Arc<DashMap<Kmer, usize>> = Arc::new(DashMap::new()); // it's a Counter
    let mut query_abundance_table     : Arc<DashMap<Kmer, usize>> = Arc::new(DashMap::new()); // it's a Counter

    let mut query_matches     : Arc<DashMap<String, Vec<QueryMatch>>> = Arc::new(DashMap::new()); // it's a Counter

    //let mut bloom : RacyBloom = RacyBloom::new(Bloom::new_with_rate(if use_bf {100_000_000} else {1}, 1e-7)); // a bf to avoid putting stuff into kmer_table too early
    let mut bloom         : RacyBloom = RacyBloom::new(Bloom::new(if use_bf {500_000_000} else {1}, 1_000_000_000_000_000)); // same bf but making sure we use only 1 hash function for speed
    static node_index     : AtomicUsize = AtomicUsize::new(0); // associates a unique integer to each k-min-mer
    let mut kmer_seqs     : HashMap<Kmer, String> = HashMap::new(); // associate a k-min-mer to an arbitrary sequence from the reads
    let mut kmer_seqs_lens: HashMap<Kmer, Vec<u32>> = HashMap::new(); // associate a k-min-mer to the lengths of all its sequences from the reads
    let mut kmer_origin   : HashMap<Kmer, String> = HashMap::new(); // remember where in the read/refgenome the kmer comes from

    // delete all previous sequence files
    for path in glob(&format!("{}*.sequences", output_prefix.to_str().unwrap()).to_string()).expect("Failed to read glob pattern.") {
        let path = path.unwrap();
        let path = path.to_str().unwrap(); // rust really requires me to split the let statement in two..
        println!("Removing old sequences file: {}.", &path);
        fs::remove_file(path);
    }
    let mut seq_write = |file: &mut SeqFileType, s| {write!(file, "{}", s);};
    let mut sequences_files : Arc<DashMap<ThreadIdType, SeqFileType>> = Arc::new(DashMap::new());
    let create_sequences_file = |thread_id: ThreadIdType| -> SeqFileType {
        let seq_path = PathBuf::from(format!("{}.{}.sequences", output_prefix.to_str().unwrap(), thread_id));
        let mut file = match File::create(&seq_path) {
            Err(why) => panic!("Couldn't create file: {}.", why.description()),
            Ok(file) => file,
        };
        //let mut sequences_file = BufWriter::new(file);
        let mut sequences_file = SeqFileType::new(file); // regular lz4f
        //let mut sequences_file = WriteCompressor::new(&mut file, PreferencesBuilder::new().compression_level(CLEVEL_HIGH).build()).unwrap();  // too slow
        seq_write(&mut sequences_file, format!("# k = {}\n",k));
        seq_write(&mut sequences_file, format!("# l = {}\n",l));
        seq_write(&mut sequences_file, "# Structure of remaining of the file:\n".to_string());
        seq_write(&mut sequences_file, "# [node name]\t[list of minimizers]\t[sequence of node]\t[abundance]\t[origin]\t[shift]\n".to_string());
        sequences_file
    };

    let mut add_ref_kminmer = |ref_kminmer: Kmer, seq: Option<&str>, seq_reversed: &bool, origin: &str, shift: &(usize, usize), sequences_file: &mut SeqFileType, thread_id: usize, read_seq: Option<&str>, read_offsets: Option<(usize, usize, usize)>|
    {
        //println!("Seq: {:?}, read_seq: {:?}", seq, read_seq);
        let seq = if *seq_reversed {utils::revcomp(&seq.unwrap())} else {seq.unwrap().to_string()};
        let seqlen = read_seq.unwrap().len() as u32;
        let seq_line = format!("{}\t{}\t{}\t{}\t{:?}", ref_kminmer.print_as_string(), seq, "*", origin, shift);
        seq_write(sequences_file, format!("{}\n", seq_line));
        let mut contains_key = kmer_table.contains_key(&ref_kminmer);
        if contains_key {
            let mut entry_vec = kmer_table.get_mut(&ref_kminmer).unwrap();
            entry_vec.push(HashEntry{origin: origin.to_string(), seq: seq, seqlen: seqlen, shift: *shift, seq_rev: *seq_reversed}); 
        }
        else {
            let mut entry_vec = Vec::<HashEntry>::new();
            entry_vec.push(HashEntry{origin: origin.to_string(), seq: seq, seqlen: seqlen, shift: *shift, seq_rev: *seq_reversed});
            kmer_table.insert(ref_kminmer.clone(), entry_vec);
        }
        let mut contains_abund = ref_abundance_table.contains_key(&ref_kminmer);
        if contains_abund {*ref_abundance_table.get_mut(&ref_kminmer).unwrap() += 1;}
        else {ref_abundance_table.insert(ref_kminmer.clone(), 1);}
    };
    let ref_process_read_aux_sync = |ref_str: &[u8], ref_id: &str| -> Option<(Vec<(Kmer, String, bool, String, (usize, usize))>, Option<Read>, Option<ReadSync>)> {
        let thread_id :usize =  thread_id::get();
        let mut output : Vec<(Kmer, String, bool, String, (usize, usize))> = Vec::new();
        let seq = std::str::from_utf8(ref_str).unwrap(); // see https://docs.rs/seq_io/0.4.0-alpha.0/seq_io/fasta/index.html
            // those two next lines do a string copy of the read sequence, because in the case of a
            // reference, we'll need to remove newlines. also, that sequence will be moved to the
            // Read object. could be optimized later
        let seq_for_ref = seq.replace("\n", "").replace("\r", "");
        let seq = seq_for_ref; 
        let mut ref_obj = ReadSync::extract(&ref_id, seq, &params, &minimizer_to_int, &int_to_minimizer, &uhs_bloom, &lcp_bloom);
        if sequences_files.get(&thread_id).is_none() {
            sequences_files.insert(thread_id, create_sequences_file(thread_id));
        }
        let mut sequences_file = sequences_files.get_mut(&thread_id).unwrap();
        if ref_obj.transformed.len() > k {
            for i in 0..(ref_obj.transformed.len() - k + 1) {
                let mut ref_kminmer : Kmer = Kmer::make_from(&ref_obj.transformed[i..i+k]);
                // uncomment the line below to track where the kmer is coming from (but not needed in production)
                //let origin = "*".to_string(); 
                let mut seq_reversed = false;
                if revcomp_aware { 
                    let (ref_kminmer_norm, reversed) = ref_kminmer.normalize(); 
                    ref_kminmer = ref_kminmer_norm;
                    seq_reversed = reversed;
                } 
                
                let minimizers_pos = &ref_obj.minimizers_pos;
                let position_of_first_minimizer = match seq_reversed {
                    true => minimizers_pos[i+k-1].1 + params.l - 1,
                    false => minimizers_pos[i].0
                };
                let position_of_last_minimizer = match seq_reversed {
                    true => minimizers_pos[i].0 + params.l - 1,
                    false => minimizers_pos[i+k-1].1
                };
                let mut shift = (minimizers_pos[i].0, minimizers_pos[i+k-1].1 + params.l - 1);
                let ref_offsets = (minimizers_pos[i].0 as usize, (minimizers_pos[i+k-1].1 as usize + l - 1), (minimizers_pos[i+k-1].1 + l - minimizers_pos[i].0));
                //println!("Shift: {:?}", shift);
                let mut span_seq = &ref_obj.seq[shift.0..shift.1 + 1];
                //println!("Span_seq: {:?}, len: {:?}", span_seq, span_seq.len());
                //start of kminmer : pos[i]
                //end of kminmer : pos [i+k-1]+l
                add_ref_kminmer(ref_kminmer, Some(span_seq), &seq_reversed, &ref_id, &shift, &mut sequences_file, thread_id, Some(&ref_obj.seq), Some(ref_offsets));
            }
        }
        Some((output, None, Some(ref_obj)))
    };
    // worker thread
    let ref_process_read_aux = |ref_str: &[u8], ref_id: &str| -> Option<(Vec<(Kmer, String, bool, String, (usize, usize))>, Option<Read>, Option<ReadSync>)> {
        let thread_id :usize =  thread_id::get();
        let mut output : Vec<(Kmer, String, bool, String, (usize, usize))> = Vec::new();
        let seq = std::str::from_utf8(ref_str).unwrap(); // see https://docs.rs/seq_io/0.4.0-alpha.0/seq_io/fasta/index.html
            // those two next lines do a string copy of the read sequence, because in the case of a
            // reference, we'll need to remove newlines. also, that sequence will be moved to the
            // Read object. could be optimized later
            let seq_for_ref = seq.replace("\n", "").replace("\r", "");
            let seq = seq_for_ref; 
        let mut ref_obj = Read::extract(&ref_id, seq, &params, &minimizer_to_int, &int_to_minimizer, &uhs_bloom, &lcp_bloom);
        if sequences_files.get(&thread_id).is_none() {
            sequences_files.insert(thread_id, create_sequences_file(thread_id));
        }
        let mut sequences_file = sequences_files.get_mut(&thread_id).unwrap();
        if ref_obj.transformed.len() > k {
            for i in 0..(ref_obj.transformed.len() - k + 1) {
                let mut ref_kminmer : Kmer = Kmer::make_from(&ref_obj.transformed[i..i+k]);
                // uncomment the line below to track where the kmer is coming from (but not needed in production)
                //let origin = "*".to_string(); 
                let mut seq_reversed = false;
                if revcomp_aware { 
                    let (ref_kminmer_norm, reversed) = ref_kminmer.normalize(); 
                    ref_kminmer = ref_kminmer_norm;
                    seq_reversed = reversed;
                } 
                
                let minimizers_pos = &ref_obj.minimizers_pos;
                let position_of_first_minimizer = match seq_reversed {
                    true => minimizers_pos[i+k-1] + params.l - 1,
                    false => minimizers_pos[i]
                };
                let position_of_last_minimizer = match seq_reversed {
                    true => minimizers_pos[i] + params.l - 1,
                    false => minimizers_pos[i+k-1]
                };
                let mut shift = (minimizers_pos[i], minimizers_pos[i+k-1] + params.l - 1);
                let ref_offsets = (ref_obj.minimizers_pos[i] as usize, (ref_obj.minimizers_pos[i+k-1] as usize + l - 1), (ref_obj.minimizers_pos[i+k-1] + l - ref_obj.minimizers_pos[i]));
                //println!("Shift: {:?}", shift);
                let mut span_seq = &ref_obj.seq[shift.0..shift.1 + 1];
            // println!("Span_seq: {:?}, len: {:?}", span_seq, span_seq.len());
                //start of kminmer : pos[i]
                //end of kminmer : pos [i+k-1]+l
                add_ref_kminmer(ref_kminmer, Some(span_seq), &seq_reversed, &ref_id, &shift, &mut sequences_file, thread_id, Some(&ref_obj.seq), Some(ref_offsets));
            }
        }
        Some((output, Some(ref_obj), None))
    };
    let ref_process_read_fasta = |record: seq_io::fasta::RefRecord, found: &mut Option<(Vec<(Kmer, String, bool, String, (usize, usize))>, Option<Read>, Option<ReadSync>)>| {
        let ref_str = record.seq(); 
        let ref_id = record.id().unwrap().to_string();
        if params.syncmer != 0 {*found = ref_process_read_aux_sync(&ref_str, &ref_id);}
        else {*found = ref_process_read_aux(&ref_str, &ref_id);}
    };
    let ref_process_read_fastq = |record: seq_io::fastq::RefRecord, found: &mut Option<(Vec<(Kmer, String, bool, String, (usize, usize))>,Option<Read>, Option<ReadSync>)>| {
        let ref_str = record.seq(); 
        let ref_id = record.id().unwrap().to_string();
        *found = ref_process_read_aux(&ref_str, &ref_id);
    };

    // parallel fasta parsing, with a main thread that writes to disk and populates hash tables
    let mut ref_main_thread = |found: &Option<(Vec<(Kmer, String, bool, String, (usize, usize))>, Option<Read>, Option<ReadSync>)>| { // runs in main thread
        nb_reads += 1;
        //println!("Received read in main thread, nb kmers: {}", vec.len());
        let debug_only_display_read_and_minimizers = false;
        if debug_only_display_read_and_minimizers {
            // debug: just displays the read id and the list of minimizers
            let (vec, ref_obj, ref_obj_sync) = found.as_ref().unwrap();
            if ref_obj.is_some(){println!("{} {}", ref_obj.as_ref().unwrap().id.to_string(), ref_obj.as_ref().unwrap().transformed.to_vec().iter().join(" "));}
            else {println!("{} {}", ref_obj_sync.as_ref().unwrap().id.to_string(), ref_obj_sync.as_ref().unwrap().transformed.to_vec().iter().join(" "));}
        }
        else {
            let (vec, ref_obj, ref_obj_sync) = found.as_ref().unwrap();
        }
        // Some(value) will stop the reader, and the value will be returned.
        // In the case of never stopping, we need to give the compiler a hint about the
        // type parameter, thus the special 'turbofish' notation is needed.
        None::<()>
    };

    let buf = get_reader(&ref_filename);
    println!("Parsing reference sequence...");
    if ref_fasta_reads {
        let reader = seq_io::fasta::Reader::new(buf);
        read_process_fasta_records(reader, threads as u32, queue_len, ref_process_read_fasta, |record, found| {ref_main_thread(found)});
    }
    else {
        let reader = seq_io::fastq::Reader::new(buf);
        read_process_fastq_records(reader, threads as u32, queue_len, ref_process_read_fastq, |record, found| {ref_main_thread(found)});
    }

    println!("Indexed {} reference k-min-mers.", kmer_table.len());

    for mut sequences_file in sequences_files.iter_mut() {sequences_file.flush().unwrap();}

    let mut query_kminmer_lookup =|inp_kminmer: Kmer, seq: Option<&str>, seq_reversed: &bool, origin: &str, shift: &(usize, usize), sequences_file: &mut SeqFileType, thread_id: usize, read_seq: Option<&str>, read_offsets: Option<(usize, usize, usize)>| -> Option<QueryMatch>
    {
        let mut contains_abund = query_abundance_table.contains_key(&inp_kminmer);
        if contains_abund {*query_abundance_table.get_mut(&inp_kminmer).unwrap() += 1;}
        else {query_abundance_table.insert(inp_kminmer.clone(), 1);}
        let seq = if *seq_reversed {utils::revcomp(&seq.unwrap())} else {seq.unwrap().to_string()};
        let seqlen = read_seq.unwrap().len() as u32;
        let mut contains_key = kmer_table.contains_key(&inp_kminmer);
        if contains_key {
            let entry_vec = kmer_table.get_mut(&inp_kminmer).unwrap();
            let seq_line = format!("{}\t{}\t{}\t{}\t{:?}", inp_kminmer.print_as_string(), seq, "*", origin, shift);
            seq_write(sequences_file, format!("{}\n", seq_line));
            return Some(QueryMatch{kminmer: inp_kminmer, target: entry_vec.to_vec(), query: HashEntry{origin: origin.to_string(), seq: seq, seqlen: seqlen, shift: *shift, seq_rev: *seq_reversed}});
        }
        else {return None};
    };
        // worker thread
    let query_process_read_aux = |seq_str: &[u8], seq_id: &str| -> Option<(Vec<QueryMatch>, Option<Read>, Option<ReadSync>)> {
        let thread_id :usize =  thread_id::get();
        let mut output : Vec<(Kmer, String, bool, String, (usize, usize))> = Vec::new();
        let seq = std::str::from_utf8(seq_str).unwrap(); // see https://docs.rs/seq_io/0.4.0-alpha.0/seq_io/fasta/index.html
        // those two next lines do a string copy of the read sequence, because in the case of a
        // reference, we'll need to remove newlines. also, that sequence will be moved to the
        // Read object. could be optimized later
        let mut read_obj = Read::extract(&seq_id, seq.to_string(), &params, &minimizer_to_int, &int_to_minimizer, &uhs_bloom, &lcp_bloom);
        let mut lookups = Vec::<QueryMatch>::new();
        if sequences_files.get(&thread_id).is_none() {
            sequences_files.insert(thread_id, create_sequences_file(thread_id));
        }
        let mut sequences_file = sequences_files.get_mut(&thread_id).unwrap();
        if read_obj.transformed.len() > k {
            for i in 0..(read_obj.transformed.len() - k + 1) {
                let mut kminmer : Kmer = Kmer::make_from(&read_obj.transformed[i..i+k]);
                let mut seq_reversed = false;
                if revcomp_aware { 
                    let (kminmer_norm, reversed) = kminmer.normalize(); 
                    kminmer = kminmer_norm;
                    seq_reversed = reversed;
                } 
                let minimizers_pos = &read_obj.minimizers_pos;
                let position_of_first_minimizer = match seq_reversed {
                    true => minimizers_pos[i+k-1] + params.l - 1,
                    false => minimizers_pos[i]
                };
                let position_of_last_minimizer = match seq_reversed {
                    true => minimizers_pos[i] + params.l - 1,
                    false => minimizers_pos[i+k-1]
                };
                let mut shift = (minimizers_pos[i], minimizers_pos[i+k-1] + params.l - 1);
                let read_offsets = (read_obj.minimizers_pos[i] as usize, (read_obj.minimizers_pos[i+k-1] as usize + l - 1), (read_obj.minimizers_pos[i+k-1] + l - read_obj.minimizers_pos[i]));
                //println!("Shift: {:?}", shift);
                let mut span_seq = &read_obj.seq[shift.0..shift.1 + 1];
                //println!("Span_seq: {:?}, len: {:?}", span_seq, span_seq.len());
                //start of kminmer : pos[i]
                //end of kminmer : pos [i+k-1]+l
                let lookup = query_kminmer_lookup(kminmer, Some(span_seq), &seq_reversed, &seq_id, &shift, &mut sequences_file, thread_id, Some(&read_obj.seq), None);
                match lookup {
                    Some(x) => lookups.push(x),
                    None => continue
                }
            }
        }
        if lookups.len() > 0 {
            return Some((lookups, Some(read_obj), None))
        }
        else {return None}
    };
    let query_process_read_aux_sync = |seq_str: &[u8], seq_id: &str| -> Option<(Vec<QueryMatch>, Option<Read>, Option<ReadSync>)> {
        let thread_id :usize =  thread_id::get();
        let mut output : Vec<(Kmer, String, bool, String, (usize, usize))> = Vec::new();
        let seq = std::str::from_utf8(seq_str).unwrap(); // see https://docs.rs/seq_io/0.4.0-alpha.0/seq_io/fasta/index.html
        // those two next lines do a string copy of the read sequence, because in the case of a
        // reference, we'll need to remove newlines. also, that sequence will be moved to the
        // Read object. could be optimized later
        let mut read_obj = ReadSync::extract(&seq_id, seq.to_string(), &params, &minimizer_to_int, &int_to_minimizer, &uhs_bloom, &lcp_bloom);
        let mut lookups = Vec::<QueryMatch>::new();
        if sequences_files.get(&thread_id).is_none() {
            sequences_files.insert(thread_id, create_sequences_file(thread_id));
        }
        let mut sequences_file = sequences_files.get_mut(&thread_id).unwrap();
        if read_obj.transformed.len() > k {
            for i in 0..(read_obj.transformed.len() - k + 1) {
                let mut kminmer : Kmer = Kmer::make_from(&read_obj.transformed[i..i+k]);
                let mut seq_reversed = false;
                if revcomp_aware { 
                    let (kminmer_norm, reversed) = kminmer.normalize(); 
                    kminmer = kminmer_norm;
                    seq_reversed = reversed;
                } 
                let minimizers_pos = &read_obj.minimizers_pos;
                let position_of_first_minimizer = match seq_reversed {
                    true => minimizers_pos[i+k-1].1 + params.l - 1,
                    false => minimizers_pos[i].0
                };
                let position_of_last_minimizer = match seq_reversed {
                    true => minimizers_pos[i].1 + params.l - 1,
                    false => minimizers_pos[i+k-1].1
                };
                let mut shift = (minimizers_pos[i].0, minimizers_pos[i+k-1].1 + params.l - 1);
                let read_offsets = (read_obj.minimizers_pos[i].0 as usize, (read_obj.minimizers_pos[i+k-1].1 as usize + l - 1), (read_obj.minimizers_pos[i+k-1].1 + l - read_obj.minimizers_pos[i].0));
                //println!("Shift: {:?}", shift);
                let mut span_seq = &read_obj.seq[shift.0..shift.1 + 1];
                //println!("Span_seq: {:?}, len: {:?}", span_seq, span_seq.len());
                //start of kminmer : pos[i]
                //end of kminmer : pos [i+k-1]+l
                let lookup = query_kminmer_lookup(kminmer, Some(span_seq), &seq_reversed, &seq_id, &shift, &mut sequences_file, thread_id, Some(&read_obj.seq), None);
                match lookup {
                    Some(x) => lookups.push(x),
                    None => continue
                };
            }
        }
        if lookups.len() > 0 {
            return Some((lookups, None, Some(read_obj)))
        }
        else {return None}
    };
    let query_process_read_fasta = |record: seq_io::fasta::RefRecord, found: &mut Option<(Vec<QueryMatch>, Option<Read>, Option<ReadSync>)>| {
        let seq_str = record.seq(); 
        let seq_id = record.id().unwrap().to_string();
        if params.syncmer != 0 {*found = query_process_read_aux_sync(&seq_str, &seq_id);}
        else {*found = query_process_read_aux(&seq_str, &seq_id);}
    };
    let query_process_read_fastq = |record: seq_io::fastq::RefRecord, found: &mut Option<(Vec<QueryMatch>, Option<Read>, Option<ReadSync>)>| {
        let seq_str = record.seq(); 
        let seq_id = record.id().unwrap().to_string();
        if params.syncmer != 0 {*found = query_process_read_aux_sync(&seq_str, &seq_id);}
        else {*found = query_process_read_aux(&seq_str, &seq_id);}
    };


    // parallel fasta parsing, with a main thread that writes to disk and populates hash tables
    let mut main_thread = |found: &Option<(Vec<QueryMatch>, Option<Read>, Option<ReadSync>)>| { // runs in main thread
        nb_reads += 1;
        if found.is_some() {
            //println!("Received read in main thread, nb kmers: {}", vec.len());
            let debug_only_display_matches = false;
            let (query_lookups, read_obj, read_sync_obj) = found.as_ref().unwrap();
            if debug_only_display_matches {
                for qmatch in query_lookups.iter() {
                    println!("K-min-mer {:?}, target(s) {:?}, query {:?}", qmatch.kminmer.print_as_string(), qmatch.target, qmatch.query);
                } // kminmer: Kmer, target: Vec<HashEntry>, query: HashEntry, unique: bool};
            }
            if read_obj.is_some() {query_matches.insert(read_obj.as_ref().unwrap().id.to_string(), query_lookups.to_vec());}
            else if read_sync_obj.is_some(){query_matches.insert(read_sync_obj.as_ref().unwrap().id.to_string(), query_lookups.to_vec());}
        }
        // Some(value) will stop the reader, and the value will be returned.
        // In the case of never stopping, we need to give the compiler a hint about the
        // type parameter, thus the special 'turbofish' notation is needed.
        None::<()>
    };
    
    let buf = get_reader(&filename);
    println!("Parsing input sequences...");
    if fasta_reads {
        let reader = seq_io::fasta::Reader::new(buf);
        read_process_fasta_records(reader, threads as u32, queue_len, query_process_read_fasta, |record, found| {main_thread(found)});
    }
    else {
        let reader = seq_io::fastq::Reader::new(buf);
        read_process_fastq_records(reader, threads as u32, queue_len, query_process_read_fastq, |record, found| {main_thread(found)});
    }

    pb.finish_print("Finished lookup.");

    for mut sequences_file in sequences_files.iter_mut() {sequences_file.flush().unwrap();}
    let path = format!("{}{}", output_prefix.to_str().unwrap(), ".paf");
    let mut file = match File::create(&path) {
        Err(why) => panic!("Couldn't create {}: {}", path, why.description()),
        Ok(file) => file,
    };
    let query_matches_view = Arc::try_unwrap(query_matches).unwrap().into_read_only();
    let mut match_count = 0;
    for (id, entry) in query_matches_view.iter() {
        for qmatch in entry.iter() {
            let mut query_abund = *query_abundance_table.get(&qmatch.kminmer).unwrap();
            let mut ref_abund = *ref_abundance_table.get(&qmatch.kminmer).unwrap();
            //if query_abund == 1 {continue;}
            if query_abund == 1 || ref_abund > 1 {continue;}
            let query = &qmatch.query;
            let mut target = &qmatch.target[0];
            let mut query_span = query.shift.1 - query.shift.0;
            let mut target_span = target.shift.1 - target.shift.0;
            let ori = paf_output::determine_orientation(query, target);
            let mut block_len = target_span;
            println!("Query origin: {}\tTarget origin: {}\tQuery abundance: {:?}\tTarget abundance: {:?}\tQuery span: {}\tTarget span: {}", query.origin, target.origin, query_abund, ref_abund, query_span, target_span);
            let paf_line = format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n", query.origin, query.seqlen, query.shift.0, query.shift.1, ori, target.origin, target.seqlen, target.shift.0, target.shift.1, block_len, block_len, "255");
            write!(file, "{}", paf_line).expect("Error writing line.");
            match_count += 1;
            //target = &qmatch.target[1];}
            //let ori = paf_output::determine_orientation(query, target);
            //let paf_line = format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n", query.origin, query.seqlen, query.shift.0, query.shift.1, ori, target.origin, target.seqlen, target.shift.0, target.shift.1, residue_match, block_len, "255");
            //let (cigar, block_len, residue_match) = paf_output::align(query, target);
            //let paf_line = format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tcg:Z:{}\n", query.origin, query.seqlen, query.shift.0, query.shift.1, ori, target.origin, target.seqlen, target.shift.0, target.shift.1, residue_match, block_len, "255", cigar);
        }
    }
    println!("Number of reads: {}", nb_reads);
    println!("Found {} matches.", match_count);
    let duration = start.elapsed();
    println!("Total execution time: {:?}", duration);
    println!("Maximum RSS: {:?}GB", (get_memory_rusage() as f32) / 1024.0 / 1024.0 / 1024.0);
}

