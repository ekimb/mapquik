// stats.rs
// Contains the "Stats" struct, to collect statistics on the potential reference locations of reads
//
// heuristic: collect all potential reference locations given by the kminmers of the read, 
// sort those locations, find the number of 'jumps' as defined by two consecutive locations whose 
// distance is > 2x read length

use std::{fs::File, io::Write, mem::MaybeUninit, sync::Mutex};
use crate::{Entry};
use fxhash::{hash32};

static ENABLED: bool = false;

// some rust black magic to have multithreaded static file handle
// from https://stackoverflow.com/a/72187853
static mut STATS_FILE: MaybeUninit<Mutex<File>> = MaybeUninit::uninit();

#[derive(Clone, Debug, PartialEq, Default)]
pub struct Stats {
    pub q_id: String,
    pub ref_loci: Vec<(u32,usize)>, // reference locis for analyzed read
}

impl Stats {

    pub fn init(nb_threads: usize, output_prefix: &str)
    {
        if ENABLED
        {
            if nb_threads != 1
            {
                //panic!("If stats are enabled, mapquik can only be run in 1 thread");
                //seems code works with multithreads..
            }

            // Stats file generation
            let stats_path = format!("{}{}", output_prefix, ".read_stats");
            unsafe
            {
            STATS_FILE = match File::create(&stats_path) {
                Err(why) => panic!("Couldn't create {}: {}", stats_path, why.to_string()),
                Ok(stats_file) =>  MaybeUninit::new(Mutex::new(stats_file)),
            };
            }
            println!("Stats module initialized.");
        }
    }

    // A new Stats object for a read
    pub fn new(q_id: &str) -> Self {
        if ENABLED
        {
            let s = Stats {
                q_id: q_id.to_string(),
                ref_loci: vec![]
            };
            s
        }
        else { Stats::default() }
    }
    
    // Add a new potential reference location
    pub fn add(&mut self, r: &Entry)
    {
        if ENABLED
        {
            self.ref_loci.push((hash32(&r.id),r.start));
        }
    }

    // Compute number of jumps and write to stats file
    pub fn finalize(&mut self)
    {
        if ENABLED
        {
            self.ref_loci.sort();
            let mut prev: (u32, usize) = (0,0);
            let dist = 48000; // expected minimal distance between unrelated regions (2x average read length)
            let mut nb_loci = 0;
            for (a,b) in &self.ref_loci
            {
                if *a != prev.0  || *b-prev.1 > dist
                {
                    nb_loci += 1;
                }
                prev = (*a,*b);
            }
            let stats_line = format!("{}: {}",self.q_id,nb_loci);
            unsafe
            {
            write!(STATS_FILE.assume_init_mut().lock().unwrap(), "{}\n", stats_line).expect("Error writing stats line.");
            }
        }
    }

}
