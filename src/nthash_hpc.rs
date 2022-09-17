// homopolymer compression version of ntHash1,
// which simulatenously HPC and record original offsets of sequences
//
//
// single-file version
// adapted from ntHash1 and Luiz Irber's Rust crate port
// kept only the canonical version (no 'forward')


const MAXIMUM_K_SIZE: usize = u32::max_value() as usize;

#[derive(thiserror::Error, Debug)]
pub enum Error {
    #[error("K size {ksize} is out of range for the given sequence size {seq_size}")]
    KSizeOutOfRange { ksize: usize, seq_size: usize },
    #[error("K size {0} cannot exceed the size of a u32 {MAXIMUM_K_SIZE}")]
    KSizeTooBig(usize),
}

pub type Result<T> = std::result::Result<T, Error>;

const H_LOOKUP: [u64; 256] = {
    let mut lookup = [1; 256];
    lookup[b'A' as usize] = 0x3c8b_fbb3_95c6_0474;
    lookup[b'C' as usize] = 0x3193_c185_62a0_2b4c;
    lookup[b'G' as usize] = 0x2032_3ed0_8257_2324;
    lookup[b'T' as usize] = 0x2955_49f5_4be2_4456;
    lookup[b'N' as usize] = 0;
    lookup
};

const RC_LOOKUP: [u64; 256] = {
    let mut lookup = [1; 256];
    lookup[b'A' as usize] = 0x2955_49f5_4be2_4456;
    lookup[b'C' as usize] = 0x2032_3ed0_8257_2324;
    lookup[b'G' as usize] = 0x3193_c185_62a0_2b4c;
    lookup[b'T' as usize] = 0x3c8b_fbb3_95c6_0474;
    lookup[b'N' as usize] = 0;
    lookup
};

#[inline(always)]
fn h(c: u8) -> u64 {
    let val = H_LOOKUP[c as usize];
    if val == 1 {
        panic!("Non-ACGTN nucleotide encountered! {}", c as char)
    }
    val
}

#[inline(always)]
fn rc(nt: u8) -> u64 {
    let val = RC_LOOKUP[nt as usize];
    if val == 1 {
        panic!("Non-ACGTN nucleotide encountered! {}", nt as char)
    }
    val
}

/// An efficient iterator for returning rev-comp aware hashes in HPC space, keeping only those
/// below a hash bound. That's a very specialized application, useful for rust-mdbg.
///
/// Since it implements the `Iterator` trait it also
/// exposes many other useful methods. In this example we use `collect` to
/// generate all hashes and put them in a `Vec<u64>`.
/// ```
///     # use nthash::Result;
///     use nthash::NtHashHPCIterator;
///
///     # fn main() -> Result<()> {
///     let seq = b"ACTGC";
///     let hash_bound = 0.1;
///     let iter = NtHashHPCIterator::new(seq, 3, hash_bound)?;
///     let hashes: Vec<u64> = iter.collect();
///     assert_eq!(hashes,
///                vec![0x9b1eda9a185413ce, 0x9f6acfa2235b86fc, 0xd4a29bf149877c5c]);
///     # Ok(())
///     # }
/// ```

#[derive(Debug)]
pub struct NtHashHPCIterator<'a> {
    seq: &'a [u8],
    k: usize,
    fh: u64,
    rh: u64,
    current_idx: usize,
    current_idx_plus_k: usize,
    seq_len: usize,
    hash_bound: u64,
}

impl<'a> NtHashHPCIterator<'a> {
    /// Creates a new NtHashHPCIterator with internal state properly initialized.
    pub fn new(seq: &'a [u8], k: usize, hash_bound: u64) -> Result<NtHashHPCIterator<'a>> {
        let seq_len = seq.len();
        if k > seq_len {
            return Err(Error::KSizeOutOfRange {
                ksize: k,
                seq_size: seq.len(),
            });
        }
        if k > MAXIMUM_K_SIZE {
            return Err(Error::KSizeTooBig(k));
        }
        let mut fh = 0;
        let mut j = 0; // current position in non-HPC space
        let mut i = 0; // current position in HPC space
        let mut prev;
        let mut prev_j = 0;
        let mut v;

        // pre-compute hash of first kmer in HPC space
        while i < k && j < seq_len
        {
            v = seq[j];
            fh ^= h(v).rotate_left((k - i - 1) as u32);
            i += 1;
            // HPC on the fly
            prev = v;
            prev_j = j;
            while j < seq_len &&  seq[j] == prev { j += 1};
        }
        // if sequence is shorter than k in HPC space: hash will be incorrect and i<k , but that's ok, 
        // hash is only returned later in the iterator where proper length checks are performed
        assert!( (j >= seq_len && i < k) || (j < seq_len && i==k) );

        // at this point just assume we're at position i==k in HPC space, j in non-HPC space,
        // read the sequence backwards to get hash of reverse complement
        i -= 1; 
        j = prev_j;
        let k_in_non_hpc_space = j;
        let mut rh = 0;
        while i >= 0 && j >= 0
        {
            v = seq[j];
            rh ^= rc(v).rotate_left(i as u32);
            if i == 0 { break;}
            i -= 1;
            // HPC on the fly
            prev = v;
            while j > 0 && seq[j] == prev { j -= 1};
        }

        Ok(NtHashHPCIterator {
            seq,
            k,
            fh,
            rh,
            current_idx: 0,
            current_idx_plus_k: k_in_non_hpc_space,
            seq_len,
            hash_bound
        })
    }
}

impl<'a> Iterator for NtHashHPCIterator<'a> {
    type Item = (usize,u64);


    fn next(&mut self) -> Option<(usize,u64)> {
        let mut hash = 0;
        let mut prev_current_idx = 0;
        loop
        {
            if self.current_idx_plus_k >= self.seq_len - 1 {
                return None;
            };

            if self.current_idx != 0 {
                let i = self.current_idx - 1;
                let seqi = self.seq[i];
                let seqk = self.seq[self.current_idx_plus_k];

                self.fh = self.fh.rotate_left(1) ^ h(seqi).rotate_left(self.k as u32) ^ h(seqk);

                self.rh = self.rh.rotate_right(1)
                    ^ rc(seqi).rotate_right(1)
                    ^ rc(seqk).rotate_left(self.k as u32 - 1);
            }
            
            hash = u64::min(self.rh, self.fh);

            // update (current,current+k) pointers to next positions in HPC space
            let mut prev = self.seq[self.current_idx_plus_k];
            while self.current_idx_plus_k < self.seq_len  && self.seq[self.current_idx_plus_k] == prev 
            {
                self.current_idx_plus_k += 1;
            }
            prev = self.seq[self.current_idx];
            prev_current_idx = self.current_idx;
            while  self.current_idx < self.seq_len && self.seq[self.current_idx] == prev 
            {
                self.current_idx += 1;
            }

            if hash <= self.hash_bound { break; }
        }
        Some((prev_current_idx, hash))
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let size_hint = (self.seq_len - self.k + 1) * (self.hash_bound as u64) as usize; // rough estimation
        (size_hint, Some(size_hint)) 
    }
}
