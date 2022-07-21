use std::borrow::Borrow;

lazy_static! {
    static ref COMPLEMENT: [u8; 256] = {
        let mut comp = [0; 256];
        for (v, a) in comp.iter_mut().enumerate() {
            *a = v as u8;
        }
        for (&a, &b) in b"ACGT".iter().zip(b"TGCA".iter()) {
            comp[a as usize] = b;
            comp[a as usize + 32] = b + 32;  // lowercase variants
        }
        comp
    };
}

pub fn complement(a: u8) -> u8 {
    COMPLEMENT[a as usize]
}

pub fn revcomp<C, T>(text: T) -> Vec<u8>
where
    C: Borrow<u8>,
    T: IntoIterator<Item = C>,
    T::IntoIter: DoubleEndedIterator,
{
    text.into_iter()
        .rev()
        .map(|a| complement(*a.borrow()))
        .collect()
}

pub fn canonical(kmer: &[u8]) -> Vec<u8> {
    let rev = revcomp(kmer);
    if rev < kmer.to_vec() {return rev;}
    return kmer.to_vec();
}
