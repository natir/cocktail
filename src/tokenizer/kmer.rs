//! This module provides iterator to produce kmer from DNA sequence

/* standard use */

/* crates use */

/* project use */
use crate::kmer;

/// An iterator that takes a DNA sequence and produces kmers, in the forward orientation and 2bit form.
///
/// # Example
///
/// ```
/// use cocktail::tokenizer::kmer::Forward;
///
/// let tokenizer = Forward::new(b"GTACTGTGCCCGTGTTACTTAGTAAGCGTGAAAGGTGCGTGTTTCCGAGA", 5);
///
/// for kmer in tokenizer {
///     // ... do what you want ...
/// }
pub struct Forward<'a> {
    kmer_mask: u64,
    seq: &'a [u8],
    pos: usize,
    kmer: u64,
}

impl<'a> Forward<'a> {
    /// Create a new Forward on seq DNA kmer size is equal to k
    pub fn new(seq: &'a [u8], k: u8) -> Self {
        Forward {
            kmer_mask: (1 << (k * 2)) - 1,
            seq,
            pos: (k - 1) as usize,
            kmer: kmer::seq2bit(unsafe { seq.get_unchecked(0..((k - 1) as usize)) }),
        }
    }
}

impl<'a> Iterator for Forward<'a> {
    type Item = u64;

    fn next(&mut self) -> Option<Self::Item> {
        if self.pos == self.seq.len() {
            None
        } else {
            self.kmer = unsafe {
                ((self.kmer << 2) & self.kmer_mask)
                    | kmer::nuc2bit(*self.seq.get_unchecked(self.pos))
            };

            self.pos += 1;

            Some(self.kmer)
        }
    }
}

/// An iterator that takes a DNA sequence and produces kmers, in the canonical orientation and 2bit form.
///
/// # Example
///
/// ```
/// use cocktail::tokenizer::kmer::Canonical;
///
/// let tokenizer = Canonical::new(b"GTACTGTGCCCGTGTTACTTAGTAAGCGTGAAAGGTGCGTGTTTCCGAGA", 5);
///
/// for kmer in tokenizer {
///     // ... do what you want ...
/// }
pub struct Canonical<'a> {
    move_bit: u8,
    kmer_mask: u64,
    seq: &'a [u8],
    pos: usize,
    kmers: [u64; 2],
}

impl<'a> Canonical<'a> {
    /// Create a new Canonical tokenizer on seq DNA, kmer size is equal to k
    pub fn new(seq: &'a [u8], k: u8) -> Self {
        let forward = unsafe { kmer::seq2bit(seq.get_unchecked(0..((k - 1) as usize))) };

        Canonical {
            move_bit: (k - 1) * 2,
            kmer_mask: (1 << (k * 2)) - 1,
            seq,
            pos: (k - 1) as usize,
            kmers: [forward, kmer::revcomp(forward, k)],
        }
    }
}

impl<'a> Iterator for Canonical<'a> {
    type Item = u64;

    fn next(&mut self) -> Option<Self::Item> {
        if self.pos == self.seq.len() {
            None
        } else {
            unsafe {
                let nuc = kmer::nuc2bit(*self.seq.get_unchecked(self.pos));
                self.pos += 1;

                *self.kmers.get_unchecked_mut(0) =
                    ((self.kmers.get_unchecked(0) << 2) & self.kmer_mask) | nuc;
                *self.kmers.get_unchecked_mut(1) =
                    (*self.kmers.get_unchecked(1) >> 2) ^ ((nuc ^ 0b10) << self.move_bit);

                if kmer::parity_even(*self.kmers.get_unchecked(0)) {
                    Some(*self.kmers.get_unchecked(0))
                } else {
                    Some(*self.kmers.get_unchecked(1))
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn forward() {
        assert_eq!(
            vec![108, 433, 710, 795],
            Forward::new(b"ACTGACTG", 5).collect::<Vec<u64>>()
        );
    }

    #[test]
    fn forward_equal_k() {
        assert_eq!(vec![108], Forward::new(b"ACTGA", 5).collect::<Vec<u64>>());
    }

    #[test]
    fn forward_hash() {
        assert_eq!(
            vec![54, 457, 114, 397],
            Forward::new(b"ACTGACTG", 5)
                .map(|x| crate::kmer::remove_first_bit(crate::kmer::canonical(x, 5)))
                .collect::<Vec<u64>>()
        );
    }

    #[test]
    fn canonical() {
        assert_eq!(
            vec![108, 915, 228, 795],
            Canonical::new(b"ACTGACTG", 5).collect::<Vec<u64>>()
        );
    }

    #[test]
    fn canonical_only_check() {
        let token = Canonical::new(b"ACTGACTG", 5);

        for cano in token {
            assert!(kmer::parity_even(cano));
        }
    }
}
