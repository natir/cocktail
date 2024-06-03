//! An iterator than produce RLE kmer

/* standard use */

/* crates use */

/* project use */

/* module declaration */
use crate::kmer;
use crate::rle;

/// An iterator that takes a DNA sequence and produces kmers, in the forward orientation and 2bit form, homopolymer are compacted.
///
/// # Example
///
/// ```
/// use cocktail::tokenizer::rle::Forward;
///
/// let tokenizer = Forward::new(b"GTACTGTGCCCGTGTTACTTAGTAAGCGTGAAAGGTGCGTGTTTCCGAGA", 5);
///
/// for rle_kmer in tokenizer {
///     // ... do what you want ...
/// }
pub struct Forward {
    kmer_mask: u64,
    seq: Box<[u8]>,
    pos: usize,
    kmer: u64,
}

impl Forward {
    /// Create a new TokenizerRLE on seq DNA kmer size is equal to k
    pub fn new(seq: &[u8], k: u8) -> Self {
        Forward {
            kmer_mask: (1 << (k * 2)) - 1,
            seq: rle::seq2rle(seq),
            pos: (k - 1) as usize,
            kmer: kmer::seq2bit(unsafe { seq.get_unchecked(0..((k - 1) as usize)) }),
        }
    }
}

impl Iterator for Forward {
    type Item = u64;

    fn next(&mut self) -> Option<Self::Item> {
        if self.pos == self.seq.len() {
            None
        } else {
            self.kmer = ((self.kmer << 2) & self.kmer_mask) | rle::rle2bit(self.seq[self.pos]);

            self.pos += 1;

            Some(self.kmer)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn basic() {
        assert_eq!(
            vec![108, 433, 710, 795],
            Forward::new(b"ACTGACTG", 5).collect::<Vec<u64>>()
        );
    }

    #[test]
    fn hash() {
        assert_eq!(
            vec![54, 457, 114, 397],
            Forward::new(b"ACTGACTG", 5)
                .map(|x| kmer::remove_first_bit(kmer::canonical(x, 5)))
                .collect::<Vec<u64>>()
        );
    }

    #[test]
    fn canonical() {
        assert_eq!(
            vec![108, 915, 228, 795],
            Forward::new(b"ACTGACTG", 5)
                .map(|x| kmer::canonical(x, 5))
                .collect::<Vec<u64>>()
        );
    }
}
