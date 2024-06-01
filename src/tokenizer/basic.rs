//! A forward tokenizer

/* standard use */

/* crates use */

/* local use */
use crate::kmer;

/// An iterator that takes a DNA sequence and produces kmers, in the forward orientation and 2bit form.
///
/// # Example
///
/// ```
/// use cocktail::tokenizer::Tokenizer;
///
/// let tokenizer = Tokenizer::new(b"GTACTGTGCCCGTGTTACTTAGTAAGCGTGAAAGGTGCGTGTTTCCGAGA", 5);
///
/// for kmer in tokenizer {
///     // ... do what you want ...
/// }
pub struct Tokenizer<'a> {
    kmer_mask: u64,
    seq: &'a [u8],
    pos: usize,
    kmer: u64,
}

impl<'a> Tokenizer<'a> {
    /// Create a new Tokenizer on seq DNA kmer size is equal to k
    pub fn new(seq: &'a [u8], k: u8) -> Self {
        Tokenizer {
            kmer_mask: (1 << (k * 2)) - 1,
            seq,
            pos: (k - 1) as usize,
            kmer: kmer::seq2bit(unsafe { seq.get_unchecked(0..((k - 1) as usize)) }),
        }
    }
}

impl<'a> Iterator for Tokenizer<'a> {
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

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn basic() {
        assert_eq!(
            vec![108, 433, 710, 795],
            Tokenizer::new(b"ACTGACTG", 5).collect::<Vec<u64>>()
        );
    }

    #[test]
    fn seqlen_equal_k() {
        assert_eq!(vec![108], Tokenizer::new(b"ACTGA", 5).collect::<Vec<u64>>());
    }

    #[test]
    fn hash() {
        assert_eq!(
            vec![54, 457, 114, 397],
            Tokenizer::new(b"ACTGACTG", 5)
                .map(|x| crate::kmer::remove_first_bit(crate::kmer::canonical(x, 5)))
                .collect::<Vec<u64>>()
        );
    }

    #[test]
    fn canonical() {
        assert_eq!(
            vec![108, 915, 228, 795],
            Tokenizer::new(b"ACTGACTG", 5)
                .map(|x| crate::kmer::canonical(x, 5))
                .collect::<Vec<u64>>()
        );
    }
}
