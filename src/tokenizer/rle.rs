/*
Copyright (c) 2020 Pierre Marijon <pmarijon@mpi-inf.mpg.de>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
 */

/* local use */
use crate::{kmer, rle};

/// An iterator that takes a DNA sequence and produces kmers, in the forward orientation and 2bit form, homopolymer are compacted.
///
/// # Example
///
/// ```
/// use cocktail::tokenizer::TokenizerRLE;
///
/// let tokenizer = TokenizerRLE::new(b"GTACTGTGCCCGTGTTACTTAGTAAGCGTGAAAGGTGCGTGTTTCCGAGA", 5);
///
/// for rle_kmer in tokenizer {
///     // ... do what you want ...
/// }
pub struct TokenizerRLE {
    kmer_mask: u64,
    seq: Box<[u8]>,
    pos: usize,
    kmer: u64,
}

impl TokenizerRLE {
    /// Create a new TokenizerRLE on seq DNA kmer size is equal to k
    pub fn new(seq: &[u8], k: u8) -> Self {
        TokenizerRLE {
            kmer_mask: (1 << (k * 2)) - 1,
            seq: rle::seq2rle(seq),
            pos: (k - 1) as usize,
            kmer: kmer::seq2bit(&seq[0..((k - 1) as usize)]),
        }
    }
}

impl Iterator for TokenizerRLE {
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
mod test {
    use super::*;

    #[test]
    fn basic() {
        assert_eq!(
            vec![108, 433, 710, 795],
            TokenizerRLE::new(b"ACTGACTG", 5).collect::<Vec<u64>>()
        );
    }

    #[test]
    fn hash() {
        assert_eq!(
            vec![54, 457, 114, 397],
            TokenizerRLE::new(b"ACTGACTG", 5)
                .map(|x| crate::kmer::remove_first_bit(crate::kmer::canonical(x, 5)))
                .collect::<Vec<u64>>()
        );
    }

    #[test]
    fn canonical() {
        assert_eq!(
            vec![108, 915, 228, 795],
            TokenizerRLE::new(b"ACTGACTG", 5)
                .map(|x| crate::kmer::canonical(x, 5))
                .collect::<Vec<u64>>()
        );
    }
}
