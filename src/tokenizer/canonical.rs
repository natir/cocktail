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
use crate::kmer;

/// An iterator that takes a DNA sequence and produces kmers, in the forward orientation and 2bit form.
///
/// # Example
///
/// ```
/// use cocktail::tokenizer::Canonical;
///
/// let tokenizer = Canonical::new(b"GTACTGTGCCCGTGTTACTTAGTAAGCGTGAAAGGTGCGTGTTTCCGAGA", 5);
///
/// for kmer in tokenizer {
///     // ... do what you want ...
/// }
pub struct Canonical<'a> {
    k: u8,
    kmer_mask: u64,
    seq: &'a [u8],
    pos: usize,
    forward: u64,
    reverse: u64,
    forward_cano: bool,
}

impl<'a> Canonical<'a> {
    pub fn new(seq: &'a [u8], k: u8) -> Self {
        let forward = kmer::seq2bit(&seq[0..((k - 1) as usize)]);

        Canonical {
            k,
            kmer_mask: (1 << (k * 2)) - 1,
            seq,
            pos: (k - 1) as usize,
            forward,
            reverse: kmer::revcomp(forward, k),
            forward_cano: kmer::parity_even(forward),
        }
    }
}

impl<'a> Iterator for Canonical<'a> {
    type Item = u64;

    fn next(&mut self) -> Option<Self::Item> {
        if self.pos == self.seq.len() {
            None
        } else {
            let nuc = kmer::nuc2bit(self.seq[self.pos]);
            self.pos += 1;

            let chg_cano = ((nuc ^ self.forward >> ((self.k - 1) * 2)).count_ones() % 2) != 0;

            self.forward = ((self.forward << 2) & self.kmer_mask) | nuc;
            self.reverse = self.reverse >> 2 ^ ((nuc ^ 0b10) << ((self.k - 1) * 2));

            self.forward_cano = if chg_cano {
                !self.forward_cano
            } else {
                self.forward_cano
            };

            if self.forward_cano {
                Some(self.forward)
            } else {
                Some(self.reverse)
            }
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn basic() {
        assert_eq!(
            vec![108, 915, 228, 795],
            Canonical::new(b"ACTGACTG", 5).collect::<Vec<u64>>()
        );
    }

    #[test]
    fn only_canonical() {
        let token = Canonical::new(b"ACTGACTG", 5);

        for cano in token {
            assert!(kmer::parity_even(cano));
        }
    }
}
