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

/* standard use */
use std::ops::BitXor;

/* local use */
use crate::kmer;

/// An iterator that takes a DNA sequence and produces kmers (in the forward direction and 2bit form) and the associated minimizer.
///
/// # Example
///
/// ```
/// use cocktail::tokenizer::TokenizerMini;
///
/// let tokenizer = TokenizerMini::new(b"GTACTGTGCCCGTGTTACTTAGTAAGCGTGAAAGGTGCGTGTTTCCGAGA", 5, 3);
///
/// for (kmer, minimizer) in tokenizer {
///     // ... do what you want ...
/// }
pub struct TokenizerMini<'a> {
    kmer_mask: u64,
    seq: &'a [u8],
    pos: usize,
    kmer: u64,
    minimizer: MinimizerRing,
}

impl<'a> TokenizerMini<'a> {
    /// Create a new TokenizerMini on seq DNA kmer size is equal to k
    pub fn new(seq: &'a [u8], k: u8, m: u8) -> Self {
        let kmer = kmer::seq2bit(&seq[0..((k - 1) as usize)]);

        Self {
            kmer_mask: (1 << (k * 2)) - 1,
            seq,
            pos: (k - 1) as usize,
            kmer,
            minimizer: MinimizerRing::new(k, m, kmer),
        }
    }
}

impl<'a> Iterator for TokenizerMini<'a> {
    type Item = (u64, u64);

    fn next(&mut self) -> Option<Self::Item> {
        if self.pos == self.seq.len() {
            None
        } else {
            self.kmer = ((self.kmer << 2) & self.kmer_mask) | kmer::nuc2bit(self.seq[self.pos]);

            self.minimizer.add_kmer(self.kmer);

            self.pos += 1;

            Some((self.kmer, self.minimizer.get_mini().0))
        }
    }
}

/// A struct to get minimizer of sucessive kmer
///
/// At initialization all subkmer with weight is compute and store in a ring buffer.
/// When the next kmer is add only the new subkmer and is weight is compute.
/// If the new subkmer erase the previous minimizer but is score isn't lower than previous minimizer, the ring buffer is scanned completely to find the new minimizer.
pub struct MinimizerRing {
    ring_buffer: Box<[(u64, u64)]>,
    current: usize,
    minimizer: usize,
    mask: u64,
    k: u8,
    m: u8,
}

impl MinimizerRing {
    /// Create a MinimizerRing, with kmer size equale to `k`, subkmer size equale to `m` and init ring buffer with `kmer`
    pub fn new(k: u8, m: u8, kmer: u64) -> Self {
        let mut obj = Self {
            ring_buffer: vec![(0, 0); (k - m + 1) as usize].into_boxed_slice(),
            current: 0,
            minimizer: 0,
            mask: (1 << (m * 2)) - 1,
            k,
            m,
        };

        obj.populate_buffer(kmer);

        obj
    }

    /// Reset the ring buffer with a new kmer
    pub fn populate_buffer(&mut self, mut kmer: u64) {
        let mut score = u64::max_value();
        let max_len = (self.k - self.m + 1) as usize;

        for i in 0..max_len {
            let rb_index = (max_len - i - 1) as usize;

            let mini = kmer::cannonical(kmer & self.mask, self.m);

            let local_score = MinimizerRing::get_score(mini);
            self.ring_buffer[rb_index] = (mini, local_score);

            if local_score < score {
                score = local_score;
                self.minimizer = rb_index;
            }

            kmer >>= 2;
        }

        self.current = 0;
    }

    /// Add the next kmer
    pub fn add_kmer(&mut self, kmer: u64) {
        let minimizer = kmer::cannonical(kmer & self.mask, self.m);
        let score = MinimizerRing::get_score(minimizer);

        let previous_mini = self.get_mini();
        self.ring_buffer[self.current] = (minimizer, score);

        if score < previous_mini.1 {
            self.minimizer = self.current;
        } else if self.current == self.minimizer {
            self.update_minimizer();
        }

        self.current = (self.current + 1) % self.ring_buffer.len();
    }

    /// Get a pair of value first one is the minimizer second one is his score
    pub fn get_mini(&self) -> (u64, u64) {
        self.ring_buffer[self.minimizer]
    }

    fn update_minimizer(&mut self) {
        let mut index: usize = 0;
        let mut min = u64::max_value();

        for (i, elt) in self.ring_buffer.iter().enumerate() {
            if elt.1 < min {
                index = i;
                min = elt.1;
            }
        }

        self.minimizer = index;
    }

    fn get_score(x: u64) -> u64 {
        x.rotate_left(5)
            .bitxor(x)
            .wrapping_mul(0x517c_c1b7_2722_0a95)
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn basic() {
        let seq = b"acgtgcgtgagaattgttcggtggaacaggcg";

        let mut token = TokenizerMini::new(seq, 11, 7);

        let mut kmers = Vec::new();
        let mut canos = Vec::new();
        let mut minis = Vec::new();

        while let Some((kmer, mini)) = token.next() {
            let cano = crate::kmer::cannonical(kmer, 11);

            kmers.push(kmer);
            canos.push(cano);
            minis.push(mini);
        }

        assert_eq!(
            kmers,
            [
                505779, 2023116, 3898160, 3009730, 3650314, 2018347, 3879086, 2933434, 3345129,
                797607, 3190431, 178814, 715259, 2861039, 3055548, 3833584, 2751425, 2617092,
                2079763, 4124751, 3916093, 3081463
            ]
        );

        assert_eq!(
            canos,
            [
                505779, 2023116, 2724305, 3009730, 3650314, 2018347, 3879086, 68196, 3162777,
                1839270, 1508393, 178814, 715259, 2861039, 2430724, 3833584, 3821936, 2617092,
                1811735, 4124751, 3521105, 1928852
            ]
        );

        assert_eq!(
            minis,
            [
                4561, 4561, 4561, 10641, 13066, 13066, 13066, 698, 698, 698, 698, 698, 7184, 5212,
                5212, 5212, 5212, 5212, 10565, 10565, 5865, 1271
            ]
        );
    }

    #[test]
    fn same_in_each_strand() {
        let fwd = b"CACTCCTGTCACATCATAATCGTTTGCTATT";
        let rev = b"AATAGCAAACGATTATGATGTGACAGGAGTG";

        let mut fwd_token = TokenizerMini::new(fwd, 11, 7);
        let mut fwd_canos = Vec::new();
        let mut fwd_minis = Vec::new();

        let mut rev_token = TokenizerMini::new(rev, 11, 7);
        let mut rev_canos = Vec::new();
        let mut rev_minis = Vec::new();

        while let Some((kmer, mini)) = fwd_token.next() {
            let cano = crate::kmer::cannonical(kmer, 11);

            fwd_canos.push(cano);
            fwd_minis.push(mini);
        }

        while let Some((kmer, mini)) = rev_token.next() {
            let cano = crate::kmer::cannonical(kmer, 11);

            rev_canos.push(cano);
            rev_minis.push(mini);
        }

        if !fwd_canos.iter().eq(rev_canos.iter().rev()) {
            panic!("\nleft: {:?}\nright: {:?}", fwd_canos, rev_canos);
        }

        if !fwd_minis.iter().eq(rev_minis.iter().rev()) {
            panic!("\nleft: {:?}\nright: {:?}", fwd_minis, rev_minis);
        }
    }
}
