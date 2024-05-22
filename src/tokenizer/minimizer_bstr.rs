//! Bstring minimizer tokenizer

/* std use */
use std::ops::BitXor;

/* crate use */

/* local use */
use crate::kmer;

/// An iterator that takes a DNA sequence and produces kmers bytes string (in the forward direction) and the associated minimizer in 2bit form.
///
/// # Example
///
/// ```
/// use cocktail::tokenizer::MiniBstr;
///
/// let tokenizer = MiniBstr::new(b"GTACTGTGCCCGTGTTACTTAGTAAGCGTGAAAGGTGCGTGTTTCCGAGA", 5, 3);
///
/// for (kmer, minimizer) in tokenizer {
///     // ... do what you want ...
/// }
pub struct MiniBstr<'a> {
    minimizers_mask: u64,
    seq: &'a [u8],
    pos: usize,
    minimizers: [u64; 2],
    mini_pos: usize,
    k: u64,
    m: u8,
}

impl<'a> MiniBstr<'a> {
    /// Create a new MiniBstr on seq DNA kmer size is equal to k, minimizer size is to m
    pub fn new(seq: &'a [u8], k: u64, m: u8) -> Self {
        let mut obj = Self {
            minimizers_mask: (1 << (m * 2)) - 1,
            seq,
            pos: 0,
            minimizers: [0, 0],
            mini_pos: 0,
            k,
            m,
        };

        obj.found_minimizer();

        obj
    }

    fn found_minimizer(&mut self) {
        let kmer = &self.seq[self.pos..(self.pos + self.k as usize)];

        let mut forward = unsafe { kmer::seq2bit(kmer.get_unchecked(0..self.m as usize)) };
        self.minimizers[0] = forward;
        let mut score = Self::get_score(kmer::canonical(forward, self.m));
        self.mini_pos = 0;

        for (i, nuc) in kmer[self.m as usize..].iter().enumerate() {
            let bits = kmer::nuc2bit(*nuc);
            forward = ((forward << 2) & self.minimizers_mask) | bits;

            let actual_score = Self::get_score(kmer::canonical(forward, self.m));

            if actual_score < score {
                score = actual_score;
                self.minimizers[0] = forward;
                self.mini_pos = i;
            }
        }
        self.minimizers[1] = kmer::revcomp(self.minimizers[0], self.m);
    }

    fn get_score(x: u64) -> u64 {
        x.rotate_left(5)
            .bitxor(x)
            .wrapping_mul(0x517c_c1b7_2722_0a95)
    }
}

impl<'a> Iterator for MiniBstr<'a> {
    type Item = (&'a [u8], u64);

    fn next(&mut self) -> Option<Self::Item> {
        if self.pos + self.k as usize - 1 == self.seq.len() {
            None
        } else {
            let kmer = &self.seq[self.pos..self.pos + self.k as usize];
            let mini = self.minimizers;

            self.pos += 1;

            if self.mini_pos == 0 {
                self.found_minimizer()
            } else {
                self.mini_pos -= 1;
            }

            if kmer::parity_even(mini[0]) {
                Some((kmer, mini[0]))
            } else {
                Some((kmer, mini[1]))
            }
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    const SEQ: &[u8] = b"AACTAAACGCTGTACGTACGAGCGATT";

    #[test]
    fn basic() {
        let value = MiniBstr::new(SEQ, 10, 5).collect::<Vec<(&[u8], u64)>>();

        for (kmer, mini) in &value {
            println!("{}", String::from_utf8(kmer.to_vec()).unwrap());
            println!("{}", kmer::kmer2seq(*mini, 5))
        }

        let truth: Vec<(&[u8], u64)> = vec![
            (b"AACTAAACGC", kmer::seq2bit(b"AACTA")),
            (b"ACTAAACGCT", kmer::seq2bit(b"CGTTT")),
            (b"CTAAACGCTG", kmer::seq2bit(b"CGTTT")),
            (b"TAAACGCTGT", kmer::seq2bit(b"CGTTT")),
            (b"AAACGCTGTA", kmer::seq2bit(b"CGTTT")),
            (b"AACGCTGTAC", kmer::seq2bit(b"GTACA")),
            (b"ACGCTGTACG", kmer::seq2bit(b"GTACA")),
            (b"CGCTGTACGT", kmer::seq2bit(b"GTACA")),
            (b"GCTGTACGTA", kmer::seq2bit(b"GTACA")),
            (b"CTGTACGTAC", kmer::seq2bit(b"GTACA")),
            (b"TGTACGTACG", kmer::seq2bit(b"GTACA")),
            (b"GTACGTACGA", kmer::seq2bit(b"TACGA")),
            (b"TACGTACGAG", kmer::seq2bit(b"TACGA")),
            (b"ACGTACGAGC", kmer::seq2bit(b"TACGA")),
            (b"CGTACGAGCG", kmer::seq2bit(b"TACGA")),
            (b"GTACGAGCGA", kmer::seq2bit(b"TACGA")),
            (b"TACGAGCGAT", kmer::seq2bit(b"GCGAT")),
            (b"ACGAGCGATT", kmer::seq2bit(b"GCGAT")),
        ];

        assert_eq!(value, truth);
    }
}
