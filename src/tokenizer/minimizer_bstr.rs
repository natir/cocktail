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
        let kmer_range = self.pos..(self.pos + self.k as usize);

        if kmer_range.end > self.seq.len() {
            return;
        }

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
            println!(
                "(b\"{}\", kmer::seq2bit(b\"{}\"))",
                String::from_utf8(kmer.to_vec()).unwrap(),
                kmer::kmer2seq(*mini, 5)
            );
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

    const SEQ_LONG: &[u8] = b"AGGATAGAAGCTTAAGTACAAGATAATTCCCATAGAGGAAGGGTGGTATTACAGTGCCGCCTGTTGAAAGCCCCAATCCCGCTTCAATTGTTGAGCTCAG";

    #[test]
    fn long() {
        let value = MiniBstr::new(SEQ_LONG, 10, 5).collect::<Vec<(&[u8], u64)>>();

        for (kmer, mini) in &value {
            println!(
                "(b\"{}\", kmer::seq2bit(b\"{}\")),",
                String::from_utf8(kmer.to_vec()).unwrap(),
                kmer::kmer2seq(*mini, 5)
            );
        }

        let truth: Vec<(&[u8], u64)> = vec![
            (b"AGGATAGAAG", kmer::seq2bit(b"ATCCT")),
            (b"GGATAGAAGC", kmer::seq2bit(b"CTATC")),
            (b"GATAGAAGCT", kmer::seq2bit(b"CTATC")),
            (b"ATAGAAGCTT", kmer::seq2bit(b"AAGCT")),
            (b"TAGAAGCTTA", kmer::seq2bit(b"AAGCT")),
            (b"AGAAGCTTAA", kmer::seq2bit(b"AAGCT")),
            (b"GAAGCTTAAG", kmer::seq2bit(b"AAGCT")),
            (b"AAGCTTAAGT", kmer::seq2bit(b"TTAAG")),
            (b"AGCTTAAGTA", kmer::seq2bit(b"TTAAG")),
            (b"GCTTAAGTAC", kmer::seq2bit(b"TTAAG")),
            (b"CTTAAGTACA", kmer::seq2bit(b"AGTAC")),
            (b"TTAAGTACAA", kmer::seq2bit(b"AGTAC")),
            (b"TAAGTACAAG", kmer::seq2bit(b"AGTAC")),
            (b"AAGTACAAGA", kmer::seq2bit(b"AGTAC")),
            (b"AGTACAAGAT", kmer::seq2bit(b"TCTTG")),
            (b"GTACAAGATA", kmer::seq2bit(b"TCTTG")),
            (b"TACAAGATAA", kmer::seq2bit(b"TCTTG")),
            (b"ACAAGATAAT", kmer::seq2bit(b"TCTTG")),
            (b"CAAGATAATT", kmer::seq2bit(b"TCTTG")),
            (b"AAGATAATTC", kmer::seq2bit(b"AATTA")),
            (b"AGATAATTCC", kmer::seq2bit(b"AATTA")),
            (b"GATAATTCCC", kmer::seq2bit(b"AATTA")),
            (b"ATAATTCCCA", kmer::seq2bit(b"AATTA")),
            (b"TAATTCCCAT", kmer::seq2bit(b"AATTA")),
            (b"AATTCCCATA", kmer::seq2bit(b"GGGAA")),
            (b"ATTCCCATAG", kmer::seq2bit(b"GGGAA")),
            (b"TTCCCATAGA", kmer::seq2bit(b"GGGAA")),
            (b"TCCCATAGAG", kmer::seq2bit(b"CATAG")),
            (b"CCCATAGAGG", kmer::seq2bit(b"CATAG")),
            (b"CCATAGAGGA", kmer::seq2bit(b"CATAG")),
            (b"CATAGAGGAA", kmer::seq2bit(b"AGGAA")),
            (b"ATAGAGGAAG", kmer::seq2bit(b"AGGAA")),
            (b"TAGAGGAAGG", kmer::seq2bit(b"AGGAA")),
            (b"AGAGGAAGGG", kmer::seq2bit(b"AGGAA")),
            (b"GAGGAAGGGT", kmer::seq2bit(b"AGGAA")),
            (b"AGGAAGGGTG", kmer::seq2bit(b"AGGAA")),
            (b"GGAAGGGTGG", kmer::seq2bit(b"AAGGG")),
            (b"GAAGGGTGGT", kmer::seq2bit(b"AAGGG")),
            (b"AAGGGTGGTA", kmer::seq2bit(b"TGGTA")),
            (b"AGGGTGGTAT", kmer::seq2bit(b"TGGTA")),
            (b"GGGTGGTATT", kmer::seq2bit(b"TGGTA")),
            (b"GGTGGTATTA", kmer::seq2bit(b"TGGTA")),
            (b"GTGGTATTAC", kmer::seq2bit(b"TGGTA")),
            (b"TGGTATTACA", kmer::seq2bit(b"TGTAA")),
            (b"GGTATTACAG", kmer::seq2bit(b"TGTAA")),
            (b"GTATTACAGT", kmer::seq2bit(b"TGTAA")),
            (b"TATTACAGTG", kmer::seq2bit(b"TGTAA")),
            (b"ATTACAGTGC", kmer::seq2bit(b"TGTAA")),
            (b"TTACAGTGCC", kmer::seq2bit(b"TGTAA")),
            (b"TACAGTGCCG", kmer::seq2bit(b"CAGTG")),
            (b"ACAGTGCCGC", kmer::seq2bit(b"CAGTG")),
            (b"CAGTGCCGCC", kmer::seq2bit(b"CAGTG")),
            (b"AGTGCCGCCT", kmer::seq2bit(b"CGGCA")),
            (b"GTGCCGCCTG", kmer::seq2bit(b"CGGCA")),
            (b"TGCCGCCTGT", kmer::seq2bit(b"CAGGC")),
            (b"GCCGCCTGTT", kmer::seq2bit(b"CAGGC")),
            (b"CCGCCTGTTG", kmer::seq2bit(b"CAGGC")),
            (b"CGCCTGTTGA", kmer::seq2bit(b"CAGGC")),
            (b"GCCTGTTGAA", kmer::seq2bit(b"CAACA")),
            (b"CCTGTTGAAA", kmer::seq2bit(b"CAACA")),
            (b"CTGTTGAAAG", kmer::seq2bit(b"CAACA")),
            (b"TGTTGAAAGC", kmer::seq2bit(b"CAACA")),
            (b"GTTGAAAGCC", kmer::seq2bit(b"AAGCC")),
            (b"TTGAAAGCCC", kmer::seq2bit(b"AAGCC")),
            (b"TGAAAGCCCC", kmer::seq2bit(b"AAGCC")),
            (b"GAAAGCCCCA", kmer::seq2bit(b"AAGCC")),
            (b"AAAGCCCCAA", kmer::seq2bit(b"AAGCC")),
            (b"AAGCCCCAAT", kmer::seq2bit(b"GGGCT")),
            (b"AGCCCCAATC", kmer::seq2bit(b"GGGCT")),
            (b"GCCCCAATCC", kmer::seq2bit(b"GCCCC")),
            (b"CCCCAATCCC", kmer::seq2bit(b"ATTGG")),
            (b"CCCAATCCCG", kmer::seq2bit(b"ATTGG")),
            (b"CCAATCCCGC", kmer::seq2bit(b"CCCGC")),
            (b"CAATCCCGCT", kmer::seq2bit(b"CCCGC")),
            (b"AATCCCGCTT", kmer::seq2bit(b"CCCGC")),
            (b"ATCCCGCTTC", kmer::seq2bit(b"CCCGC")),
            (b"TCCCGCTTCA", kmer::seq2bit(b"CCCGC")),
            (b"CCCGCTTCAA", kmer::seq2bit(b"CCCGC")),
            (b"CCGCTTCAAT", kmer::seq2bit(b"ATTGA")),
            (b"CGCTTCAATT", kmer::seq2bit(b"ATTGA")),
            (b"GCTTCAATTG", kmer::seq2bit(b"ATTGA")),
            (b"CTTCAATTGT", kmer::seq2bit(b"ATTGA")),
            (b"TTCAATTGTT", kmer::seq2bit(b"ATTGA")),
            (b"TCAATTGTTG", kmer::seq2bit(b"TTGTT")),
            (b"CAATTGTTGA", kmer::seq2bit(b"TTGTT")),
            (b"AATTGTTGAG", kmer::seq2bit(b"TTGTT")),
            (b"ATTGTTGAGC", kmer::seq2bit(b"TTGTT")),
            (b"TTGTTGAGCT", kmer::seq2bit(b"TTGTT")),
            (b"TGTTGAGCTC", kmer::seq2bit(b"TGAGC")),
            (b"GTTGAGCTCA", kmer::seq2bit(b"TGAGC")),
            (b"TTGAGCTCAG", kmer::seq2bit(b"TGAGC")),
        ];

        assert_eq!(value, truth);
    }
}
