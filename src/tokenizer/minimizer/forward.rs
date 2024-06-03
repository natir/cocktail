//! 2bit minimizer forward tokenizer

/* standard use */

/* crates use */

/* local use */
use crate::kmer;
use crate::tokenizer::minimizer::method;

/// An iterator that takes a DNA sequence and produces kmers (in the forward direction and 2bit form) and the associated minimizer.
///
/// # Example
///
/// ```
/// use cocktail::tokenizer::minimizer::Forward;
/// use cocktail::tokenizer::minimizer::method;
///
/// let tokenizer = Forward::<method::Random>::new(b"GTACTGTGCCCGTGTTACTTAGTAAGCGTGAAAGGTGCGTGTTTCCGAGA", 5, 3);
///
/// for (kmer, minimizer) in tokenizer {
///     // ... do what you want ...
/// }
pub struct Forward<'a, M>
where
    M: method::Method<u64>,
{
    kmer_mask: u64,
    seq: &'a [u8],
    pos: usize,
    kmer: u64,
    minimizer: M,
}

impl<'a, M> Forward<'a, M>
where
    M: method::Method<u64>,
{
    /// Create a new Forward on seq DNA kmer size is equal to k, minimizer size is equal to m
    pub fn new(seq: &'a [u8], k: u8, m: u8) -> Self {
        let kmer = unsafe { kmer::seq2bit(seq.get_unchecked(0..((k - 1) as usize))) };

        let mut minimizer = M::default();
        minimizer.init(k, m, kmer);

        Self {
            kmer_mask: (1 << (k * 2)) - 1,
            seq,
            pos: (k - 1) as usize,
            kmer,
            minimizer,
        }
    }
}

impl<'a, M> Iterator for Forward<'a, M>
where
    M: method::Method<u64>,
{
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn basic() {
        let seq = b"acgtgcgtgagaattgttcggtggaacaggcg";

        let token = Forward::<method::Random>::new(seq, 11, 7);

        let mut kmers = Vec::new();
        let mut canos = Vec::new();
        let mut minis = Vec::new();

        for (kmer, mini) in token {
            let cano = crate::kmer::canonical(kmer, 11);

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

        let fwd_token = Forward::<method::Random>::new(fwd, 11, 7);
        let mut fwd_canos = Vec::new();
        let mut fwd_minis = Vec::new();

        let rev_token = Forward::<method::Random>::new(rev, 11, 7);
        let mut rev_canos = Vec::new();
        let mut rev_minis = Vec::new();

        for (kmer, mini) in fwd_token {
            let cano = crate::kmer::canonical(kmer, 11);

            fwd_canos.push(cano);
            fwd_minis.push(mini);
        }

        for (kmer, mini) in rev_token {
            let cano = crate::kmer::canonical(kmer, 11);

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
