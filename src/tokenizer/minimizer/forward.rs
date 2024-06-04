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
/// let tokenizer = Forward::<method::Random, u64>::new(b"GTACTGTGCCCGTGTTACTTAGTAAGCGTGAAAGGTGCGTGTTTCCGAGA", 5, 3);
///
/// for (kmer, minimizer) in tokenizer {
///     // ... do what you want ...
/// }
pub struct Forward<'a, M, K>
where
    M: method::Method<K>,
{
    kmer_mask: u64,
    seq: &'a [u8],
    pos: usize,
    kmer: K,
    minimizer: M,
}

impl<'a, M> Forward<'a, M, u64>
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

impl<'a, M> Iterator for Forward<'a, M, u64>
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

impl<'a, M> Forward<'a, M, Vec<u8>>
where
    M: method::Method<Vec<u8>>,
{
    /// Create a new Forward on seq DNA kmer size is equal to k, minimizer size is equal to m
    pub fn new(seq: &'a [u8], k: u8, m: u8) -> Self {
        let mut kmer = unsafe { seq.get_unchecked(0..(k as usize - 1)).to_vec() };
        kmer.push(b'n');

        let mut minimizer = M::default();
        minimizer.init(k, m, kmer.clone());

        kmer.rotate_right(1);

        Self {
            kmer_mask: (1 << (k * 2)) - 1,
            seq,
            pos: (k - 1) as usize,
            kmer,
            minimizer,
        }
    }
}

impl<'a, M> Iterator for Forward<'a, M, Vec<u8>>
where
    M: method::Method<Vec<u8>>,
{
    type Item = (Vec<u8>, u64);

    fn next(&mut self) -> Option<Self::Item> {
        if self.pos == self.seq.len() {
            None
        } else {
            self.kmer.rotate_left(1);
            let end = self.kmer.len() - 1;
            unsafe { *self.kmer.get_unchecked_mut(end) = self.seq[self.pos] }

            self.minimizer.add_kmer(self.kmer.clone());

            self.pos += 1;

            Some((self.kmer.clone(), self.minimizer.get_mini().0))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn u64() {
        let seq = b"acgtgcgtgagaattgttcggtggaacaggcg";

        let token = Forward::<method::Random, u64>::new(seq, 11, 7);

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
    fn bytevec() {
        let seq = b"acgtgcgtgagaattgttcggtggaacaggcg";

        let token = Forward::<method::Random, Vec<u8>>::new(seq, 11, 7);

        let mut kmers = Vec::new();
        let mut minis = Vec::new();

        for (kmer, mini) in token {
            kmers.push(kmer);
            minis.push(mini);
        }

        assert_eq!(
            kmers,
            [
                b"acgtgcgtgag".to_vec(),
                b"cgtgcgtgaga".to_vec(),
                b"gtgcgtgagaa".to_vec(),
                b"tgcgtgagaat".to_vec(),
                b"gcgtgagaatt".to_vec(),
                b"cgtgagaattg".to_vec(),
                b"gtgagaattgt".to_vec(),
                b"tgagaattgtt".to_vec(),
                b"gagaattgttc".to_vec(),
                b"agaattgttcg".to_vec(),
                b"gaattgttcgg".to_vec(),
                b"aattgttcggt".to_vec(),
                b"attgttcggtg".to_vec(),
                b"ttgttcggtgg".to_vec(),
                b"tgttcggtgga".to_vec(),
                b"gttcggtggaa".to_vec(),
                b"ttcggtggaac".to_vec(),
                b"tcggtggaaca".to_vec(),
                b"cggtggaacag".to_vec(),
                b"ggtggaacagg".to_vec(),
                b"gtggaacaggc".to_vec(),
                b"tggaacaggcg".to_vec(),
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

        let fwd_token = Forward::<method::Random, u64>::new(fwd, 11, 7);
        let mut fwd_canos = Vec::new();
        let mut fwd_minis = Vec::new();

        let rev_token = Forward::<method::Random, u64>::new(rev, 11, 7);
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
