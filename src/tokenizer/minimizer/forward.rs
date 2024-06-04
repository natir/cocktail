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
    /* crate use */
    use biotest::Format as _;

    /* project use */
    use super::*;
    use crate::bytevec;

    #[test]
    fn u64() {
        let mut rng = biotest::rand();
        let generator = biotest::Sequence::builder()
            .sequence_len(50)
            .build()
            .unwrap();
        let mut seq = vec![];
        generator.record(&mut seq, &mut rng).unwrap();

        let token = Forward::<method::Random, u64>::new(&seq, 11, 7);

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
                2239527, 569501, 2278007, 723422, 2893691, 3186158, 161722, 646888, 2587555,
                1961614, 3652154, 2025704, 3908512, 3051139, 3815949, 2680885, 2334932, 951121,
                3804487, 2635039, 2151550, 217592, 870368, 3481474, 1342987, 1177645, 516278,
                2065114, 4066155, 3681710, 2143928, 187105, 748423, 2993693, 3586164, 1761747,
                2852687, 3022140, 3699954, 2216904
            ]
        );

        assert_eq!(
            canos,
            [
                1877128, 3615010, 1952328, 488082, 1170596, 292649, 161722, 646888, 2587555,
                1961614, 3652154, 2025704, 2627601, 1705476, 3572097, 2680885, 3106840, 3922438,
                2029185, 1555872, 388968, 217592, 870368, 661437, 1213935, 3449211, 862302,
                2065114, 1102469, 3681710, 2143928, 187105, 748423, 3638324, 3586164, 1761747,
                1498640, 3022140, 3699954, 2216904
            ]
        );

        assert_eq!(
            minis,
            [
                7332, 7332, 7332, 7332, 7332, 10107, 10107, 10107, 10107, 7912, 14906, 14906,
                14906, 14906, 14906, 3715, 3715, 3715, 7926, 7926, 7926, 7926, 7926, 8571, 8571,
                8571, 2583, 730, 730, 730, 730, 730, 1076, 11805, 11805, 11805, 11805, 11805,
                13554, 13554
            ]
        );
    }

    #[test]
    fn bytevec() {
        let mut rng = biotest::rand();
        let generator = biotest::Sequence::builder()
            .sequence_len(50)
            .build()
            .unwrap();
        let mut seq = vec![];
        generator.record(&mut seq, &mut rng).unwrap();

        let token = Forward::<method::Random, Vec<u8>>::new(&seq, 11, 7);

        let mut kmers = Vec::new();
        let mut minis = Vec::new();

        for (kmer, mini) in token {
            kmers.push(kmer);
            minis.push(mini);
        }

        assert_eq!(
            kmers,
            [
                b"taTATgAAtCG".to_vec(),
                b"aTATgAAtCGC".to_vec(),
                b"TATgAAtCGCg".to_vec(),
                b"ATgAAtCGCgt".to_vec(),
                b"TgAAtCGCgtG".to_vec(),
                b"gAAtCGCgtGT".to_vec(),
                b"AAtCGCgtGTT".to_vec(),
                b"AtCGCgtGTTA".to_vec(),
                b"tCGCgtGTTAG".to_vec(),
                b"CGCgtGTTAGT".to_vec(),
                b"GCgtGTTAGTT".to_vec(),
                b"CgtGTTAGTTA".to_vec(),
                b"gtGTTAGTTAa".to_vec(),
                b"tGTTAGTTAag".to_vec(),
                b"GTTAGTTAagc".to_vec(),
                b"TTAGTTAagcc".to_vec(),
                b"TAGTTAagccA".to_vec(),
                b"AGTTAagccAc".to_vec(),
                b"GTTAagccAcg".to_vec(),
                b"TTAagccAcgg".to_vec(),
                b"TAagccAcggt".to_vec(),
                b"AagccAcggtA".to_vec(),
                b"agccAcggtAa".to_vec(),
                b"gccAcggtAat".to_vec(),
                b"ccAcggtAatG".to_vec(),
                b"cAcggtAatGc".to_vec(),
                b"AcggtAatGcT".to_vec(),
                b"cggtAatGcTt".to_vec(),
                b"ggtAatGcTtg".to_vec(),
                b"gtAatGcTtgt".to_vec(),
                b"tAatGcTtgta".to_vec(),
                b"AatGcTtgtaC".to_vec(),
                b"atGcTtgtaCg".to_vec(),
                b"tGcTtgtaCgc".to_vec(),
                b"GcTtgtaCgcA".to_vec(),
                b"cTtgtaCgcAG".to_vec(),
                b"TtgtaCgcAGg".to_vec(),
                b"tgtaCgcAGgA".to_vec(),
                b"gtaCgcAGgAt".to_vec(),
                b"taCgcAGgAta".to_vec(),
            ]
        );

        assert_eq!(
            minis,
            [
                7332, 7332, 7332, 7332, 7332, 10107, 10107, 10107, 10107, 7912, 14906, 14906,
                14906, 14906, 14906, 3715, 3715, 3715, 7926, 7926, 7926, 7926, 7926, 8571, 8571,
                8571, 2583, 730, 730, 730, 730, 730, 1076, 11805, 11805, 11805, 11805, 11805,
                13554, 13554
            ]
        );
    }

    #[test]
    fn same_in_each_strand() {
        let mut rng = biotest::rand();
        let generator = biotest::Sequence::builder()
            .sequence_len(30)
            .build()
            .unwrap();
        let mut fwd = vec![];
        generator.record(&mut fwd, &mut rng).unwrap();
        let rev = bytevec::revcomp(&fwd);

        let fwd_token = Forward::<method::Random, u64>::new(&fwd, 11, 7);
        let mut fwd_canos = Vec::new();
        let mut fwd_minis = Vec::new();

        let rev_token = Forward::<method::Random, u64>::new(&rev, 11, 7);
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
