//! 2bit minimizer tokenizer

/* standard use */

/* crates use */

/* local use */
use crate::bytevec;
use crate::kmer;
use crate::tokenizer::minimizer::method;

/// An iterator that takes a DNA sequence and produces kmers (in the canonical direction and 2bit form) and the associated minimizer.
///
/// # Example
///
/// ```
/// use cocktail::tokenizer::minimizer::Canonical;
/// use cocktail::tokenizer::minimizer::method;
///
/// let tokenizer = Canonical::<method::Random, u64>::new(b"GTACTGTGCCCGTGTTACTTAGTAAGCGTGAAAGGTGCGTGTTTCCGAGA", 8, 7);
///
/// for kmer in tokenizer {
///     // ... do what you want ...
/// }
pub struct Canonical<'a, M, K>
where
    M: method::Method<K>,
{
    move_bit: u8,
    kmer_mask: u64,
    seq: &'a [u8],
    pos: usize,
    kmers: [K; 2],
    minimizer: M,
}

impl<'a, M> Canonical<'a, M, u64>
where
    M: method::Method<u64>,
{
    /// Create a new Canonical tokenizer on seq DNA, kmer size is equal to k
    pub fn new(seq: &'a [u8], k: u8, m: u8) -> Self {
        let forward = unsafe { kmer::seq2bit(seq.get_unchecked(0..((k - 1) as usize))) };

        let mut minimizer = M::default();
        minimizer.init(k, m, forward);

        Canonical {
            move_bit: (k - 1) * 2,
            kmer_mask: (1 << (k * 2)) - 1,
            seq,
            pos: (k - 1) as usize,
            kmers: [forward, kmer::revcomp(forward, k)],
            minimizer,
        }
    }
}

impl<'a, M> Iterator for Canonical<'a, M, u64>
where
    M: method::Method<u64>,
{
    type Item = (u64, u64);

    fn next(&mut self) -> Option<Self::Item> {
        if self.pos == self.seq.len() {
            None
        } else {
            unsafe {
                let nuc = kmer::nuc2bit(*self.seq.get_unchecked(self.pos));
                self.pos += 1;

                *self.kmers.get_unchecked_mut(0) =
                    ((self.kmers.get_unchecked(0) << 2) & self.kmer_mask) | nuc;
                *self.kmers.get_unchecked_mut(1) =
                    (*self.kmers.get_unchecked(1) >> 2) ^ ((nuc ^ 0b10) << self.move_bit);

                self.minimizer.add_kmer(*self.kmers.get_unchecked(0));

                if kmer::parity_even(*self.kmers.get_unchecked(0)) {
                    Some((*self.kmers.get_unchecked(0), self.minimizer.get_mini().0))
                } else {
                    Some((*self.kmers.get_unchecked(1), self.minimizer.get_mini().0))
                }
            }
        }
    }
}

impl<'a, M> Canonical<'a, M, Vec<u8>>
where
    M: method::Method<Vec<u8>>,
{
    /// Create a new Canonical tokenizer on seq DNA, kmer size is equal to k
    pub fn new(seq: &'a [u8], k: u8, m: u8) -> Self {
        let mut forward = unsafe { seq.get_unchecked(0..((k - 1) as usize)).to_vec() };
        forward.push(b'n');

        let mut minimizer = M::default();
        minimizer.init(k, m, forward.clone());

        forward.rotate_right(1);
        let reverse = bytevec::revcomp(&forward);

        Canonical {
            move_bit: (k - 1) * 2,
            kmer_mask: (1 << (k * 2)) - 1,
            seq,
            pos: (k - 1) as usize,
            kmers: [forward, reverse],
            minimizer,
        }
    }
}

impl<'a, M> Iterator for Canonical<'a, M, Vec<u8>>
where
    M: method::Method<Vec<u8>>,
{
    type Item = (Vec<u8>, u64);

    fn next(&mut self) -> Option<Self::Item> {
        if self.pos == self.seq.len() {
            None
        } else {
            unsafe {
                let nuc = *self.seq.get_unchecked(self.pos);
                self.pos += 1;

                let forward = self.kmers.get_unchecked_mut(0);
                forward.rotate_left(1);
                let end = forward.len() - 1;
                *forward.get_unchecked_mut(end) = nuc;

                let reverse = self.kmers.get_unchecked_mut(1);
                reverse.rotate_right(1);
                *reverse.get_unchecked_mut(0) = bytevec::comp(&nuc);

                self.minimizer
                    .add_kmer(self.kmers.get_unchecked(0).to_vec());

                if self.kmers.get_unchecked(0) < self.kmers.get_unchecked(1) {
                    Some((
                        self.kmers.get_unchecked(0).to_vec(),
                        self.minimizer.get_mini().0,
                    ))
                } else {
                    Some((
                        self.kmers.get_unchecked(1).to_vec(),
                        self.minimizer.get_mini().0,
                    ))
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    /* std use */

    /* 3rd party use */
    use biotest::Format;

    /* local use */
    use super::*;

    #[test]
    fn u64() {
        let mut rng = biotest::rand();
        let generator = biotest::Sequence::builder()
            .sequence_len(50)
            .build()
            .unwrap();
        let mut seq = vec![];
        generator.record(&mut seq, &mut rng).unwrap();

        let token = Canonical::<method::Random, u64>::new(&seq, 11, 7);

        let mut kmers = Vec::new();
        let mut minis = Vec::new();

        for (kmer, mini) in token {
            assert!(kmer::parity_even(kmer));
            kmers.push(kmer::kmer2seq(kmer, 11));
            minis.push(kmer::kmer2seq(mini, 7));
        }

        assert_eq!(
            kmers,
            vec![
                b"CGATTCATATA".to_vec(),
                b"GCGATTCATAT".to_vec(),
                b"CGCGATTCATA".to_vec(),
                b"ACGCGATTCAT".to_vec(),
                b"CACGCGATTCA".to_vec(),
                b"ACACGCGATTC".to_vec(),
                b"AATCGCGTGTT".to_vec(),
                b"ATCGCGTGTTA".to_vec(),
                b"TCGCGTGTTAG".to_vec(),
                b"CGCGTGTTAGT".to_vec(),
                b"GCGTGTTAGTT".to_vec(),
                b"CGTGTTAGTTA".to_vec(),
                b"TTAACTAACAC".to_vec(),
                b"CTTAACTAACA".to_vec(),
                b"GCTTAACTAAC".to_vec(),
                b"TTAGTTAAGCC".to_vec(),
                b"TGGCTTAACTA".to_vec(),
                b"GTGGCTTAACT".to_vec(),
                b"CGTGGCTTAAC".to_vec(),
                b"CCGTGGCTTAA".to_vec(),
                b"ACCGTGGCTTA".to_vec(),
                b"AAGCCACGGTA".to_vec(),
                b"AGCCACGGTAA".to_vec(),
                b"ATTACCGTGGC".to_vec(),
                b"CATTACCGTGG".to_vec(),
                b"GCATTACCGTG".to_vec(),
                b"AGCATTACCGT".to_vec(),
                b"CGGTAATGCTT".to_vec(),
                b"CAAGCATTACC".to_vec(),
                b"GTAATGCTTGT".to_vec(),
                b"TAATGCTTGTA".to_vec(),
                b"AATGCTTGTAC".to_vec(),
                b"ATGCTTGTACG".to_vec(),
                b"GCGTACAAGCA".to_vec(),
                b"GCTTGTACGCA".to_vec(),
                b"CTTGTACGCAG".to_vec(),
                b"CCTGCGTACAA".to_vec(),
                b"TGTACGCAGGA".to_vec(),
                b"GTACGCAGGAT".to_vec(),
                b"TACGCAGGATA".to_vec(),
            ]
        );

        assert_eq!(
            minis,
            vec![
                b"CGATTCA".to_vec(),
                b"CGATTCA".to_vec(),
                b"CGATTCA".to_vec(),
                b"CGATTCA".to_vec(),
                b"CGATTCA".to_vec(),
                b"TCGCGTG".to_vec(),
                b"TCGCGTG".to_vec(),
                b"TCGCGTG".to_vec(),
                b"TCGCGTG".to_vec(),
                b"CGTGTTA".to_vec(),
                b"GTTAGTT".to_vec(),
                b"GTTAGTT".to_vec(),
                b"GTTAGTT".to_vec(),
                b"GTTAGTT".to_vec(),
                b"GTTAGTT".to_vec(),
                b"AGTTAAG".to_vec(),
                b"AGTTAAG".to_vec(),
                b"AGTTAAG".to_vec(),
                b"CGTGGCT".to_vec(),
                b"CGTGGCT".to_vec(),
                b"CGTGGCT".to_vec(),
                b"CGTGGCT".to_vec(),
                b"CGTGGCT".to_vec(),
                b"TACCGTG".to_vec(),
                b"TACCGTG".to_vec(),
                b"TACCGTG".to_vec(),
                b"ATTACCG".to_vec(),
                b"AATGCTT".to_vec(),
                b"AATGCTT".to_vec(),
                b"AATGCTT".to_vec(),
                b"AATGCTT".to_vec(),
                b"AATGCTT".to_vec(),
                b"ACAAGCA".to_vec(),
                b"TGTACGC".to_vec(),
                b"TGTACGC".to_vec(),
                b"TGTACGC".to_vec(),
                b"TGTACGC".to_vec(),
                b"TGTACGC".to_vec(),
                b"GCAGGAT".to_vec(),
                b"GCAGGAT".to_vec(),
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

        let token = Canonical::<method::Random, Vec<u8>>::new(&seq, 11, 7);

        let mut kmers = Vec::new();
        let mut minis = Vec::new();

        for (kmer, mini) in token {
            kmers.push(kmer);
            minis.push(kmer::kmer2seq(mini, 7));
        }

        // taTATgAAtCGCgtGTTAGTTAagccAcggtAatGcTtgtaCgcAGgAta
        assert_eq!(
            kmers,
            vec![
                b"CGaTTcATAta".to_vec(),
                b"GCGaTTcATAt".to_vec(),
                b"TATgAAtCGCg".to_vec(),
                b"ATgAAtCGCgt".to_vec(),
                b"CacGCGaTTcA".to_vec(),
                b"ACacGCGaTTc".to_vec(),
                b"AACacGCGaTT".to_vec(),
                b"AtCGCgtGTTA".to_vec(),
                b"CTAACacGCGa".to_vec(),
                b"ACTAACacGCG".to_vec(),
                b"AACTAACacGC".to_vec(),
                b"CgtGTTAGTTA".to_vec(),
                b"gtGTTAGTTAa".to_vec(),
                b"ctTAACTAACa".to_vec(),
                b"GTTAGTTAagc".to_vec(),
                b"TTAGTTAagcc".to_vec(),
                b"TAGTTAagccA".to_vec(),
                b"AGTTAagccAc".to_vec(),
                b"GTTAagccAcg".to_vec(),
                b"TTAagccAcgg".to_vec(),
                b"TAagccAcggt".to_vec(),
                b"AagccAcggtA".to_vec(),
                b"agccAcggtAa".to_vec(),
                b"atTaccgTggc".to_vec(),
                b"CatTaccgTgg".to_vec(),
                b"cAcggtAatGc".to_vec(),
                b"AcggtAatGcT".to_vec(),
                b"aAgCatTaccg".to_vec(),
                b"caAgCatTacc".to_vec(),
                b"acaAgCatTac".to_vec(),
                b"tAatGcTtgta".to_vec(),
                b"AatGcTtgtaC".to_vec(),
                b"atGcTtgtaCg".to_vec(),
                b"gcGtacaAgCa".to_vec(),
                b"GcTtgtaCgcA".to_vec(),
                b"CTgcGtacaAg".to_vec(),
                b"TtgtaCgcAGg".to_vec(),
                b"TcCTgcGtaca".to_vec(),
                b"aTcCTgcGtac".to_vec(),
                b"taCgcAGgAta".to_vec(),
            ]
        );

        assert_eq!(
            minis,
            vec![
                b"CGATTCA".to_vec(),
                b"CGATTCA".to_vec(),
                b"CGATTCA".to_vec(),
                b"CGATTCA".to_vec(),
                b"CGATTCA".to_vec(),
                b"TCGCGTG".to_vec(),
                b"TCGCGTG".to_vec(),
                b"TCGCGTG".to_vec(),
                b"TCGCGTG".to_vec(),
                b"CGTGTTA".to_vec(),
                b"GTTAGTT".to_vec(),
                b"GTTAGTT".to_vec(),
                b"GTTAGTT".to_vec(),
                b"GTTAGTT".to_vec(),
                b"GTTAGTT".to_vec(),
                b"AGTTAAG".to_vec(),
                b"AGTTAAG".to_vec(),
                b"AGTTAAG".to_vec(),
                b"CGTGGCT".to_vec(),
                b"CGTGGCT".to_vec(),
                b"CGTGGCT".to_vec(),
                b"CGTGGCT".to_vec(),
                b"CGTGGCT".to_vec(),
                b"TACCGTG".to_vec(),
                b"TACCGTG".to_vec(),
                b"TACCGTG".to_vec(),
                b"ATTACCG".to_vec(),
                b"AATGCTT".to_vec(),
                b"AATGCTT".to_vec(),
                b"AATGCTT".to_vec(),
                b"AATGCTT".to_vec(),
                b"AATGCTT".to_vec(),
                b"ACAAGCA".to_vec(),
                b"TGTACGC".to_vec(),
                b"TGTACGC".to_vec(),
                b"TGTACGC".to_vec(),
                b"TGTACGC".to_vec(),
                b"TGTACGC".to_vec(),
                b"GCAGGAT".to_vec(),
                b"GCAGGAT".to_vec(),
            ]
        );
    }
}
