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
    use biotest::values::Generate;

    /* local use */
    use super::*;

    #[test]
    fn u64() {
        let mut rng = biotest::rand();
        let seq = biotest::values::Nucleotides::Dna
            .generate(&mut rng, 50)
            .unwrap();

        let token = Canonical::<method::Random, u64>::new(&seq, 11, 7);

        let mut kmers = Vec::new();
        let mut minis = Vec::new();

        for (kmer, mini) in token {
            assert!(kmer::parity_even(kmer));
            kmers.push(kmer::kmer2seq(kmer, 11));
            minis.push(kmer::kmer2seq(mini, 7));
        }

        println!("{}", String::from_utf8(seq).unwrap());
        // taTATgAAtCGCgtGTTAGTTAagccAcggtAatGcTtgtaCgcAGgAta
        // ATATACTTAGC
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
        let seq = biotest::values::Nucleotides::DnaUpper
            .generate(&mut rng, 50)
            .unwrap();

        let token = Canonical::<method::Random, Vec<u8>>::new(&seq, 11, 7);

        let mut kmers = Vec::new();
        let mut minis = Vec::new();

        for (kmer, mini) in token {
            kmers.push(kmer);
            minis.push(kmer::kmer2seq(mini, 7));
        }
        println!("{}", String::from_utf8(seq).unwrap());
        //GGTCTACACAAGGCCGACCAGTAGAGAAGTGCCGCAGTCAACTAGTCGGT
        assert_eq!(
            kmers,
            vec![
                b"GGTCTACACAA".to_vec(),
                b"CTTGTGTAGAC".to_vec(),
                b"CCTTGTGTAGA".to_vec(),
                b"CTACACAAGGC".to_vec(),
                b"GGCCTTGTGTA".to_vec(),
                b"ACACAAGGCCG".to_vec(),
                b"CACAAGGCCGA".to_vec(),
                b"ACAAGGCCGAC".to_vec(),
                b"CAAGGCCGACC".to_vec(),
                b"AAGGCCGACCA".to_vec(),
                b"AGGCCGACCAG".to_vec(),
                b"ACTGGTCGGCC".to_vec(),
                b"GCCGACCAGTA".to_vec(),
                b"CCGACCAGTAG".to_vec(),
                b"CGACCAGTAGA".to_vec(),
                b"CTCTACTGGTC".to_vec(),
                b"ACCAGTAGAGA".to_vec(),
                b"CCAGTAGAGAA".to_vec(),
                b"CAGTAGAGAAG".to_vec(),
                b"ACTTCTCTACT".to_vec(),
                b"CACTTCTCTAC".to_vec(),
                b"GCACTTCTCTA".to_vec(),
                b"AGAGAAGTGCC".to_vec(),
                b"CGGCACTTCTC".to_vec(),
                b"AGAAGTGCCGC".to_vec(),
                b"GAAGTGCCGCA".to_vec(),
                b"AAGTGCCGCAG".to_vec(),
                b"ACTGCGGCACT".to_vec(),
                b"GACTGCGGCAC".to_vec(),
                b"TGACTGCGGCA".to_vec(),
                b"GCCGCAGTCAA".to_vec(),
                b"CCGCAGTCAAC".to_vec(),
                b"AGTTGACTGCG".to_vec(),
                b"GCAGTCAACTA".to_vec(),
                b"CAGTCAACTAG".to_vec(),
                b"ACTAGTTGACT".to_vec(),
                b"GACTAGTTGAC".to_vec(),
                b"CGACTAGTTGA".to_vec(),
                b"CAACTAGTCGG".to_vec(),
                b"AACTAGTCGGT".to_vec(),
            ]
        );

        assert_eq!(
            minis,
            vec![
                b"GTGTAGA".to_vec(),
                b"GTGTAGA".to_vec(),
                b"GTGTAGA".to_vec(),
                b"TTGTGTA".to_vec(),
                b"TTGTGTA".to_vec(),
                b"ACACAAG".to_vec(),
                b"AAGGCCG".to_vec(),
                b"AAGGCCG".to_vec(),
                b"GCCGACC".to_vec(),
                b"CCGACCA".to_vec(),
                b"CCGACCA".to_vec(),
                b"CCGACCA".to_vec(),
                b"CCGACCA".to_vec(),
                b"CCGACCA".to_vec(),
                b"CTACTGG".to_vec(),
                b"CTACTGG".to_vec(),
                b"CTACTGG".to_vec(),
                b"CTACTGG".to_vec(),
                b"CAGTAGA".to_vec(),
                b"CTCTACT".to_vec(),
                b"CACTTCT".to_vec(),
                b"CACTTCT".to_vec(),
                b"CACTTCT".to_vec(),
                b"CACTTCT".to_vec(),
                b"CACTTCT".to_vec(),
                b"TGCCGCA".to_vec(),
                b"CTGCGGC".to_vec(),
                b"CTGCGGC".to_vec(),
                b"CTGCGGC".to_vec(),
                b"CTGCGGC".to_vec(),
                b"TTGACTG".to_vec(),
                b"TTGACTG".to_vec(),
                b"TTGACTG".to_vec(),
                b"TTGACTG".to_vec(),
                b"TTGACTG".to_vec(),
                b"TCAACTA".to_vec(),
                b"ACTAGTC".to_vec(),
                b"ACTAGTC".to_vec(),
                b"ACTAGTC".to_vec(),
                b"ACTAGTC".to_vec(),
            ]
        );
    }
}
