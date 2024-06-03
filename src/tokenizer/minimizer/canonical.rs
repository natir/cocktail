//! 2bit minimizer tokenizer

/* standard use */

/* crates use */

/* local use */
use crate::kmer;
use crate::tokenizer::minimizer::method;

/// An iterator that takes a DNA sequence and produces kmers (in the canonical direction and 2bit form) and the associated minimizer.
///
/// # Example
///
/// ```
/// use cocktail::tokenizer::kmer::Canonical;
///
/// let tokenizer = Canonical::new(b"GTACTGTGCCCGTGTTACTTAGTAAGCGTGAAAGGTGCGTGTTTCCGAGA", 5);
///
/// for kmer in tokenizer {
///     // ... do what you want ...
/// }
pub struct Canonical<'a, M>
where
    M: method::Method<u64>,
{
    move_bit: u8,
    kmer_mask: u64,
    seq: &'a [u8],
    pos: usize,
    kmers: [u64; 2],
    minimizer: M,
}

impl<'a, M> Canonical<'a, M>
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

impl<'a, M> Iterator for Canonical<'a, M>
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

#[cfg(test)]
mod tests {
    /* std use */

    /* 3rd party use */
    use biotest::values::Generate;

    /* local use */
    use super::*;

    #[test]
    fn basic() {
        let mut rng = biotest::rand();
        let seq = biotest::values::Nucleotides::Dna
            .generate(&mut rng, 100)
            .unwrap();

        let token = Canonical::<method::Random>::new(&seq, 11, 7);

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
                b"ACGCAGGATAT".to_vec(),
                b"GATATCCTGCG".to_vec(),
                b"GCAGGATATCG".to_vec(),
                b"CAGGATATCGA".to_vec(),
                b"TTCGATATCCT".to_vec(),
                b"GGATATCGAAT".to_vec(),
                b"AATTCGATATC".to_vec(),
                b"TAATTCGATAT".to_vec(),
                b"TATCGAATTAT".to_vec(),
                b"TATAATTCGAT".to_vec(),
                b"CTATAATTCGA".to_vec(),
                b"CGAATTATAGA".to_vec(),
                b"GAATTATAGAT".to_vec(),
                b"AATTATAGATG".to_vec(),
                b"ATTATAGATGG".to_vec(),
                b"ACCATCTATAA".to_vec(),
                b"AACCATCTATA".to_vec(),
                b"ATAGATGGTTG".to_vec(),
                b"GCAACCATCTA".to_vec(),
                b"AGCAACCATCT".to_vec(),
                b"GATGGTTGCTC".to_vec(),
                b"ATGGTTGCTCA".to_vec(),
                b"ATGAGCAACCA".to_vec(),
                b"GGTTGCTCATG".to_vec(),
                b"ACATGAGCAAC".to_vec(),
                b"TTGCTCATGTC".to_vec(),
                b"TGCTCATGTCT".to_vec(),
                b"CAGACATGAGC".to_vec(),
                b"CTCATGTCTGC".to_vec(),
                b"TCATGTCTGCT".to_vec(),
                b"CAGCAGACATG".to_vec(),
                b"ATGTCTGCTGG".to_vec(),
                b"ACCAGCAGACA".to_vec(),
                b"GTCTGCTGGTA".to_vec(),
                b"GTACCAGCAGA".to_vec(),
                b"AGTACCAGCAG".to_vec(),
                b"TGCTGGTACTG".to_vec(),
                b"GCTGGTACTGT".to_vec(),
                b"CTGGTACTGTG".to_vec(),
                b"TGGTACTGTGC".to_vec(),
                b"TGCACAGTACC".to_vec(),
                b"TTGCACAGTAC".to_vec(),
                b"TTTGCACAGTA".to_vec(),
                b"ACTGTGCAAAA".to_vec(),
                b"CTGTGCAAAAG".to_vec(),
                b"CCTTTTGCACA".to_vec(),
                b"GTGCAAAAGGG".to_vec(),
                b"TGCAAAAGGGG".to_vec(),
                b"TCCCCTTTTGC".to_vec(),
                b"CTCCCCTTTTG".to_vec(),
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
                b"GCAGGAT".to_vec(),
                b"GCAGGAT".to_vec(),
                b"GCAGGAT".to_vec(),
                b"TCGATAT".to_vec(),
                b"TTCGATA".to_vec(),
                b"TTCGATA".to_vec(),
                b"TCGAATT".to_vec(),
                b"TCGAATT".to_vec(),
                b"TCGAATT".to_vec(),
                b"TCGAATT".to_vec(),
                b"TCGAATT".to_vec(),
                b"TCTATAA".to_vec(),
                b"ATCTATA".to_vec(),
                b"ATCTATA".to_vec(),
                b"TAGATGG".to_vec(),
                b"TAGATGG".to_vec(),
                b"TAGATGG".to_vec(),
                b"TAGATGG".to_vec(),
                b"TAGATGG".to_vec(),
                b"TGGTTGC".to_vec(),
                b"TGGTTGC".to_vec(),
                b"TGGTTGC".to_vec(),
                b"TGGTTGC".to_vec(),
                b"ATGAGCA".to_vec(),
                b"ATGAGCA".to_vec(),
                b"GACATGA".to_vec(),
                b"GACATGA".to_vec(),
                b"GACATGA".to_vec(),
                b"GCAGACA".to_vec(),
                b"GCAGACA".to_vec(),
                b"GCAGACA".to_vec(),
                b"GCAGACA".to_vec(),
                b"GCAGACA".to_vec(),
                b"CTGCTGG".to_vec(),
                b"CTGCTGG".to_vec(),
                b"CTGCTGG".to_vec(),
                b"TACCAGC".to_vec(),
                b"TACCAGC".to_vec(),
                b"TACTGTG".to_vec(),
                b"ACTGTGC".to_vec(),
                b"ACTGTGC".to_vec(),
                b"ACTGTGC".to_vec(),
                b"GTGCAAA".to_vec(),
                b"GTGCAAA".to_vec(),
                b"GTGCAAA".to_vec(),
                b"GTGCAAA".to_vec(),
                b"GTGCAAA".to_vec(),
                b"CCTTTTG".to_vec(),
                b"CCTTTTG".to_vec(),
                b"CCTTTTG".to_vec(),
            ]
        );
    }
}
