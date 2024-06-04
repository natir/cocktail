//! A set of function to convert small sequence (less than 32 nucleotide) in 2 bit representation.
//! - A or a -> 00
//! - C or c -> 01
//! - T or t -> 10
//! - G or g -> 11
//!
//! We use the second and thrid bit of each value provide, if you provide no ACTG value this function silently convert to A, C, T or G, for exemple N or n is convert in G.
//!
//! With this coding and if kmer size is odd, if the popcount of forward is odd the popcount of reverse is even. In this library if a kmer have even popcount is the canonical kmer.
//!
//! If we work only with canonical kmer, we can remove one bit at any extremity. To reconstruct lost bit, if result have even popcount we add a 0, if it's ood we add 1.
//!
//! This 2bit coding is inspired by https://cs.stackexchange.com/questions/82644/compact-mapping-from-an-involuted-set

/* standard use */

/* crates use */

/* project use */

/// Convert a sequence in 2 bit representation if suseq is larger than 32 only the last 32 nuc is store
#[inline(always)]
pub fn seq2bit(subseq: &[u8]) -> u64 {
    let mut kmer: u64 = 0;

    for n in subseq {
        kmer <<= 2;
        kmer |= nuc2bit(*n);
    }

    kmer
}

/// Convert a nucleotide in 2bit representation
#[inline(always)]
pub fn nuc2bit(nuc: u8) -> u64 {
    (nuc as u64 >> 1) & 0b11
}

/// Convert a 2 bit repersentation in String
#[inline(always)]
pub fn kmer2seq(mut kmer: u64, k: u8) -> Vec<u8> {
    let mut buffer: [u8; 31] = [0; 31];

    for i in (0..k).rev() {
        buffer[i as usize] = bit2nuc(kmer & 0b11);

        kmer >>= 2;
    }

    buffer[..k as usize].to_vec()
}

/// Convert the 2bit representation of a nucleotide in nucleotide
#[inline(always)]
pub fn bit2nuc(bit: u64) -> u8 {
    match bit {
        0 => b'A',
        1 => b'C',
        2 => b'T',
        3 => b'G',
        _ => b'G',
    }
}

/// Take a kmer and return the canonical form
#[inline(always)]
pub fn canonical(kmer: u64, k: u8) -> u64 {
    if parity_even(kmer) {
        kmer
    } else {
        revcomp(kmer, k)
    }
}

/// Return true if the kmer parity is even
#[inline(always)]
pub fn parity_even(kmer: u64) -> bool {
    kmer.count_ones() % 2 == 0
}

/// Return the reverse complement of kmer
#[inline(always)]
pub fn revcomp(kmer: u64, k: u8) -> u64 {
    rev(speed_comp(kmer), k)
}

#[inline(always)]
fn speed_comp(kmer: u64) -> u64 {
    kmer ^ 0b1010_1010_1010_1010_1010_1010_1010_1010_1010_1010_1010_1010_1010_1010_1010_1010
}

/// Return the complement of kmer
#[inline(always)]
pub fn comp(kmer: u64, k: u8) -> u64 {
    speed_comp(kmer) & ((1 << (2 * k)) - 1)
}

/// Return true if the right bit of kmer is 1
#[inline(always)]
pub fn get_first_bit(kmer: u64) -> bool {
    kmer & 1 != 0
}

/// Return the kmer without the rightest bit of kmer
#[inline(always)]
pub fn remove_first_bit(kmer: u64) -> u64 {
    kmer >> 1
}

/// Take a subseq and return the canonical kmer with out the rightest bit
#[inline(always)]
pub fn hash(subseq: &[u8], k: u8) -> u64 {
    remove_first_bit(canonical(seq2bit(subseq), k))
}

/// Return the reverse of kmer
#[inline(always)]
pub fn rev(mut kmer: u64, k: u8) -> u64 {
    // Thank to needtail people ! :)
    kmer = (kmer >> 2 & 0x3333_3333_3333_3333) | (kmer & 0x3333_3333_3333_3333) << 2;
    kmer = (kmer >> 4 & 0x0F0F_0F0F_0F0F_0F0F) | (kmer & 0x0F0F_0F0F_0F0F_0F0F) << 4;
    kmer = (kmer >> 8 & 0x00FF_00FF_00FF_00FF) | (kmer & 0x00FF_00FF_00FF_00FF) << 8;
    kmer = (kmer >> 16 & 0x0000_FFFF_0000_FFFF) | (kmer & 0x0000_FFFF_0000_FFFF) << 16;
    kmer = (kmer >> 32 & 0x0000_0000_FFFF_FFFF) | (kmer & 0x0000_0000_FFFF_FFFF) << 32;

    kmer >> (64 - k * 2)
}

/// Return the cardinality of canonical kmer set for a given kmer size
#[inline(always)]
pub fn get_kmer_space_size(k: u8) -> u64 {
    1 << (k * 2)
}

/// Return the cardinality of canonical hash set for a given kmer size
#[inline(always)]
pub fn get_hash_space_size(k: u8) -> u64 {
    1 << (k * 2 - 1)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn seq2bit_() {
        // TAGGC -> 1000111101
        assert_eq!(seq2bit(b"TAGGC"), 0b1000111101);

        // GCCTA -> 110101100
        assert_eq!(seq2bit(b"GCCTA"), 0b1101011000);
    }

    #[test]
    fn bit2seq_() {
        // 1000111101 -> TAGGC
        assert_eq!(kmer2seq(0b1000111101, 5), b"TAGGC".to_vec());

        // 110101100 -> GCCTA
        assert_eq!(kmer2seq(0b1101011000, 5), b"GCCTA".to_vec());

        assert_eq!(
            kmer2seq(0b1101011000, 31),
            b"AAAAAAAAAAAAAAAAAAAAAAAAAAGCCTA".to_vec()
        );
    }

    #[test]
    fn canonical_() {
        // TAGGC -> 1000111101 canonical TAGGC -> 1000111101
        assert_eq!(canonical(0b1000111101, 5), 0b1000111101);

        // GCCTA -> 1101011000 canonical TAGGC -> 1000111101
        assert_eq!(canonical(0b1101011000, 5), 0b1000111101);
    }

    #[test]
    fn parity_even_() {
        assert!(parity_even(0b1111));
        assert!(!parity_even(0b1110));
    }

    #[test]
    fn revcomp_() {
        // TAGGC -> 1000111101 revcomp GCCTA -> 1101011000
        assert_eq!(0b1000111101, revcomp(0b1101011000, 5))
    }

    #[test]
    fn comp_() {
        // TAGGC -> 1000111101 comp 0001001011
        assert_eq!(comp(0b1000111101, 5), 0b0010010111);
    }

    #[test]
    fn rev_() {
        // TAGGC -> 1000111101 rev CGGAT -> 0111110010
        let var = 0b1000111101;

        assert_eq!(498, rev(var, 5));
    }

    #[test]
    fn hash_() {
        // TAGGC -> 100011110
        assert_eq!(hash(b"TAGGC", 5), 0b100011110);

        // GCCTA -> 110101100
        assert_eq!(hash(b"GCCTA", 5), 0b100011110);
    }

    #[test]
    fn kmer_space_size() {
        assert_eq!(get_kmer_space_size(5), 1024);
        assert_eq!(get_kmer_space_size(15), 1073741824);
    }
}
