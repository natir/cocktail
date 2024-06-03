//! A set of function to convert DNA sequence in Run Length Encoding (RLE) and work with it.
//! RLE is a small compression method where repetition of lettre are replace by the number of repetition and the lettre
//!
//! Exemple:
//! ACCTGGGAAT -> 1A2CT3G2A1T
//!
//! A RLE symbol are encode on a byte.
//! The two lowweight bits encode the nucletoide value, see [kmer](crate::kmer) module for more details.
//! - A or a -> 00
//! - C or c -> 01
//! - T or t -> 10
//! - G or g -> 11
//!
//! The six heavyweight bits encode the number of repetition plus one

/* standard use */

/* crates use */

/* project use */

/// Convert a sequence in rle representation
#[inline(always)]
pub fn seq2rle(seq: &[u8]) -> Box<[u8]> {
    let mut iterator = seq.iter();

    if let Some(mut prev) = iterator.next() {
        let mut rles = Vec::with_capacity(seq.len());

        let mut occurence = 0;

        for next in iterator {
            if prev != next || occurence == 63 {
                rles.push(build_rle(*prev, occurence));

                prev = next;
                occurence = 0;
            } else {
                occurence += 1;
            }
        }

        rles.push(build_rle(*prev, occurence));

        rles.into_boxed_slice()
    } else {
        Vec::new().into_boxed_slice()
    }
}

#[inline(always)]
fn build_rle(nuc: u8, occurence: u8) -> u8 {
    (crate::kmer::nuc2bit(nuc) as u8) | (occurence << 2)
}

/// Convert a rle sequence in 2 bit representation repetition is ignore. If sequence is larger than 32 only the last 32 nuc is store.
#[inline(always)]
pub fn rle2kmer(subrel: &[u8]) -> u64 {
    let mut kmer: u64 = 0;

    let mut iterator = subrel.iter();
    if let Some(mut prev) = iterator.next().map(|x| rle2bit(*x)) {
        kmer <<= 2;
        kmer |= prev;

        for n in iterator.map(|x| rle2bit(*x)) {
            if n != prev {
                kmer <<= 2;
                kmer |= n;

                prev = n;
            }
        }
    }

    kmer
}

/// Convert a rle in 2bit representation repetition is ignore.
#[inline(always)]
pub fn rle2bit(rle: u8) -> u64 {
    (rle & 0b11) as u64
}

/// Convert a rle in String.
#[inline(always)]
pub fn rle2seq(rles: &[u8]) -> String {
    let mut seq = Vec::with_capacity(rles.len());

    for rle in rles.iter() {
        for _ in 0..(((rle >> 2) + 1u8) as usize) {
            seq.push(crate::kmer::bit2nuc((rle & 0b11) as u64))
        }
    }

    unsafe { String::from_utf8_unchecked(seq) }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn build_rle_() {
        assert_eq!(build_rle(b'A', 0), 0b00000000);
        assert_eq!(build_rle(b'C', 0), 0b00000001);
        assert_eq!(build_rle(b'T', 0), 0b00000010);
        assert_eq!(build_rle(b'G', 0), 0b00000011);

        assert_eq!(build_rle(b'A', 18), 0b01001000);
        assert_eq!(build_rle(b'C', 21), 0b01010101);
        assert_eq!(build_rle(b'T', 15), 0b00111110);
        assert_eq!(build_rle(b'G', 62), 0b11111011);

        assert_eq!(build_rle(b'A', 63), 0b11111100);
        assert_eq!(build_rle(b'C', 63), 0b11111101);
        assert_eq!(build_rle(b'T', 63), 0b11111110);
        assert_eq!(build_rle(b'G', 63), 0b11111111);
    }

    #[test]
    fn seq2rle_() {
        assert_eq!(seq2rle(b"ACTG"), vec![0, 1, 2, 3].into_boxed_slice());
        assert_eq!(seq2rle(b"ACCTG"), vec![0, 5, 2, 3].into_boxed_slice());
        assert_eq!(seq2rle(b"ACCTG"), vec![0, 5, 2, 3].into_boxed_slice());
        assert_eq!(seq2rle(b"AACCTTGGG"), vec![4, 5, 6, 11].into_boxed_slice());
        assert_eq!(
            seq2rle(b"TAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAC"),
            vec![2, 252, 0, 1].into_boxed_slice()
        );
    }

    #[test]
    fn rle2kmer_() {
        assert_eq!(
            rle2kmer(&seq2rle(b"ACTGAG")),
            crate::kmer::seq2bit(b"ACTGAG")
        );
        assert_eq!(
            rle2kmer(&seq2rle(b"ACCTTTGAAAG")),
            crate::kmer::seq2bit(b"ACTGAG")
        );
        assert_eq!(
            rle2kmer(&seq2rle(
                b"TAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAC"
            )),
            crate::kmer::seq2bit(b"TAC")
        );
    }

    #[test]
    fn rle2bit_() {
        assert_eq!(rle2bit(seq2rle(b"AAA")[0]), 0);
        assert_eq!(rle2bit(seq2rle(b"C")[0]), 1);
        assert_eq!(rle2bit(seq2rle(b"TTTT")[0]), 2);
        assert_eq!(rle2bit(seq2rle(b"GGGGGGGGG")[0]), 3);
    }

    #[test]
    fn rle2seq_() {
        assert_eq!(rle2seq(&seq2rle(b"ACCGTTAGcATG")), "ACCGTTAGCATG");
        assert_eq!(
            rle2seq(&seq2rle(
                b"TAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAC"
            )),
            "TAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAC"
        );
    }
}
