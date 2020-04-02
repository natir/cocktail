/*
Copyright (c) 2019 Pierre Marijon <pmarijon@mpi-inf.mpg.de>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
 */

#[inline(always)]
pub fn seq2bit(subseq: &[u8]) -> u64 {
    let mut kmer: u64 = 0;

    for n in subseq {
        kmer <<= 2;
        kmer |= nuc2bit(*n);
    }

    kmer
}

#[inline(always)]
pub fn nuc2bit(nuc: u8) -> u64 {
    (nuc as u64 >> 1) & 0b11
}

static mut KMER2SEQ_BUFFER: [u8; 31] = [0; 31];

#[inline(always)]
pub fn kmer2seq(mut kmer: u64, k: u8) -> String {
    for i in (0..k).rev() {
        unsafe {
            KMER2SEQ_BUFFER[i as usize] = bit2nuc(kmer & 0b11);
        }

        kmer >>= 2;
    }

    unsafe { String::from_utf8_unchecked((&KMER2SEQ_BUFFER[..k as usize]).to_vec()) }
}

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

#[inline(always)]
pub fn cannonical(kmer: u64, k: u8) -> u64 {
    if parity_even(kmer) {
        kmer
    } else {
        revcomp(kmer, k)
    }
}

#[inline(always)]
pub fn parity_even(kmer: u64) -> bool {
    kmer.count_ones() % 2 == 0
}

pub fn revcomp(kmer: u64, k: u8) -> u64 {
    rev(speed_comp(kmer), k)
}

fn speed_comp(kmer: u64) -> u64 {
    kmer ^ 0b1010_1010_1010_1010_1010_1010_1010_1010_1010_1010_1010_1010_1010_1010_1010_1010
}
#[inline(always)]
pub fn comp(kmer: u64, k: u8) -> u64 {
    speed_comp(kmer) & ((1 << (2 * k)) - 1)
}

#[inline(always)]
pub fn get_first_bit(kmer: u64) -> bool {
    kmer & 1 != 0
}

#[inline(always)]
pub fn remove_first_bit(kmer: u64) -> u64 {
    kmer >> 1
}

#[inline(always)]
pub fn hash(subseq: &[u8], k: u8) -> u64 {
    remove_first_bit(cannonical(seq2bit(subseq), k))
}

#[inline(always)]
pub fn rev(kmer: u64, k: u8) -> u64 {
    loop_rev(kmer, k)
}

mod lookup_table;

#[inline(always)]
pub fn loop_rev(mut kmer: u64, k: u8) -> u64 {
    let nb_bit = k * 2;
    let mut reverse: u64 = 0;

    let nb_block = 1 + nb_bit / 8; //odd kmer never fit in 8 block

    for i in 0..nb_block {
        reverse ^= (lookup_table::REVERSE_2_LOOKUP[(kmer & 255) as u8 as usize] as u64)
            << ((nb_block - i - 1) * 8);

        kmer >>= 8;
    }

    reverse >> (nb_block * 8 - nb_bit)
}

#[inline(always)]
pub fn unrool_rev(mut kmer: u64, k: u8) -> u64 {
    let nb_bit = k * 2;
    let mut reverse: u64 = 0;

    reverse ^= (lookup_table::REVERSE_2_LOOKUP[(kmer & 255) as u8 as usize] as u64) << 56;
    kmer >>= 8;

    reverse ^= (lookup_table::REVERSE_2_LOOKUP[(kmer & 255) as u8 as usize] as u64) << 48;
    kmer >>= 8;

    reverse ^= (lookup_table::REVERSE_2_LOOKUP[(kmer & 255) as u8 as usize] as u64) << 40;
    kmer >>= 8;

    reverse ^= (lookup_table::REVERSE_2_LOOKUP[(kmer & 255) as u8 as usize] as u64) << 32;
    kmer >>= 8;

    reverse ^= (lookup_table::REVERSE_2_LOOKUP[(kmer & 255) as u8 as usize] as u64) << 24;
    kmer >>= 8;

    reverse ^= (lookup_table::REVERSE_2_LOOKUP[(kmer & 255) as u8 as usize] as u64) << 16;
    kmer >>= 8;

    reverse ^= (lookup_table::REVERSE_2_LOOKUP[(kmer & 255) as u8 as usize] as u64) << 8;
    kmer >>= 8;

    reverse ^= lookup_table::REVERSE_2_LOOKUP[(kmer & 255) as u8 as usize] as u64;

    reverse >> (64 - nb_bit)
}

#[inline(always)]
pub fn get_kmer_space_size(k: u8) -> u64 {
    1 << (k * 2)
}

#[inline(always)]
pub fn get_hash_space_size(k: u8) -> u64 {
    1 << (k * 2 - 1)
}

#[cfg(test)]
mod test {
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
        assert_eq!(kmer2seq(0b1000111101, 5), "TAGGC");

        // 110101100 -> GCCTA
        assert_eq!(kmer2seq(0b1101011000, 5), "GCCTA");

        assert_eq!(
            kmer2seq(0b1101011000, 31),
            "AAAAAAAAAAAAAAAAAAAAAAAAAAGCCTA"
        );
    }

    #[test]
    fn cannonical_() {
        // TAGGC -> 1000111101 cannonical TAGGC -> 1000111101
        assert_eq!(cannonical(0b1000111101, 5), 0b1000111101);

        // GCCTA -> 1101011000 cannonical TAGGC -> 1000111101
        assert_eq!(cannonical(0b1101011000, 5), 0b1000111101);
    }

    #[test]
    fn parity_even_() {
        assert_eq!(parity_even(0b1111), true);
        assert_eq!(parity_even(0b1110), false);
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
