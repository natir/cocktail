/*
Copyright (c) 2020 Pierre Marijon <pmarijon@mpi-inf.mpg.de>

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

//! This module contain function for binding a C and Python binding is avaible

use crate::kmer;

/// Binding for [kmer::seq2bit] in Python the name is seq2bit and the parameter len isn't present
#[no_mangle]
pub extern fn cocktail_seq2bit(c_subseq: *const std::os::raw::c_char, len: usize) -> u64 {
    let subseq = unsafe { std::slice::from_raw_parts(c_subseq as *const u8, len) };

    kmer::seq2bit(subseq)
}

/// Binding for [kmer::nuc2bit] in Python the name is nuc2bit
#[no_mangle]
pub extern fn cocktail_nuc2bit(nuc: u8) -> u64 {
    kmer::nuc2bit(nuc)
}

/// Binding for [kmer::kmer2seq] in Python the name is kmer2seq
#[no_mangle]
pub extern fn cocktail_kmer2seq(kmer: u64, k: u8) -> *const std::os::raw::c_char {
    let c_string = kmer::kmer2seq(kmer, k);

    unsafe { std::ffi::CString::from_vec_unchecked(c_string.into_bytes()).into_raw() }
}

/// Binding for [kmer::bit2nuc] in Python the name is bit2nuc
#[no_mangle]
pub extern fn cocktail_bit2nuc(bit: u64) -> u8 {
    kmer::bit2nuc(bit)
}

/// Binding for [kmer::cannonical] in Python the name is cannonical
#[no_mangle]
pub extern fn cocktail_cannonical(kmer: u64, k: u8) -> u64 {
    kmer::cannonical(kmer, k)
}

/// Binding for [kmer::parity_even] in Python the name is parity_even
#[no_mangle]
pub extern fn cocktail_parity_even(kmer: u64) -> bool {
    kmer::parity_even(kmer)
}

/// Binding for [kmer::revcomp] in Python the name is revcomp
#[no_mangle]
pub extern fn cocktail_revcomp(kmer: u64, k: u8) -> u64 {
    kmer::revcomp(kmer, k)
}

/// Binding for [kmer::comp] in Python the name is comp
#[no_mangle]
pub extern fn cocktail_comp(kmer: u64, k: u8) -> u64 {
    kmer::comp(kmer, k)
}

/// Binding for [kmer::get_first_bit] in Python the name is get_first_bit
#[no_mangle]
pub extern fn cocktail_get_first_bit(kmer: u64) -> bool {
    kmer::get_first_bit(kmer)
}

/// Binding for [kmer::remove_first_bit] in Python the name is remove_first_bit
#[no_mangle]
pub extern fn cocktail_remove_first_bit(kmer: u64) -> u64 {
    kmer::remove_first_bit(kmer)
}

/// Binding for [kmer::hash] in Python the name is hash 
#[no_mangle]
pub extern fn cocktail_hash(c_subseq: *const std::os::raw::c_char, k: u8) -> u64 {
    let subseq = unsafe { std::slice::from_raw_parts(c_subseq as *const u8, k as usize) };

    kmer::hash(subseq, k)
}

/// Binding for [kmer::rev] in Python the name is rev
#[no_mangle]
pub extern fn cocktail_rev(kmer: u64, k: u8) -> u64 {
    kmer::rev(kmer, k)
}

/// Binding for [kmer::get_kmer_space_size] in Python the name is get_kmer_space_size
#[no_mangle]
pub extern fn cocktail_get_kmer_space_size(k: u8) -> u64 {
    kmer::get_kmer_space_size(k)
}

/// Binding for [kmer::get_hash_space_size] in Python the name is get_hash_space_size
#[no_mangle]
pub extern fn cocktail_get_hash_space_size(k: u8) -> u64 {
    kmer::get_hash_space_size(k)
}
