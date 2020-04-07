#!/usr/bin/env python

from ._lowlevel import ffi, lib

from .minimizer_ring import MinimizerRing

def seq2bit(subseq: str):
    return lib.cocktail_seq2bit(subseq.encode("utf-8"), len(subseq))

def nuc2bit(nuc: chr):
    return lib.cocktail_nuc2bit(ord(nuc))

def kmer2seq(kmer: int, k: int) -> str:
    return ffi.string(lib.cocktail_kmer2seq(kmer, k))

def bit2nuc(bit: int) -> chr:
    return chr(lib.cocktail_bit2nuc(bit))

def cannonical(kmer: int, k: int) -> int:
    return lib.cocktail_cannonical(kmer, k)

def parity_even(kmer: int) -> bool:
    return lib.cocktail_parity_even(kmer)

def revcomp(kmer: int, k: int) -> int:
    return lib.cocktail_revcomp(kmer, k)

def comp(kmer: int, k: int) -> int:
    return lib.cocktail_comp(kmer, k)

def get_first_bit(kmer: int) -> bool:
    return lib.cocktail_get_first_bit(kmer)

def remove_first_bit(kmer: int) -> int:
    return lib.cocktail_remove_first_bit(kmer)

def hash(subseq: str) -> int:
    return lib.cocktail_hash(subseq.encode("utf-8"), len(subseq))

def rev(kmer: int, k: int) -> int:
    return lib.cocktail_rev(kmer, k)

def get_kmer_space_size(k: int) -> int:
    return lib.cocktail_get_kmer_space_size(k)

def get_hash_space_size(k: int) -> int:
    return lib.cocktail_get_hash_space_size(k)
