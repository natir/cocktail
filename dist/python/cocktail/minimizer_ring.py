#!/usr/bin/env python

from ._lowlevel import ffi, lib

from .utils import RustObject

class MinimizerRing(RustObject):
    __dealloc_func__ = lib.cocktail_minimizerring_free

    def __init__(self, k: int, m: int, kmer: int):
        self._objptr = lib.cocktail_minimizerring_new(k, m, kmer)

    def populate_buffer(self, kmer: int):
        self._methodcall(lib.cocktail_minimizerring_populate_buffer, kmer)

    def add_kmer(self, kmer: int):
        self._methodcall(lib.cocktail_minimizerring_add_kmer, kmer)

    def get_mini(self) -> int:
        return self._methodcall(lib.cocktail_minimizerring_get_mini)
