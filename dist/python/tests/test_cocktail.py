import sys

sys.path.insert(0, sys.path.pop())

import cocktail

def test_nuc2bit():
    assert cocktail.nuc2bit('A') == 0
    assert cocktail.nuc2bit('C') == 1
    assert cocktail.nuc2bit('T') == 2
    assert cocktail.nuc2bit('G') == 3

def test_bit2nuc():
    assert cocktail.bit2nuc(0) == 'A'
    assert cocktail.bit2nuc(1) == 'C'
    assert cocktail.bit2nuc(2) == 'T'
    assert cocktail.bit2nuc(3) == 'G'

def test_kmer2seq_seq2kmer():
    assert cocktail.seq2bit("ACTGC") == 109
    assert cocktail.kmer2seq(246, 5) == b'AGGCT'

def test_cannonical():
    kmer = cocktail.seq2bit("ACTGC")
    assert cocktail.parity_even(kmer) == False

    cano = cocktail.cannonical(kmer, 5)
    assert cocktail.parity_even(cano) == True
    assert cano == 846

    assert cano == cocktail.revcomp(kmer, 5)

def test_space_size():
    assert cocktail.get_kmer_space_size(5) == 1024
    assert cocktail.get_hash_space_size(5) == 512
