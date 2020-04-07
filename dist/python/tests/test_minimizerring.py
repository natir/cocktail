import sys

sys.path.insert(0, sys.path.pop())

import cocktail

def test_minizerring():
    miniring = cocktail.MinimizerRing(5, 3, cocktail.seq2bit("ACTGT"))

    assert miniring.get_mini() == 6
    
    miniring.add_kmer(cocktail.seq2bit("CTGTA"))
    assert miniring.get_mini() == 46

    miniring.add_kmer(cocktail.seq2bit("TGTAG"))
    assert miniring.get_mini() == 24

    miniring.add_kmer(cocktail.seq2bit("GTAGA"))
    assert miniring.get_mini() == 12

    miniring.add_kmer(cocktail.seq2bit("TAGAA"))
    assert miniring.get_mini() == 12

    miniring.add_kmer(cocktail.seq2bit("AGAAA"))
    assert miniring.get_mini() == 0
