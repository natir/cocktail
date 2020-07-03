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

#include <stdio.h>

#include "cocktail.h"

int main(void) {

  printf("A in 2bit %lu\n", cocktail_nuc2bit('A'));
  printf("C in 2bit %lu\n", cocktail_nuc2bit('C'));
  printf("T in 2bit %lu\n", cocktail_nuc2bit('T'));
  printf("G in 2bit %lu\n\n", cocktail_nuc2bit('G'));
  
  printf("00 is %c\n", cocktail_bit2nuc(0));  
  printf("01 is %c\n", cocktail_bit2nuc(1));
  printf("10 is %c\n", cocktail_bit2nuc(2));
  printf("11 is %c\n\n", cocktail_bit2nuc(3));
  
  uint64_t kmer = cocktail_seq2bit("ACTGC", 5);
  printf("kmer ACTGT in 2bit %lu\n", kmer);
  printf("246 is kmer %s\n\n", cocktail_kmer2seq(246, 5));

  uint64_t cano = cocktail_canonical(kmer, 5);
  printf("kmer ACTGT parity %i\n", cocktail_parity_even(kmer));
  printf("kmer ACTGT revcomp parity %i\n", cocktail_parity_even(cano));
  printf("kmer ACTGT canonical %lu\n", cano);
  printf("kmer ACTGT revcomp %lu\n\n", cocktail_revcomp(kmer, 5));

  printf("kmer ACTGT canonical first bit %i\n", cocktail_get_first_bit(cano));
  printf("kmer ACTGT canonical without first bit %lu\n", cocktail_remove_first_bit(cano));
  printf("kmer ACTGT hash %lu\n\n", cocktail_hash("ACTGC", 5));

  printf("kmer space %lu\n", cocktail_get_kmer_space_size(5));
  printf("hash space %lu\n\n", cocktail_get_hash_space_size(5));

  MinimizerRing* miniring = cocktail_minimizerring_new(5, 3, cocktail_seq2bit("ACTGT", 5));

  printf("minimizer of ACTGT is %lu\n", cocktail_minimizerring_get_mini(miniring));
  
  cocktail_minimizerring_add_kmer(miniring, cocktail_seq2bit("CTGTA", 5));
  printf("minimizer of CTGTA is %lu\n", cocktail_minimizerring_get_mini(miniring));

  cocktail_minimizerring_add_kmer(miniring, cocktail_seq2bit("TGTAG", 5));
  printf("minimizer of TGTAG is %lu\n", cocktail_minimizerring_get_mini(miniring));

  cocktail_minimizerring_add_kmer(miniring, cocktail_seq2bit("GTAGA", 5));
  printf("minimizer of GTAGA is %lu\n", cocktail_minimizerring_get_mini(miniring));

  cocktail_minimizerring_add_kmer(miniring, cocktail_seq2bit("TAGAA", 5));
  printf("minimizer of TAGAA is %lu\n", cocktail_minimizerring_get_mini(miniring));
   
  cocktail_minimizerring_add_kmer(miniring, cocktail_seq2bit("AGAAA", 5));
  printf("minimizer of AGAAA is %lu\n", cocktail_minimizerring_get_mini(miniring));
  
  cocktail_minimizerring_free(miniring);
}
