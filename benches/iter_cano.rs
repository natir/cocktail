use cocktail::kmer;

pub struct Lexi<'a> {
    move_bit: u8,
    kmer_mask: u64,
    seq: &'a [u8],
    pos: usize,
    forward: u64,
    reverse: u64,
    forward_cano: bool,
}

impl<'a> Lexi<'a> {
    pub fn new(seq: &'a [u8], k: u8) -> Self {
        let forward = kmer::seq2bit(&seq[0..((k - 1) as usize)]);

        Lexi {
            move_bit: (k - 1) * 2,
            kmer_mask: (1 << (k * 2)) - 1,
            seq,
            pos: (k - 1) as usize,
            forward,
            reverse: kmer::revcomp(forward, k),
            forward_cano: kmer::parity_even(forward),
        }
    }
}

impl<'a> Iterator for Lexi<'a> {
    type Item = u64;

    fn next(&mut self) -> Option<Self::Item> {
        if self.pos == self.seq.len() {
            None
        } else {
            let nuc = kmer::nuc2bit(self.seq[self.pos]);
            self.pos += 1;

            self.forward = ((self.forward << 2) & self.kmer_mask) | nuc;
            self.reverse = (self.reverse >> 2) ^ ((nuc ^ 0b10) << self.move_bit);

            if self.forward < self.reverse {
                Some(self.forward)
            } else {
                Some(self.reverse)
            }
        }
    }
}
