//! Select minimizer with random scoring method

/* standard use */
use std::ops::BitXor;

/* crates use */

/* local use */
use crate::kmer;
use crate::tokenizer::minimizer::method;

/// A struct to get minimizer of sucessive kmer
///
/// At initialization all subkmer with weight is compute and store in a ring buffer.
/// When the next kmer is add only the new subkmer and is weight is compute.
/// If the new subkmer erase the previous minimizer but is score isn't lower than previous minimizer, the ring buffer is scanned completely to find the new minimizer.
#[derive(core::default::Default)]
pub struct Random {
    ring_buffer: Box<[(u64, u64)]>,
    current: usize,
    minimizer: usize,
    mask: u64,
    k: u8,
    m: u8,
}

impl Random {
    fn update_minimizer(&mut self) {
        let mut index: usize = 0;
        let mut min = u64::max_value();

        for (i, elt) in self.ring_buffer.iter().enumerate() {
            if elt.1 < min {
                index = i;
                min = elt.1;
            }
        }

        self.minimizer = index;
    }

    fn get_score(x: u64) -> u64 {
        x.rotate_left(5)
            .bitxor(x)
            .wrapping_mul(0x517c_c1b7_2722_0a95)
    }
}

impl method::Method<u64> for Random {
    /// Create a MinimizerRing, with kmer size equale to `k`, subkmer size equale to `m` and init ring buffer with `kmer`
    fn init(&mut self, k: u8, m: u8, mut kmer: u64) {
        self.ring_buffer = vec![(0, 0); (k - m + 1) as usize].into_boxed_slice();
        self.current = 0;
        self.minimizer = 0;
        self.mask = (1 << (m * 2)) - 1;
        self.k = k;
        self.m = m;

        // Populate buffer
        let mut score = u64::max_value();
        let max_len = (self.k - self.m + 1) as usize;

        for i in 0..max_len {
            let rb_index = max_len - i - 1;

            let mini = kmer::canonical(kmer & self.mask, self.m);

            let local_score = Random::get_score(mini);
            self.ring_buffer[rb_index] = (mini, local_score);

            if local_score < score {
                score = local_score;
                self.minimizer = rb_index;
            }

            kmer >>= 2;
        }

        self.current = 0;
    }

    /// Add the next kmer
    fn add_kmer(&mut self, kmer: u64) {
        let minimizer = kmer::canonical(kmer & self.mask, self.m);
        let score = Random::get_score(minimizer);

        let previous_mini = <Random as method::Method<u64>>::get_mini(self);
        self.ring_buffer[self.current] = (minimizer, score);

        if score < previous_mini.1 {
            self.minimizer = self.current;
        } else if self.current == self.minimizer {
            self.update_minimizer();
        }

        self.current = (self.current + 1) % self.ring_buffer.len();
    }

    /// Get a pair of value first one is the minimizer second one is his score
    fn get_mini(&self) -> (u64, u64) {
        self.ring_buffer[self.minimizer]
    }
}

impl method::Method<Vec<u8>> for Random {
    /// Create a MinimizerRing, with kmer size equale to `k`, subkmer size equale to `m` and init ring buffer with `kmer`
    fn init(&mut self, k: u8, m: u8, kmer: Vec<u8>) {
        self.ring_buffer = vec![(0, 0); (k - m + 1) as usize].into_boxed_slice();
        self.current = 0;
        self.minimizer = 0;
        self.mask = (1 << (m * 2)) - 1;
        self.k = k;
        self.m = m;

        // Populate buffer
        let mut score = u64::max_value();
        let max_len = (self.k - self.m + 1) as usize;

        for i in 0..max_len - 1 {
            let rb_index = i + 1;
            let mini = kmer::canonical(kmer::seq2bit(&kmer[i..i + self.m as usize]), self.m);

            let local_score = Random::get_score(mini);
            self.ring_buffer[rb_index] = (mini, local_score);

            if local_score < score {
                score = local_score;
                self.minimizer = rb_index;
            }
        }

        self.current = 0;
    }

    /// Add the next kmer
    fn add_kmer(&mut self, kmer: Vec<u8>) {
        let minimizer = kmer::canonical(kmer::seq2bit(&kmer[(self.k - self.m) as usize..]), self.m);
        let score = Random::get_score(minimizer);

        let previous_mini = <Random as method::Method<Vec<u8>>>::get_mini(self);
        self.ring_buffer[self.current] = (minimizer, score);

        if score < previous_mini.1 {
            self.minimizer = self.current;
        } else if self.current == self.minimizer {
            self.update_minimizer();
        }

        self.current = (self.current + 1) % self.ring_buffer.len();
    }

    /// Get a pair of value first one is the minimizer second one is his score
    fn get_mini(&self) -> (u64, u64) {
        self.ring_buffer[self.minimizer]
    }
}
