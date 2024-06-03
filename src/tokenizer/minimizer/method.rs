//! Minimizer method
/* standard use */

/* crates use */

/* local use */

/* module declaration */
pub mod random;

/* reexport */
pub use random::Random;

/// Method
pub trait Method<T>: core::default::Default {
    fn init(&mut self, k: u8, m: u8, init_kmer: T);
    fn add_kmer(&mut self, kmer: T);
    fn get_mini(&self) -> (T, u64);
}
