//! This module provides iterator to produce kmer from DNA sequence

/* standard use */

/* crates use */

/* project use */

/* module declaration */
pub mod basic;
pub use basic::Tokenizer;

pub mod minimizer;
pub use minimizer::{MinimizerRing, TokenizerMini};

pub mod minimizer_bstr;
pub use minimizer_bstr::MiniBstr;

pub mod rle;
pub use rle::TokenizerRLE;

pub mod canonical;
pub use canonical::Canonical;
