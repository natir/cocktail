//! C binding

/* standard use */

/* crates use */

/* project use */
use crate::tokenizer::MinimizerRing;

/// Create a cocktail MinimizerRing. See [MinimizerRing::new()]. In python MinimizerRing is  an object, this function is call in default constructor.
#[no_mangle]
pub extern "C" fn cocktail_minimizerring_new(k: u8, m: u8, kmer: u64) -> *mut MinimizerRing {
    Box::into_raw(Box::new(MinimizerRing::new(k, m, kmer)))
}

/// Free a cocktail minimizer ring
/// # Safety
/// Pointer is free don't use it after
#[no_mangle]
pub unsafe extern "C" fn cocktail_minimizerring_free(miniring: *mut MinimizerRing) {
    if miniring.is_null() {
        return;
    }

    drop(Box::from_raw(miniring));
}

/// Reset the ring buffer. See [MinimizerRing::populate_buffer()]. In python it's populate_buffer methode of MinimizerRing
#[no_mangle]
pub extern "C" fn cocktail_minimizerring_populate_buffer(miniring: &mut MinimizerRing, kmer: u64) {
    miniring.populate_buffer(kmer);
}

/// Add the next kmer. See [MinimizerRing::add_kmer()]. In python it's add_kmer methode of MinimizerRing
#[no_mangle]
pub extern "C" fn cocktail_minimizerring_add_kmer(miniring: &mut MinimizerRing, kmer: u64) {
    miniring.add_kmer(kmer);
}

/// Get the actual minimizer. See [MinimizerRing::get_mini()]. In python it's get_mini methode of MinimizerRing
#[no_mangle]
pub extern "C" fn cocktail_minimizerring_get_mini(miniring: &mut MinimizerRing) -> u64 {
    miniring.get_mini().0
}
