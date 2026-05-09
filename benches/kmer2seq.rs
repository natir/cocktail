static mut KMER2SEQ_BUFFER: [u8; 31] = [0; 31];

#[inline(always)]
pub fn static_buffer(mut kmer: u64, k: u8) -> Vec<u8> {
    for i in (0..k).rev() {
        unsafe {
            KMER2SEQ_BUFFER[i as usize] = cocktail::kmer::bit2nuc(kmer & 0b11);
        }

        kmer >>= 2;
    }

    unsafe { KMER2SEQ_BUFFER[..k as usize].to_vec() }
}

#[inline(always)]
pub fn local_buffer(mut kmer: u64, k: u8) -> Vec<u8> {
    let mut buffer: [u8; 31] = [0; 31];

    for i in (0..k).rev() {
        buffer[i as usize] = cocktail::kmer::bit2nuc(kmer & 0b11);

        kmer >>= 2;
    }

    buffer[..k as usize].to_vec()
}

#[inline(always)]
pub fn dyn_local_buffer(mut kmer: u64, k: u8) -> Vec<u8> {
    let mut buffer = vec![0; k as usize].into_boxed_slice();

    for i in (0..k).rev() {
        buffer[i as usize] = cocktail::kmer::bit2nuc(kmer & 0b11);

        kmer >>= 2;
    }

    buffer.to_vec()
}
