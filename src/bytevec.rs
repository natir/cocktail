//! Collection of function

/* standard use */

/* crates use */

/* project use */

#[inline(always)]
pub fn comp(nuc: &u8) -> u8 {
    if nuc & 2 != 0 {
        nuc ^ 4
    } else {
        nuc ^ 21
    }
}

#[inline(always)]
pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(comp).collect()
}

#[inline(always)]
pub fn canonical(forward: &[u8]) -> Vec<u8> {
    let reverse = revcomp(forward);
    if forward < &reverse {
        forward.to_vec()
    } else {
        reverse
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn uppercase() {
        assert_eq!(b"ACTGA".to_vec(), revcomp(b"TCAGT"))
    }

    #[test]
    fn lowercase() {
        assert_eq!(b"actga".to_vec(), revcomp(b"tcagt"))
    }

    #[test]
    fn cano() {
        assert_eq!(canonical(b"ACgTA"), b"ACgTA".to_vec());

        assert_eq!(canonical(b"GatCC"), b"GGatC".to_vec());
    }
}
