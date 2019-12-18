/*
Copyright (c) 2019 Pierre Marijon <pmarijon@mpi-inf.mpg.de>

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

use bv;

fn get_k(filesize: u64) -> u8 {
    (((filesize as f64).log2() as u8) / 2) + 2
}

pub fn read_solidity_bitfield<R>(mut reader: R, filesize: u64) -> (u8, bv::BitVec<u8>)
where
    R: std::io::BufRead,
{
    let k = get_k(filesize);

    let mut data: Vec<u8> = vec![0u8; (crate::kmer::get_kmer_space_size(k) / 8) as usize];
    reader
        .read_exact(&mut data)
        .expect("Error durring reading of data");

    (k, bv::BitVec::from_bits(&data))
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn found_k_value() {
        assert_eq!(get_k(4), 3);
        assert_eq!(get_k(64), 5);
    }

    #[test]
    fn read_bitfield() {
        let input = std::io::Cursor::new([1, 1, 1, 1]);

        let (k, bitfield) = read_solidity_bitfield(std::io::BufReader::new(input), 4);

        assert_eq!(k, 3);

        assert_eq!(
            bitfield.as_slice(),
            [
                true, false, false, false, false, false, false, false, true, false, false, false,
                false, false, false, false, true, false, false, false, false, false, false, false,
                true, false, false, false, false, false, false, false
            ]
        )
    }
}
