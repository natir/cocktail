/*
Copyright (c) 2020 Pierre Marijon <pierre.marijon@hhu.de>

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

//! Found the first local minimum of a series of u128
pub fn found_first_local_min(data: Box<[u128]>) -> Option<usize> {
    for (i, d) in data.windows(2).enumerate() {
        if d[1] > d[0] {
            return Some(i);
        }
    }

    None
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn found_abundance() {
        let data = vec![
            992273316, 64106898, 6792586, 1065818, 220444, 62400, 36748, 54062, 100806, 178868,
            287058, 424184, 568742, 705680, 805332, 871544, 874546, 827252, 744428, 636722, 523488,
            418036, 320506, 237956, 170642, 118046, 77290, 48320, 30500, 21096, 15632, 12758,
            11838, 10888, 10402, 9872, 9018, 7960, 7236, 6304, 5276, 4524, 3714, 3056, 2628, 2018,
            1578, 1256, 1036, 906, 708, 716, 592, 476, 540, 520, 446, 388, 316, 264, 258, 200, 230,
            172, 164, 184, 154, 162, 126, 124, 126, 156, 152, 98, 116, 108, 134, 116, 88, 124, 96,
            94, 96, 72, 52, 56, 68, 50, 54, 66, 54, 28, 44, 48, 30, 42, 48, 32, 38, 34, 44, 30, 32,
            28, 18, 34, 20, 28, 26, 28, 28, 32, 22, 16, 10, 26, 8, 26, 14, 14, 30, 6, 32, 38, 26,
            26, 16, 30, 20, 38, 20, 22, 22, 28, 14, 16, 20, 20, 20, 10, 12, 14, 12, 10, 18, 16, 16,
            12, 18, 2, 14, 6, 12, 8, 0, 6, 2, 4, 2, 0, 0, 2, 4, 2, 2, 6, 6, 0, 0, 2, 0, 2, 4, 0, 2,
            2, 6, 2, 0, 0, 0, 2, 2, 2, 2, 2, 0, 2, 2, 0, 2, 2, 0, 2, 2, 2, 0, 0, 2, 4, 2, 0, 2, 0,
            2, 2, 2, 0, 2, 2, 2, 2, 2, 0, 0, 0, 2, 0, 0, 2, 2, 2, 2, 4, 0, 2, 4, 4, 0, 2, 0, 0, 2,
            2, 0, 0, 0, 0, 0, 4, 2, 0, 2, 0, 0, 0, 2, 0, 4, 2, 0, 4, 2, 0, 0, 284,
        ];

        assert_eq!(found_first_local_min(data.into_boxed_slice()), Some(6));
    }

    #[test]
    #[should_panic]
    fn failled_found_abundance() {
        let data = vec![1, 1, 1, 1, 1, 1];

        found_first_local_min(data.into_boxed_slice()).unwrap();
    }
}
