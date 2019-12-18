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

pub static REVERSE_2_LOOKUP: &[u8] = &[
    0, 64, 128, 192, 16, 80, 144, 208, 32, 96, 160, 224, 48, 112, 176, 240, 4, 68, 132, 196, 20,
    84, 148, 212, 36, 100, 164, 228, 52, 116, 180, 244, 8, 72, 136, 200, 24, 88, 152, 216, 40, 104,
    168, 232, 56, 120, 184, 248, 12, 76, 140, 204, 28, 92, 156, 220, 44, 108, 172, 236, 60, 124,
    188, 252, 1, 65, 129, 193, 17, 81, 145, 209, 33, 97, 161, 225, 49, 113, 177, 241, 5, 69, 133,
    197, 21, 85, 149, 213, 37, 101, 165, 229, 53, 117, 181, 245, 9, 73, 137, 201, 25, 89, 153, 217,
    41, 105, 169, 233, 57, 121, 185, 249, 13, 77, 141, 205, 29, 93, 157, 221, 45, 109, 173, 237,
    61, 125, 189, 253, 2, 66, 130, 194, 18, 82, 146, 210, 34, 98, 162, 226, 50, 114, 178, 242, 6,
    70, 134, 198, 22, 86, 150, 214, 38, 102, 166, 230, 54, 118, 182, 246, 10, 74, 138, 202, 26, 90,
    154, 218, 42, 106, 170, 234, 58, 122, 186, 250, 14, 78, 142, 206, 30, 94, 158, 222, 46, 110,
    174, 238, 62, 126, 190, 254, 3, 67, 131, 195, 19, 83, 147, 211, 35, 99, 163, 227, 51, 115, 179,
    243, 7, 71, 135, 199, 23, 87, 151, 215, 39, 103, 167, 231, 55, 119, 183, 247, 11, 75, 139, 203,
    27, 91, 155, 219, 43, 107, 171, 235, 59, 123, 187, 251, 15, 79, 143, 207, 31, 95, 159, 223, 47,
    111, 175, 239, 63, 127, 191, 255,
];
