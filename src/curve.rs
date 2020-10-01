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

//! Based on Kmergenie we assume kmer spectrum is a mixture of Pareto law and some Gaussians law
//! Erroneous kmer follow Pareto law, Gaussians law represente true and repetitive kmer
//! We use this property to found the threshold to remove many Erroneous kmer and keep Many True kmer
pub mod threshold {

    //! The first local minimum match with the intersection of Pareto and Gaussians
    pub fn first_local_min(data: &[u64]) -> Option<u8> {
        for (i, d) in data.windows(2).enumerate() {
            if d[1] > d[0] {
                return Some(i as u8);
            }
        }

        None
    }

    pub fn rarefaction(data: &[u64], limit: f64) -> Option<u8> {
        //! More we remove kmer less we remove Erroneous kmer when remove less than n percent of total kmer
        let mut cumulative_kmer = 0;

        for (index, value) in data.iter().enumerate() {
            cumulative_kmer += index as u64 * value;
            if (*value as f64 / cumulative_kmer as f64) < limit {
                return Some(index as u8);
            }
        }

        None
    }
}

pub fn draw_hist<W>(data: &[u64], mut out: W, point: Option<usize>) -> std::io::Result<()>
where
    W: std::io::Write,
{
    //! Draw kmer spectrum in console
    let shape = get_shape();

    let factor = (*data.iter().max().unwrap() as f64).log(10.0) / shape.1 as f64;

    let normalized: Box<[f64]> = data
        .iter()
        .map(|y| {
            if *y == 0 {
                0.0
            } else {
                (*y as f64).log(10.0) / factor
            }
        })
        .collect();

    for h in (1..=shape.1).rev() {
        for w in 0..shape.0 {
            if normalized[w] >= h as f64 {
                let delta = normalized[w] - h as f64;
                if delta > 1.0 {
                    write!(out, "\u{258c}")?;
                } else if delta > 0.5 {
                    write!(out, "\u{2596}")?;
                } else {
                    write!(out, " ")?;
                }
            } else {
                write!(out, " ")?;
            }
        }

        writeln!(out)?;
    }

    let mut last_line = vec![b' '; shape.0];
    for x in (0..shape.0).step_by(5) {
        last_line[x] = b'5'
    }
    last_line[0] = b'0';
    if let Some(pos) = point {
        last_line[pos] = b'*';
    }

    writeln!(out, "{}", std::str::from_utf8(&last_line).unwrap())?;

    Ok(())
}

#[allow(dead_code)]
fn term_size() -> (usize, usize) {
    match term_size::dimensions() {
        Some((w, h)) => (w, h),
        None => (80, 24),
    }
}

#[cfg(test)]
fn get_shape() -> (usize, usize) {
    (256, 48)
}

#[cfg(not(test))]
fn get_shape() -> (usize, usize) {
    let term_size = term_size();

    (
        core::cmp::min(256, term_size.0),
        core::cmp::min(48, term_size.1),
    )
}

#[cfg(test)]
mod tests {

    use super::*;

    static DATA: [u64; 256] = [
        992273316, 64106898, 6792586, 1065818, 220444, 62400, 36748, 54062, 100806, 178868, 287058,
        424184, 568742, 705680, 805332, 871544, 874546, 827252, 744428, 636722, 523488, 418036,
        320506, 237956, 170642, 118046, 77290, 48320, 30500, 21096, 15632, 12758, 11838, 10888,
        10402, 9872, 9018, 7960, 7236, 6304, 5276, 4524, 3714, 3056, 2628, 2018, 1578, 1256, 1036,
        906, 708, 716, 592, 476, 540, 520, 446, 388, 316, 264, 258, 200, 230, 172, 164, 184, 154,
        162, 126, 124, 126, 156, 152, 98, 116, 108, 134, 116, 88, 124, 96, 94, 96, 72, 52, 56, 68,
        50, 54, 66, 54, 28, 44, 48, 30, 42, 48, 32, 38, 34, 44, 30, 32, 28, 18, 34, 20, 28, 26, 28,
        28, 32, 22, 16, 10, 26, 8, 26, 14, 14, 30, 6, 32, 38, 26, 26, 16, 30, 20, 38, 20, 22, 22,
        28, 14, 16, 20, 20, 20, 10, 12, 14, 12, 10, 18, 16, 16, 12, 18, 2, 14, 6, 12, 8, 0, 6, 2,
        4, 2, 0, 0, 2, 4, 2, 2, 6, 6, 0, 0, 2, 0, 2, 4, 0, 2, 2, 6, 2, 0, 0, 0, 2, 2, 2, 2, 2, 0,
        2, 2, 0, 2, 2, 0, 2, 2, 2, 0, 0, 2, 4, 2, 0, 2, 0, 2, 2, 2, 0, 2, 2, 2, 2, 2, 0, 0, 0, 2,
        0, 0, 2, 2, 2, 2, 4, 0, 2, 4, 4, 0, 2, 0, 0, 2, 2, 0, 0, 0, 0, 0, 4, 2, 0, 2, 0, 0, 0, 2,
        0, 4, 2, 0, 4, 2, 0, 0, 284,
    ];

    #[test]
    fn first_local_min() {
        assert_eq!(threshold::first_local_min(&DATA[..]), Some(6));
    }

    #[test]
    fn failled_first_local_min() {
        let tmp = (0..256).map(|_| 1).collect::<Box<[u64]>>();

        assert_eq!(threshold::first_local_min(&tmp), None);
    }

    #[test]
    fn rarefaction() {
        assert_eq!(threshold::rarefaction(&DATA[..], 0.01), Some(4));
    }

    #[test]
    fn failled_rarefaction() {
        let tmp = (0..256).map(|_| 1).collect::<Box<[u64]>>();

        assert_eq!(threshold::rarefaction(&tmp, 0.00001), None);
    }

    #[test]
    fn test_draw_hist() {
        let mut output = vec![];
        draw_hist(&DATA[..], &mut output, Some(6)).unwrap();

        let good_output = "                                                                                                                                                                                                                                                                
▖                                                                                                                                                                                                                                                               
▌                                                                                                                                                                                                                                                               
▌                                                                                                                                                                                                                                                               
▌                                                                                                                                                                                                                                                               
▌                                                                                                                                                                                                                                                               
▌                                                                                                                                                                                                                                                               
▌▖                                                                                                                                                                                                                                                              
▌▌                                                                                                                                                                                                                                                              
▌▌                                                                                                                                                                                                                                                              
▌▌                                                                                                                                                                                                                                                              
▌▌                                                                                                                                                                                                                                                              
▌▌                                                                                                                                                                                                                                                              
▌▌▌                                                                                                                                                                                                                                                             
▌▌▌                                                                                                                                                                                                                                                             
▌▌▌                                                                                                                                                                                                                                                             
▌▌▌                                                                                                                                                                                                                                                             
▌▌▌▌          ▖▖▖▖                                                                                                                                                                                                                                              
▌▌▌▌        ▖▌▌▌▌▌▌▖▖                                                                                                                                                                                                                                           
▌▌▌▌       ▌▌▌▌▌▌▌▌▌▌▖                                                                                                                                                                                                                                          
▌▌▌▌▖     ▌▌▌▌▌▌▌▌▌▌▌▌▌▖                                                                                                                                                                                                                                        
▌▌▌▌▌    ▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▖                                                                                                                                                                                                                                       
▌▌▌▌▌   ▖▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌                                                                                                                                                                                                                                      
▌▌▌▌▌▖  ▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌                                                                                                                                                                                                                                     
▌▌▌▌▌▌ ▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▖                                                                                                                                                                                                                                    
▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▖                                                                                                                                                                                                                                   
▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌                                                                                                                                                                                                                                  
▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▖▖▖                                                                                                                                                                                                                              
▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▖▖                                                                                                                                                                                                                         
▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▖▖                                                                                                                                                                                                                      
▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▖                                                                                                                                                                                                                    
▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▖                                                                                                                                                                                                                  
▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▖                                                                                                                                                                                                                
▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▖                                                                                                                                                                                                              
▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▖ ▖                                                                                                                                                                                                         
▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▖                                                                                                                                                                                                      
▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▖▖ ▖                                                                                                                                                                                                ▌
▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▖▖▌▖▖   ▖▖                                                                                                                                                                                      ▌
▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▖▌▖▌▌ ▌▖▖▖                                                                                                                                                                            ▌
▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▖  ▖  ▖                                                                                                                                                                     ▌
▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌ ▖▖ ▖▖   ▖                                                                                                                                                          ▌
▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▖▌▌▖▌▌▌▌▌▌▖▌▖ ▌ ▖▖▖▖▌   ▖ ▖  ▖ ▌▌▖▖ ▖ ▌   ▖                                                                                                                         ▌
▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▖▌▖▌▌▌▌▌▌  ▌ ▌  ▌ ▌▌▌▌ ▌▖▌▖▌▌▌  ▖▖▖     ▖   ▖                                                                                                          ▌
▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌ ▌ ▌▌▌▌ ▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌ ▖▌▖ ▌▌▌▖▌ ▌ ▖                                                                                                      ▌
▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▖▌▌▌▌ ▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌ ▌ ▌▖                                                                                                     ▌
▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌ ▌▌▌▌ ▌         ▌▌         ▌                                                                              ▌
▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌ ▌▌▌▌ ▌ ▌    ▌  ▌▌     ▌   ▌                      ▌                       ▌  ▌▌           ▌        ▌  ▌   ▌
▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▖▌▌▌▌ ▌▖▌▖  ▖▌▖▖▌▌  ▖ ▖▌ ▖▖▌▖   ▖▖▖▖▖ ▖▖ ▖▖ ▖▖▖  ▖▌▖ ▖ ▖▖▖ ▖▖▖▖▖   ▖  ▖▖▖▖▌ ▖▌▌ ▖  ▖▖     ▌▖ ▖   ▖ ▌▖ ▌▖  ▌
0    5*   5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5    5
";

        assert_eq!(good_output, std::str::from_utf8(&output).unwrap());
    }
}
