/* crates use */

use rand::prelude::IndexedRandom as _;

/* project use */

/* criterion use */
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};

mod kmer2seq;

fn kmer2seq(c: &mut Criterion) {
    assert_eq!(
        kmer2seq::static_buffer(7271, 5),
        cocktail::kmer::kmer2seq(7271, 5)
    );
    assert_eq!(
        kmer2seq::static_buffer(7271, 5),
        kmer2seq::local_buffer(7271, 5)
    );
    assert_eq!(
        kmer2seq::static_buffer(7271, 5),
        kmer2seq::dyn_local_buffer(7271, 5)
    );

    let mut g = c.benchmark_group("kmer2seq");

    for k in (1..32).step_by(2) {
        g.bench_with_input(BenchmarkId::new("actual implementation", k), &k, |b, &k| {
            b.iter(|| std::hint::black_box(cocktail::kmer::kmer2seq(std::hint::black_box(7271), k)))
        });

        g.bench_with_input(BenchmarkId::new("static buffer", k), &k, |b, &k| {
            b.iter(|| std::hint::black_box(kmer2seq::static_buffer(std::hint::black_box(7271), k)))
        });

        g.bench_with_input(BenchmarkId::new("local buffer", k), &k, |b, &k| {
            b.iter(|| std::hint::black_box(kmer2seq::local_buffer(std::hint::black_box(7271), k)))
        });

        g.bench_with_input(BenchmarkId::new("dynamic local buffer", k), &k, |b, &k| {
            b.iter(|| {
                std::hint::black_box(kmer2seq::dyn_local_buffer(std::hint::black_box(7271), k))
            })
        });
    }
}

mod iter_cano;

fn tokenize_canonical(c: &mut Criterion) {
    for k in (11..19).step_by(2) {
        let mut g = c.benchmark_group(format!("canonical kmer iteration k={k}"));

        let mut rng = rand::rng();
        let vals = [b'A', b'C', b'G', b'T'];

        for i in 7..14 {
            let len = 1 << i;

            let seq = (0..len)
                .map(|_| *vals.choose(&mut rng).unwrap())
                .collect::<Vec<u8>>();

            g.bench_with_input(BenchmarkId::new("only forward", len), &seq, |b, seq| {
                std::hint::black_box(b.iter(|| {
                    cocktail::tokenizer::basic::Tokenizer::new(
                        std::hint::black_box(seq),
                        std::hint::black_box(k),
                    )
                    .map(|x| cocktail::kmer::canonical(x, k))
                    .collect::<Vec<u64>>()
                }))
            });

            g.bench_with_input(BenchmarkId::new("forward reverse", len), &seq, |b, seq| {
                b.iter(|| {
                    std::hint::black_box(
                        cocktail::tokenizer::kmer::Canonical::new(
                            std::hint::black_box(seq),
                            std::hint::black_box(k),
                        )
                        .collect::<Vec<u64>>(),
                    )
                })
            });

            g.bench_with_input(
                BenchmarkId::new("forward reverse lexi", len),
                &seq,
                |b, seq| {
                    b.iter(|| {
                        std::hint::black_box(
                            iter_cano::Lexi::new(
                                std::hint::black_box(seq),
                                std::hint::black_box(k),
                            )
                            .collect::<Vec<u64>>(),
                        )
                    })
                },
            );
        }
    }
}

fn setup(c: &mut Criterion) {
    tokenize_canonical(c);
    kmer2seq(c);
}

criterion_group!(benches, setup);

criterion_main!(benches);
