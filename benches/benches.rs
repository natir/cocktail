/* project use */
use cocktail;

/* criterion use */
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};

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
            b.iter(|| black_box(cocktail::kmer::kmer2seq(black_box(7271), k)))
        });

        g.bench_with_input(BenchmarkId::new("static buffer", k), &k, |b, &k| {
            b.iter(|| black_box(kmer2seq::static_buffer(black_box(7271), k)))
        });

        g.bench_with_input(BenchmarkId::new("local buffer", k), &k, |b, &k| {
            b.iter(|| black_box(kmer2seq::local_buffer(black_box(7271), k)))
        });

        g.bench_with_input(BenchmarkId::new("dynamic local buffer", k), &k, |b, &k| {
            b.iter(|| black_box(kmer2seq::dyn_local_buffer(black_box(7271), k)))
        });
    }
}

use rand::seq::SliceRandom;
use rand::Rng;

fn tokenize_canonical(c: &mut Criterion) {
    let mut g = c.benchmark_group("canonical kmer iteration");
    let mut rng = rand::thread_rng();
    let vals = [b'A', b'C', b'G', b'T'];

    for i in 3..20 {
        let len = 1 << i;

        let seq = (0..len)
            .map(|_| *vals.choose(&mut rng).unwrap())
            .collect::<Vec<u8>>();

        g.bench_with_input(BenchmarkId::new("after", len), &seq, |b, seq| {
            b.iter(|| {
                cocktail::tokenizer::Tokenizer::new(black_box(&seq), black_box(5))
                    .map(|x| cocktail::kmer::cannonical(x, 5))
                    .collect::<Vec<u64>>()
            })
        });

        g.bench_with_input(BenchmarkId::new("durring", len), &seq, |b, seq| {
            b.iter(|| {
                cocktail::tokenizer::Canonical::new(black_box(&seq), black_box(5))
                    .collect::<Vec<u64>>()
            })
        });
    }
}

fn setup(c: &mut Criterion) {
    tokenize_canonical(c);
    kmer2seq(c);
}

criterion_group!(benches, setup);

criterion_main!(benches);
