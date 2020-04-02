/* project use */
use cocktail;

/* criterion use */
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};

fn rev(c: &mut Criterion) {
    let mut group = c.benchmark_group("rev");
    for k in (1..32).step_by(2) {
        group.bench_with_input(BenchmarkId::new("pub", k), &k, |b, &k| {
            b.iter(|| cocktail::kmer::rev(black_box(1445814085), black_box(k)))
        });
        group.bench_with_input(BenchmarkId::new("loop", k), &k, |b, &k| {
            b.iter(|| cocktail::kmer::loop_rev(black_box(1445814085), black_box(k)))
        });
        group.bench_with_input(BenchmarkId::new("unrool", k), &k, |b, &k| {
            b.iter(|| cocktail::kmer::unrool_rev(black_box(1445814085), black_box(k)))
        });
    }
}

criterion_group!(benches, rev);

criterion_main!(benches);
