extern crate test;

use test::Bencher;

use my_project::{KZG, KZGFK, KZGOriginal, KZGTabDFK, Scalar};

// Adjust the imports based on your actual project structure

#[bench]
fn bench_commit_original(b: &mut Bencher) {
    let kzg = KZGOriginal::new(4);
    let data = vec![Scalar::generator(); 4];

    b.iter(|| {
        kzg.commit(&data)
    });
}

#[bench]
fn bench_commit_fk(b: &mut Bencher) {
    let kzg = KZGFK::new(4);
    let data = vec![Scalar::generator(); 4];

    b.iter(|| {
        kzg.commit(&data)
    });
}

#[bench]
fn bench_commit_tabdfk(b: &mut Bencher) {
    let kzg = KZGTabDFK::new(4);
    let data = vec![Scalar::generator(); 4];

    b.iter(|| {
        kzg.commit(&data)
    });
}

// Add other benchmarks for `open`, `verify`, and `update` similarly for all implementations
