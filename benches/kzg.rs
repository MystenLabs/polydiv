use criterion::measurement::Measurement;
use criterion::{criterion_group, criterion_main, BenchmarkGroup, Criterion};
use fastcrypto::groups::{GroupElement, Scalar};
use fastcrypto_kzg::kzg_deriv::KZGDeriv;
use fastcrypto_kzg::kzg_fk::KZGFK;
use fastcrypto_kzg::kzg_original::KZGOriginal;
use fastcrypto_kzg::kzg_tabdfk::KZGTabDFK;
use fastcrypto_kzg::KZG;
use rand::{thread_rng, Rng};

// Adjust the imports based on your actual project structure

fn kzg_single<K: KZG, M: Measurement>(name: &str, c: &mut BenchmarkGroup<M>) {
    let input_sizes = [128]; //, 128, 1024, 8192, 65536];

    for size in input_sizes {
        c.bench_function(format!("{}/new/{}", name, size), |b| {
            b.iter(|| K::new(size).unwrap());
        });
        let kzg = K::new(size).unwrap();

        let mut rng = thread_rng();
        let data = vec![<<K as KZG>::G as GroupElement>::ScalarType::rand(&mut rng); size];

        c.bench_function(format!("{}/commit/{}", name, size), |b| {
            b.iter(|| kzg.commit(&data));
        });

        let mut commitment = kzg.commit(&data);

        // Pick a random index to open
        let index = rng.gen_range(0..size);

        // Create an opening
        c.bench_function(format!("{}/open/{}", name, size), |b| {
            b.iter(|| kzg.open(&data, index));
        });
        let mut open_value = kzg.open(&data, index);

        // Pick a new index to update¨
        let mut index_j;
        loop {
            index_j = rng.gen_range(0..size);
            if index_j != index {
                break;
            }
        }

        // Set a new value for v_i
        let new_v_index_j = <<K as KZG>::G as GroupElement>::ScalarType::rand(&mut rng);

        // Update the commitment
        c.bench_function(format!("{}/update/{}", name, size), |b| {
            b.iter(|| kzg.update(&mut commitment, index_j, &data[index_j], &new_v_index_j));
        });
        let new_commitment = kzg.update(&mut commitment, index_j, &data[index_j], &new_v_index_j);

        // Update the opening
        c.bench_function(format!("{}/update_open_j/{}", name, size), |b| {
            b.iter(|| {
                kzg.update_open_j(
                    &mut open_value,
                    index,
                    index_j,
                    &data[index_j],
                    &new_v_index_j,
                )
            });
        });
        let new_opening = kzg.update_open_j(
            &mut open_value,
            index,
            index_j,
            &data[index_j],
            &new_v_index_j,
        );

        // Verify the opening
        c.bench_function(format!("{}/verify/{}", name, size), |b| {
            b.iter(|| kzg.verify(index, &data[index], &new_commitment, &new_opening));
        });
    }
}

fn kzg(c: &mut Criterion) {
    let mut group = c.benchmark_group("KZG".to_string());

    kzg_single::<KZGDeriv, _>("KZGDeriv", &mut group);
    kzg_single::<KZGOriginal, _>("KZGOriginal", &mut group);
    kzg_single::<KZGTabDFK, _>("KZGTabDFK", &mut group);
    kzg_single::<KZGFK, _>("KZGFK", &mut group);
}

criterion_group! {
    name = kzg_benchmarks;
    config = Criterion::default().sample_size(100);
    targets = kzg
}

criterion_main!(kzg_benchmarks);
