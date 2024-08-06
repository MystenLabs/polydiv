use criterion::measurement::Measurement;
use criterion::{criterion_group, criterion_main, BenchmarkGroup, Criterion};
use fastcrypto::groups::bls12381::Scalar as BLSScalar;
use fastcrypto::groups::{GroupElement, Scalar};
use fastcrypto_kzg::kzg_deriv::KZGDeriv;
use fastcrypto_kzg::kzg_fk::KZGFK;
use fastcrypto_kzg::kzg_original::KZGOriginal;
use fastcrypto_kzg::kzg_tabdfk::KZGTabDFK;
use fastcrypto_kzg::KZG;
use rand::{thread_rng, Rng};

// Adjust the imports based on your actual project structure

fn kzg_single<K: KZG, M: Measurement>(name: &str, c: &mut BenchmarkGroup<M>) {
    let input_sizes = [16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192];

    for &size in &input_sizes {
        c.bench_function(format!("{}/new/{}", name, size), |b| {
            b.iter(|| K::new(size).unwrap());
        });
        let kzg = K::new(size).unwrap();

        let mut rng = thread_rng();

        // Generate data for commit
        let commit_data: Vec<<K::G as GroupElement>::ScalarType> = (0..size)
            .map(|_| <<K::G as GroupElement>::ScalarType as Scalar>::rand(&mut rng))
            .collect();

        // Generate data for open_all
        let open_all_data: Vec<BLSScalar> = (0..size).map(|_| BLSScalar::rand(&mut rng)).collect();

        c.bench_function(format!("{}/commit/{}", name, size), |b| {
            b.iter(|| kzg.commit(&commit_data));
        });
        let mut commitment = kzg.commit(&commit_data);

        // Pick a random index to open
        let index = rng.gen_range(0..size);

        // Create an opening
        c.bench_function(format!("{}/open/{}", name, size), |b| {
            b.iter(|| kzg.open(&commit_data, index));
        });
        let mut open_value = kzg.open(&commit_data, index);

        // create all openings
        c.bench_function(format!("{}/open_all/{}", name, size), |b| {
            b.iter(|| kzg.open_all(&open_all_data));
        });
        let mut open_values = kzg.open_all(&open_all_data);

        // Pick a new index to update
        let mut index_j;
        loop {
            index_j = rng.gen_range(0..size);
            if index_j != index {
                break;
            }
        }

        // Set a new value for v_i
        let new_v_index_j = <<K::G as GroupElement>::ScalarType as Scalar>::rand(&mut rng);

        // Update the commitment
        c.bench_function(format!("{}/update/{}", name, size), |b| {
            b.iter(|| {
                kzg.update(
                    &mut commitment,
                    index_j,
                    &commit_data[index_j],
                    &new_v_index_j,
                )
            });
        });
        let new_commitment = kzg.update(
            &mut commitment,
            index_j,
            &commit_data[index_j],
            &new_v_index_j,
        );

        // Update the opening
        c.bench_function(format!("{}/update_open_j/{}", name, size), |b| {
            b.iter(|| {
                kzg.update_open_j(
                    &mut open_value,
                    index,
                    index_j,
                    &commit_data[index_j],
                    &new_v_index_j,
                )
            });
        });
        let new_opening = kzg.update_open_j(
            &mut open_value,
            index,
            index_j,
            &commit_data[index_j],
            &new_v_index_j,
        );

        // Verify the opening
        c.bench_function(format!("{}/verify/{}", name, size), |b| {
            b.iter(|| kzg.verify(index, &commit_data[index], &new_commitment, &new_opening));
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
    config = Criterion::default().sample_size(10);
    targets = kzg
}

criterion_main!(kzg_benchmarks);
