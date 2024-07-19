use std::ops::Mul;
use fastcrypto::error::{FastCryptoResult, FastCryptoError};
use fastcrypto::groups::bls12381::{G1Element, G2Element, Scalar};
use fastcrypto::groups::{GroupElement, MultiScalarMul, Pairing, Scalar as OtherScalar};
use itertools::iterate;
use rand::thread_rng;
use rayon::prelude::*;
use std::sync::{Mutex, Arc};

use crate::fft::{BLS12381Domain, FFTDomain};
use crate::KZG;

/// Adds three vectors element-wise
fn add_vectors(v1: Vec<G1Element>, v2: Vec<G1Element>, v3: Vec<G1Element>) -> Vec<G1Element> {
    v1.into_iter().zip(v2.into_iter()).zip(v3.into_iter())
        .map(|((a, b), c)| a + b + c)
        .collect()
}

/// Performs sparse matrix-vector multiplication
fn sparse_c_matrix_vector_multiply(vector: &Vec<Scalar>, n: usize) -> Vec<Scalar> {
    let mut result = vec![Scalar::zero(); n];

    // Handle the special element in the first row, last column
    result[0] += -((Scalar::from((n - 1) as u128)) / (Scalar::from(2 as u128))).unwrap() * vector[n - 1];

    // Handle the diagonal elements below the main diagonal
    for i in 1..n {
        result[i] += (-Scalar::from(i as u128) + (Scalar::from((n + 1) as u128) / Scalar::from(2u128)).unwrap()) * vector[i - 1];
    }

    result
}

/// Multiplies two matrices
fn matrix_multiply(a: &Vec<Vec<Scalar>>, b: &Vec<Vec<Scalar>>) -> Vec<Vec<Scalar>> {
    let rows_a = a.len();
    let cols_a = a[0].len();
    let cols_b = b[0].len();

    let mut result = vec![vec![Scalar::zero(); cols_b]; rows_a];

    for i in 0..rows_a {
        for j in 0..cols_b {
            for k in 0..cols_a {
                result[i][j] += a[i][k] * b[k][j];
            }
        }
    }

    result
}

/// Multiplies a matrix of scalars by a vector of scalars
fn multiply_matrix_scalar_vector(matrix: &[Vec<Scalar>], v: &[Scalar]) -> Vec<Scalar> {
    let n = matrix.len();
    let mut result = vec![Scalar::zero(); n];

    for i in 0..n {
        for j in 0..v.len() {
            if matrix[i][j] != Scalar::zero() {
                result[i] = result[i] + matrix[i][j].mul(v[j]);
            }
        }
    }

    result
}

/// Multiplies a diagonal matrix by a vector
fn multiply_d_matrix_by_vector(vector: &[Scalar]) -> Vec<Scalar> {
    let n = vector.len();
    let mut result = vec![Scalar::zero(); n];

    for i in 0..n - 1 {
        result[i] = Scalar::from((i + 1) as u128) * vector[i + 1];
    }

    result
}

/// Multiplies a matrix of scalars by a vector of group elements
fn multiply_matrix_vector(matrix: &[Vec<Scalar>], v: &[G1Element]) -> Vec<G1Element> {
    let n = v.len();
    let mut result = vec![G1Element::zero(); n];

    for i in 0..n {
        let scalars = &matrix[i];
        result[i] = G1Element::multi_scalar_mul(&scalars, &v).unwrap();
    }

    result
}

/// Struct for KZG commitment scheme using derived elements
#[derive(Clone)]
pub struct KZGDeriv {
    domain: BLS12381Domain,
    n: usize,
    omega: Scalar,
    g2_tau: G2Element,
    w_vec: Vec<G1Element>,
    check_w_vec: Vec<G1Element>,
    u_vec: Vec<G1Element>,
    tau_powers_g1: Vec<G1Element>,
}

impl KZGDeriv {
    /// Computes the omega^index element
    fn element(&self, index: usize) -> Scalar {
        if index > self.n {
            return self.element(index & self.n);
        } else if index == 0 {
            return Scalar::generator();
        } else if index == 1 {
            return self.omega;
        }

        let half = self.element(index / 2);
        if index % 2 == 0 {
            return half * half;
        }
        half * half * self.omega
    }
}

impl KZG for KZGDeriv {
    type G = G1Element;

    /// Creates a new KZGDeriv instance with a random tau
    fn new(n: usize) -> FastCryptoResult<Self> {
        let domain = BLS12381Domain::new(n)?;

        let tau = Scalar::rand(&mut thread_rng());
        let g2_tau = G2Element::generator() * tau;

        // Compute tau^i for i = 0 to n-1
        let tau_powers_g: Vec<Scalar> =
            iterate(Scalar::generator(), |g| g * tau).take(n).collect();

        let g = G1Element::generator();
        let w_vec: Vec<G1Element> = domain
            .ifft(&tau_powers_g)
            .iter()
            .map(|s| g.mul(s))
            .collect();

        let mut omega_i = Scalar::generator();
        let mut u_vec = (0..n)
            .map(|i| {
                let l_i_minus_1 = w_vec[i] - G1Element::generator();
                let denom = tau - omega_i;
                let u_i = (l_i_minus_1 / denom).unwrap();

                if i < n - 1 {
                    omega_i *= domain.element(1);
                }
                u_i
            })
            .collect();

        let omega = domain.element(1);

        let tau_powers_g1: Vec<G1Element> = tau_powers_g.iter().map(|x| g.mul(x)).collect();

        let mut check_w_vec = tau_powers_g1.clone();
        domain.ifft_in_place_g1(&mut check_w_vec);

        Ok(Self {
            domain,
            n,
            omega,
            g2_tau,
            w_vec,
            check_w_vec,
            u_vec,
            tau_powers_g1
        })
    }

    /// Commits to a vector using the KZG commitment scheme
    fn commit(&self, v: &[Scalar]) -> G1Element {
        G1Element::multi_scalar_mul(&v, &self.w_vec).unwrap()
    }

    /// Opens a KZG commitment at a specific index
    fn open(&self, v: &[Scalar], index: usize) -> G1Element {
        let scalars = Arc::new(Mutex::new(vec![Scalar::zero(); v.len()]));
        let v_prime = Arc::new(Mutex::new(Scalar::zero()));

        // Initialize omega_i and omega_j_minus_i
        let omega_i = self.element(index);
        let omega_base = Scalar::generator();
        let omega_j_minus_i_base = self.element(self.n - index);

        (0..v.len()).into_par_iter().for_each(|j| {
            let mut local_omega_j = omega_base;
            let mut local_omega_j_minus_i = omega_j_minus_i_base;

            // Compute local omegas
            for _ in 0..j {
                local_omega_j *= self.omega;
                local_omega_j_minus_i *= self.omega;
            }

            if j != index {
                let diff_inverse = (omega_i - local_omega_j).inverse().unwrap();
                let mut v_prime_guard = v_prime.lock().unwrap();
                *v_prime_guard += v[j] * local_omega_j_minus_i * diff_inverse;

                let mut scalars_guard = scalars.lock().unwrap();
                scalars_guard[j] = (v[index] - v[j]) * diff_inverse;
            } else {
                let mut v_prime_guard = v_prime.lock().unwrap();
                *v_prime_guard += v[j]
                    * (Scalar::from((v.len() - 1) as u128) / (Scalar::from(2u128) * omega_i))
                        .unwrap();
            }
        });

        let v_prime = v_prime.lock().unwrap().clone();
        let mut scalars_guard = scalars.lock().unwrap();
        scalars_guard[index] = v_prime;

        G1Element::multi_scalar_mul(&scalars_guard, &self.w_vec[..scalars_guard.len()]).unwrap()
    }

    /// Opens a KZG commitment at multiple indices
    fn open_all(&self, v: &[Scalar], indices: Vec<usize>) -> Vec<G1Element> {
        let n = v.len();

        // Compute ColEDivHat matrix
        let mut c_matrix = vec![vec![Scalar::zero(); n]; n];
        for i in 0..n {
            for j in 0..n {
                if i == 0 && j == n - 1 {
                    c_matrix[i][j] = -((Scalar::from((n - 1) as u128)) / (Scalar::from(2 as u128))).unwrap();
                } else if j + 1 == i && i >= 1 && i < n {
                    c_matrix[i][j] = -Scalar::from(i as u128) + (Scalar::from((n + 1) as u128) / Scalar::from(2u128)).unwrap();
                } else {
                    c_matrix[i][j] = Scalar::zero();
                }
            }
        }

        // Compute D matrix
        let mut d_matrix = vec![vec![Scalar::zero(); n]; n];
        for i in 0..n {
            for j in 0..n {
                if i == 0 && j == n - 1 {
                    d_matrix[i][j] = ((Scalar::from((n - 1) as u128)) / (Scalar::from(2 as u128))).unwrap();
                } else if j + 1 == i && i >= 1 && i < n {
                    d_matrix[i][j] = Scalar::from(i as u128) - (Scalar::from((n + 1) as u128) / Scalar::from(2u128)).unwrap();
                } else {
                    d_matrix[i][j] = Scalar::zero();
                }
            }
        }

        // Compute tau * Dhatv
        let idftv = self.domain.ifft(&v);
        let d_msm_idftv: Vec<Scalar> = multiply_d_matrix_by_vector(&idftv);
        let dhatv = self.domain.fft(&d_msm_idftv);
        let result1: Vec<G1Element> = self.w_vec.iter()
            .zip(dhatv.iter())
            .map(|(a, b)| a.mul(*b))
            .collect();

        // Compute ColEDiv.tau*v
        let mut powtau = self.w_vec.clone();
        self.domain.fft_in_place_g1(&mut powtau);
        let mut col_hat_dft_tau: Vec<G1Element> = multiply_matrix_vector(&c_matrix, &powtau);
        self.domain.ifft_in_place_g1(&mut col_hat_dft_tau);
        let result2: Vec<G1Element> = col_hat_dft_tau.iter()
            .zip(v.iter())
            .map(|(a, b)| a.mul(*b))
            .collect();

        // Compute diadiv.powtau*v
        let mut mult: Vec<G1Element> = self.w_vec.iter()
            .zip(v.iter())
            .map(|(a, b)| a.mul(*b))
            .collect();

        self.domain.fft_in_place_g1(&mut mult);
        let mut diadiv_idft_tau_v: Vec<G1Element> = multiply_matrix_vector(&d_matrix, &mult);
        self.domain.ifft_in_place_g1(&mut diadiv_idft_tau_v);

        let result3 = diadiv_idft_tau_v.clone();

        let result = add_vectors(result1.clone(), result2.clone(), result3.clone());

        result
    }

    /// Verifies a KZG opening
    fn verify(
        &self,
        index: usize,
        v_i: &Scalar,
        commitment: &G1Element,
        open_i: &G1Element,
    ) -> bool {
        let lhs = *commitment - G1Element::generator() * v_i;
        let rhs = self.g2_tau - G2Element::generator() * self.element(index);

        lhs.pairing(&G2Element::generator()) == open_i.pairing(&rhs)
    }

    fn update(
        &self,
        commitment: &mut G1Element,
        index: usize,
        old_v_i: &Scalar,
        new_v_i: &Scalar,
    ) -> G1Element {
        *commitment + self.w_vec[index].mul(new_v_i - old_v_i)
    }

    fn update_open_i(
        &self,
        open: &mut G1Element,
        index: usize,
        old_v_i: &Scalar,
        new_v_i: &Scalar,
    ) -> G1Element {
        *open + self.u_vec[index].mul(new_v_i - old_v_i)
    }

    fn update_open_j(
        &self,
        open: &mut G1Element,
        index: usize,
        index_j: usize,
        old_v_j: &Scalar,
        new_v_j: &Scalar,
    ) -> G1Element {
        assert_ne!(index, index_j, "index and index_j should be different.");

        let omega_i = self.element(index);
        let omega_j = self.element(index_j);
        let omega_i_inverse = self.element(self.n - index);

        let to_mul_1 = (old_v_j - new_v_j) * (omega_i - omega_j).inverse().unwrap();
        let to_mul_2 = -to_mul_1 * omega_j * omega_i_inverse;

        *open
            + G1Element::multi_scalar_mul(
                &[to_mul_1, to_mul_2],
                &[self.w_vec[index_j], self.w_vec[index]],
            )
            .unwrap()
    }
}

#[cfg(test)]
mod tests {
    use fastcrypto::groups::bls12381::Scalar;
    use rand::Rng;

    use super::*;

    #[test]
    fn test_kzg_commit_open_verify() {
        let mut rng = rand::thread_rng();
        let n = 4;
        let kzg = KZGDeriv::new(n).unwrap();
        let v: Vec<Scalar> = (0..n).map(|_| OtherScalar::rand(&mut rng)).collect();
        let commitment = kzg.commit(&v);
        let index = rng.gen_range(0..n);
        let open_value = kzg.open(&v, index);
        let is_valid = kzg.verify(index, &v[index], &commitment, &open_value);
        assert!(is_valid, "Verification of the opening should succeed.");
    }

    #[test]
    fn test_kzg_commit_open_update_i_verify() {
        let mut rng = rand::thread_rng();
        let n = 8;
        let kzg = KZGDeriv::new(n).unwrap();
        let v: Vec<Scalar> = (0..n).map(|_| OtherScalar::rand(&mut rng)).collect();
        let mut commitment = kzg.commit(&v);
        let index = rng.gen_range(0..n);
        let mut open_value = kzg.open(&v, index);
        let new_v_index = Scalar::rand(&mut rng);
        let new_commitment = kzg.update(&mut commitment, index, &v[index], &new_v_index);
        let new_opening = kzg.update_open_i(&mut open_value, index, &v[index], &new_v_index);
        let is_valid = kzg.verify(index, &new_v_index, &new_commitment, &new_opening);
        assert!(is_valid, "Verification of the opening after updating should succeed.");
    }

    #[test]
    fn test_kzg_commit_open_update_j_verify() {
        let mut rng = rand::thread_rng();
        let n = 8;
        let kzg = KZGDeriv::new(n).unwrap();
        let v: Vec<Scalar> = (0..n).map(|_| OtherScalar::rand(&mut rng)).collect();
        let mut commitment = kzg.commit(&v);
        let index = rng.gen_range(0..n);
        let mut open_value = kzg.open(&v, index);

        let mut index_j;
        loop {
            index_j = rng.gen_range(0..n);
            if index_j != index {
                break;
            }
        }

        let new_v_index_j = Scalar::rand(&mut rng);
        let new_commitment = kzg.update(&mut commitment, index_j, &v[index_j], &new_v_index_j);
        let new_opening =
            kzg.update_open_j(&mut open_value, index, index_j, &v[index_j], &new_v_index_j);
        let is_valid = kzg.verify(index, &v[index], &new_commitment, &new_opening);
        assert!(is_valid, "Verification of the opening after updating j's value should succeed.");
    }

    #[test]
    fn test_kzg_commit_open_all() {
        let mut rng = rand::thread_rng();
        let n = 8;
        let kzg = KZGDeriv::new(n).unwrap();
        let v: Vec<Scalar> = (0..n).map(|_| OtherScalar::rand(&mut rng)).collect();
        let commitment = kzg.commit(&v);
        let indices: Vec<usize> = (0..n).collect();
        let open_values = kzg.open_all(&v, indices.clone());

        for (i, open_value) in open_values.iter().enumerate() {
            let is_valid = kzg.verify(indices[i], &v[indices[i]], &commitment, open_value);
            assert!(is_valid, "Verification of the opening should succeed for index {}", i);
        }
    }
}
