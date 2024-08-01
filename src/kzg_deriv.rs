use std::ops::Mul;
use std::sync::{Arc, Mutex};

use fastcrypto::error::FastCryptoResult;
use fastcrypto::groups::bls12381::{G1Element, G2Element, Scalar};
use fastcrypto::groups::{GroupElement, MultiScalarMul, Pairing, Scalar as OtherScalar};
use itertools::iterate;
use rand::thread_rng;
use rayon::prelude::*;

use crate::fft::{BLS12381Domain, FFTDomain};
use crate::KZG;

/// Adds three vectors element-wise
fn add_vectors(v1: Vec<G1Element>, v2: Vec<G1Element>, v3: Vec<G1Element>) -> Vec<G1Element> {
    v1.iter()
        .zip(v2.iter())
        .zip(v3.iter())
        .map(|((a, b), c)| a + b + c)
        .collect()
}

/// Performs sparse matrix-vector multiplication
fn sparse_c_matrix_vector_multiply(vector: &mut Vec<G1Element>) {
    let n = vector.len();
    let two_inverse = Scalar::from(2u128).inverse().unwrap();
    let coefficient = Scalar::from((n + 1) as u128) * two_inverse;

    vector.rotate_right(1);
    vector.iter_mut().enumerate().for_each(|(i, v)| {
        if i == 0 {
            *v = v.mul(-Scalar::from((n - 1) as u128) * two_inverse)
        } else {
            *v = v.mul(coefficient - Scalar::from(i as u128))
        }
    });
}

fn sparse_d_matrix_vector_multiply(vector: &mut Vec<G1Element>) {
    let n = vector.len();
    let two_inverse = Scalar::from(2u128).inverse().unwrap();
    let coefficient = -Scalar::from((n + 1) as u128) * two_inverse;

    vector.rotate_right(1);
    vector.iter_mut().enumerate().for_each(|(i, v)| {
        if i == 0 {
            *v = v.mul((Scalar::from((n - 1) as u128)) * two_inverse)
        } else {
            *v = v.mul(coefficient + Scalar::from(i as u128))
        }
    });
}

/// Multiplies a diagonal matrix by a vector
fn multiply_d_matrix_by_vector(vector: &[Scalar]) -> Vec<Scalar> {
    vector[1..]
        .iter()
        .enumerate()
        .map(|(i, v)| Scalar::from((i + 1) as u128) * v)
        .collect()
}

/// Struct for KZG commitment scheme using derived elements
#[derive(Clone)]
pub struct KZGDeriv {
    domain: BLS12381Domain,
    n: usize,
    omega: Scalar,
    g2_tau: G2Element,
    w_vec: Vec<G1Element>,
    u_vec: Vec<G1Element>,
    omega_powers: Vec<Scalar>,
}

impl KZGDeriv {
    /// Computes the omega^index element
    fn element(&self, index: usize) -> Scalar {
        self.domain.element(index)
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
        let tau_powers_g: Vec<Scalar> = iterate(Scalar::generator(), |g| g * tau).take(n).collect();

        let g = G1Element::generator();
        let w_vec: Vec<G1Element> = domain
            .ifft(&tau_powers_g)
            .iter()
            .map(|s| g.mul(s))
            .collect();

        let mut omega_i = Scalar::generator();
        let omega = domain.element(1);

        // compute uvec using w_vec - original implementation
        let u_vec: Vec<G1Element> = (0..n)
            .map(|i| {
                let l_i_minus_1 = w_vec[i] - G1Element::generator();
                let denom = tau - omega_i;
                let u_i = (l_i_minus_1 / denom).unwrap();

                if i < n - 1 {
                    omega_i *= omega;
                }
                u_i
            })
            .collect();

        let omega_powers = iterate(Scalar::generator(), |g| g * omega)
            .take(n)
            .collect();

        Ok(Self {
            domain,
            n,
            omega,
            g2_tau,
            w_vec,
            u_vec,
            omega_powers,
        })
    }

    /// Commits to a vector using the KZG commitment scheme
    fn commit(&self, v: &[Scalar]) -> G1Element {
        G1Element::multi_scalar_mul(&v, &self.w_vec).unwrap()
    }

    /// Opens a KZG commitment at a specific index
    fn open(&self, v: &[Scalar], index: usize) -> G1Element {
        let (mut scalars, v_prime_terms): (Vec<Scalar>, Vec<Scalar>) = (0..v.len())
            .into_par_iter()
            .map(|j| {
                if j != index {
                    let diff_inverse = (self.omega_powers[index] - self.omega_powers[j])
                        .inverse()
                        .unwrap();
                    (
                        (v[index] - v[j]) * diff_inverse,
                        v[j] * self.omega_powers
                            [(self.domain.size() + j - index) % self.domain.size()]
                            * diff_inverse,
                    )
                } else {
                    (
                        Scalar::zero(),
                        v[j] * (Scalar::from((v.len() - 1) as u128)
                            / (Scalar::from(2u128) * self.omega_powers[index]))
                            .unwrap(),
                    )
                }
            })
            .collect();

        scalars[index] = v_prime_terms.into_iter().reduce(|a, b| a + b).unwrap();

        G1Element::multi_scalar_mul(&scalars, &self.w_vec).unwrap()
    }

    /// Opens a KZG commitment at multiple indices
    fn open_all(&self, v: &[Scalar], indices: &[usize]) -> Vec<G1Element> {
        // Compute tau * Dhatv
        let idftv = self.domain.ifft(&v);
        let d_msm_idftv: Vec<Scalar> = multiply_d_matrix_by_vector(&idftv);
        let dhatv = self.domain.fft(&d_msm_idftv);
        let result1: Vec<G1Element> = self
            .w_vec
            .iter()
            .zip(dhatv.iter())
            .map(|(a, b)| a.mul(*b))
            .collect();

        // Compute ColEDiv.tau*v
        let mut powtau = self.w_vec.clone();
        self.domain.fft_in_place_group(&mut powtau);
        let mut col_hat_dft_tau = powtau;
        sparse_c_matrix_vector_multiply(&mut col_hat_dft_tau);
        self.domain.ifft_in_place_group(&mut col_hat_dft_tau);
        let result2: Vec<G1Element> = col_hat_dft_tau
            .iter()
            .zip(v.iter())
            .map(|(a, b)| a.mul(*b))
            .collect();

        // Compute diadiv.powtau*v
        let mut mult: Vec<G1Element> = self
            .w_vec
            .iter()
            .zip(v.iter())
            .map(|(a, b)| a.mul(*b))
            .collect();

        self.domain.fft_in_place_group(&mut mult);
        let mut diadiv_idft_tau_v: Vec<G1Element> = mult;
        sparse_d_matrix_vector_multiply(&mut diadiv_idft_tau_v);
        self.domain.ifft_in_place_group(&mut diadiv_idft_tau_v);

        let result3 = diadiv_idft_tau_v;

        let result = add_vectors(result1, result2, result3);

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
        open: &G1Element,
        index: usize,
        old_v_i: &Scalar,
        new_v_i: &Scalar,
    ) -> G1Element {
        *open + self.u_vec[index].mul(new_v_i - old_v_i)
    }

    fn update_open_j(
        &self,
        open: &G1Element,
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
        assert!(
            is_valid,
            "Verification of the opening after updating should succeed."
        );
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
        assert!(
            is_valid,
            "Verification of the opening after updating j's value should succeed."
        );
    }

    #[test]
    fn test_kzg_commit_open_all() {
        let mut rng = rand::thread_rng();
        let n = 8;
        let kzg = KZGDeriv::new(n).unwrap();
        let v: Vec<Scalar> = (0..n).map(|_| OtherScalar::rand(&mut rng)).collect();
        let commitment = kzg.commit(&v);
        let indices: Vec<usize> = (0..n).collect();
        let open_values = kzg.open_all(&v, &indices);

        for (i, open_value) in open_values.iter().enumerate() {
            let is_valid = kzg.verify(indices[i], &v[indices[i]], &commitment, open_value);
            assert!(
                is_valid,
                "Verification of the opening should succeed for index {}",
                i
            );
        }
    }
}
