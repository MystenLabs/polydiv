use std::ops::Mul;

use fastcrypto::error::FastCryptoResult;
use fastcrypto::groups::bls12381::{G1Element, G2Element, Scalar};
use fastcrypto::groups::{GroupElement, MultiScalarMul, Pairing, Scalar as OtherScalar};
use rand::thread_rng;

use crate::fft::{BLS12381Domain, FFTDomain};
use crate::KZG;

/// Computes the matrix-vector multiplication for testing purposes -
// this is the function that is currently used for open_all
fn compute_matrix_vector_multiplication(
    coefficients: &[Scalar],
    tau_powers: &[G1Element],
) -> Vec<G1Element> {
    let d = coefficients.len();
    let mut result = vec![G1Element::zero(); d];
    let mut matrix = vec![vec![Scalar::zero(); d]; d];
    let mut coeff = coefficients.to_vec();
    coeff.reverse();

    for i in 0..d {
        matrix[i] = coeff.clone();
        coeff.pop();
        coeff.insert(0, Scalar::zero());
    }

    for i in 0..d {
        for j in 0..d {
            result[i] += tau_powers[j].mul(matrix[i][j]);
        }
    }

    result
}

/// Multiplies a Toeplitz matrix with a vector using FFT
/// This function optimizes the multiplication by using the Toeplitz structure and FFT
pub fn multiply_toeplitz_with_v(
    coefficients: &[Scalar],
    tau_powers: &[G1Element],
) -> Vec<G1Element> {
    let d = coefficients.len();
    let domain = BLS12381Domain::new(2 * d).unwrap();
    let mut result = vec![G1Element::zero(); d];

    // Construct a_2n vector based on Toeplitz matrix properties
    let mut a_2n = vec![Scalar::zero(); 2 * d];
    a_2n[0] = coefficients[d - 1];
    a_2n[d] = coefficients[d - 1];
    for i in 0..d - 1 {
        a_2n[d + 1 + i] = coefficients[i];
    }

    // Padding the tau_powers to match the length of a_2n
    let mut padded_x = vec![G1Element::zero(); 2 * d];
    for i in 0..d {
        padded_x[i] = tau_powers[i];
    }

    // Performing FFT on both vectors
    let v_vec = domain.fft(&a_2n);
    domain.fft_in_place_group(&mut padded_x);
    let y_vec = padded_x.clone();

    // Element-wise multiplication in the frequency domain
    let mut fft_result: Vec<G1Element> = v_vec
        .iter()
        .zip(y_vec.iter())
        .map(|(a, b)| b.mul(*a))
        .collect();

    // Inverse FFT to get the result in the time domain
    domain.ifft_in_place_group(&mut fft_result);

    for i in 0..d {
        result[i] = fft_result[i];
    }

    result
}

/// Struct for KZG commitments using Feist-Khovratovich technique
#[derive(Clone)]
pub struct KZGFK {
    domain: BLS12381Domain,
    tau_powers_g1: Vec<G1Element>,
    g2_tau: G2Element,
}

impl KZG for KZGFK {
    type G = G1Element;

    /// Creates a new KZGFK instance with a random tau
    fn new(n: usize) -> FastCryptoResult<Self> {
        let domain = BLS12381Domain::new(n)?;
        let n_dom = domain.size();

        let tau = fastcrypto::groups::bls12381::Scalar::rand(&mut thread_rng());

        let tau_powers_g1: Vec<G1Element> =
            itertools::iterate(G1Element::generator(), |g| g.mul(tau))
                .take(n_dom)
                .collect();

        let g2_tau = G2Element::generator().mul(tau);

        Ok(Self {
            domain,
            tau_powers_g1,
            g2_tau,
        })
    }

    /// Commits to a vector using the KZG commitment scheme
    fn commit(&self, v: &[Scalar]) -> G1Element {
        let poly = self.domain.ifft(v);
        G1Element::multi_scalar_mul(&poly, &self.tau_powers_g1).unwrap()
    }

    /// Opens a KZG commitment at a specific index
    fn open(&self, v: &[Scalar], index: usize) -> G1Element {
        let poly = self.domain.ifft(v);
        let mut quotient_coeffs: Vec<Scalar> = vec![Scalar::zero(); poly.len() - 1];
        quotient_coeffs[poly.len() - 2] = poly[poly.len() - 1];

        let omega_i = self.domain.element(index);
        for j in (0..poly.len() - 2).rev() {
            quotient_coeffs[j] = poly[j + 1] + quotient_coeffs[j + 1] * omega_i;
        }

        G1Element::multi_scalar_mul(
            &quotient_coeffs,
            &self.tau_powers_g1[..quotient_coeffs.len()],
        )
        .unwrap()
    }

    /// Opens a KZG commitment at multiple indices
    fn open_all(&self, v: &[Scalar]) -> Vec<G1Element> {
        let poly = self.domain.ifft(v);
        let degree = poly.len() - 1;

        let mut t = self.tau_powers_g1.clone();
        t.truncate(degree);
        t.reverse();

        let test_poly = &poly[poly.len() - degree..];

        let h_test = compute_matrix_vector_multiplication(test_poly, &t);

        let mut result = vec![G1Element::zero(); poly.len()];
        for i in 0..degree {
            result[i] = h_test[i];
        }

        self.domain.fft_in_place_group(&mut result);
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
        let rhs = self.g2_tau - G2Element::generator() * self.domain.element(index);

        lhs.pairing(&G2Element::generator()) == open_i.pairing(&rhs)
    }

    fn update(
        &self,
        commitment: &mut G1Element,
        _index: usize,
        _old_v_i: &Scalar,
        _new_v_i: &Scalar,
    ) -> G1Element {
        *commitment
    }

    fn update_open_i(
        &self,
        open: &G1Element,
        _index: usize,
        _old_v_i: &Scalar,
        _new_v_i: &Scalar,
    ) -> G1Element {
        *open
    }

    fn update_open_j(
        &self,
        open: &G1Element,
        _index: usize,
        _index_j: usize,
        _old_v_j: &Scalar,
        _new_v_j: &Scalar,
    ) -> G1Element {
        *open
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
        let n = 8;
        let kzg = KZGFK::new(n).unwrap();
        let v: Vec<Scalar> = (0..n).map(|_| OtherScalar::rand(&mut rng)).collect();
        let commitment = kzg.commit(&v);
        let index = rng.gen_range(0..n);
        let open_value = kzg.open(&v, index);
        let is_valid = kzg.verify(index, &v[index], &commitment, &open_value);
        assert!(is_valid, "Verification of the opening should succeed.");
    }

    #[test]
    fn test_kzg_commit_open_all() {
        let mut rng = rand::thread_rng();
        let n = 9;
        let kzg = KZGFK::new(n).unwrap();
        let v: Vec<Scalar> = (0..n).map(|_| OtherScalar::rand(&mut rng)).collect();
        let commitment = kzg.commit(&v);
        let mut open_values = kzg.open_all(&v);

        open_values.truncate(n);

        for (i, open_value) in open_values.iter().enumerate() {
            let is_valid = kzg.verify(i, &v[i], &commitment, open_value);
            assert!(
                is_valid,
                "Verification of the opening should succeed for index {}",
                i
            );
        }
    }

    #[test]
    fn test_check_toeplitz() {
        let v_scalar = vec![
            Scalar::from(1),
            Scalar::from(2),
            Scalar::from(3),
            Scalar::from(4),
            Scalar::from(4),
            Scalar::from(4),
            Scalar::from(4),
            Scalar::from(4),
        ];

        let tau = Scalar::rand(&mut thread_rng());

        let tau_powers_g1: Vec<G1Element> =
            itertools::iterate(G1Element::generator(), |g| g.mul(tau))
                .take(8)
                .collect();

        let h_alin = multiply_toeplitz_with_v(&v_scalar, &tau_powers_g1);
        let h_mult = compute_matrix_vector_multiplication(&v_scalar, &tau_powers_g1);

        assert!(h_alin == h_mult, "Toeplitz multiplication mismatch.");
    }
}
