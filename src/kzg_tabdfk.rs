use std::ops::Mul;

use fastcrypto::error::FastCryptoResult;
use fastcrypto::groups::bls12381::{G1Element, G2Element, Scalar};
use fastcrypto::groups::{GroupElement, MultiScalarMul, Pairing, Scalar as OtherScalar};
use itertools::iterate;
use rand::thread_rng;

use crate::fft::{BLS12381Domain, FFTDomain};
use crate::KZG;

/// Computes the matrix-vector multiplication using explicit matrix representation
pub fn compute_matrix_vector_multiplication(
    coefficients: &[Scalar],
    tau_powers: &[G1Element],
) -> Vec<G1Element> {
    let d = coefficients.len();
    let mut result = vec![G1Element::zero(); d];
    let mut matrix = vec![vec![Scalar::zero(); d]; d];

    // Construct the matrix from coefficients
    let mut coeff = coefficients.to_vec();
    coeff.reverse();

    for i in 0..d {
        matrix[i] = coeff.clone();
        coeff.pop();
        coeff.insert(0, Scalar::zero());
    }

    // Perform matrix-vector multiplication
    for i in 0..d {
        for j in 0..d {
            result[i] += tau_powers[j].mul(matrix[i][j]);
        }
    }

    result
}

/// Computes the Lagrange interpolation polynomial for given roots and values
pub fn lagrange_interpolation(roots: &[Scalar], values: &[Scalar]) -> Vec<Scalar> {
    let n = roots.len();
    let mut coefficients = vec![Scalar::zero(); n];

    for i in 0..n {
        let mut basis_poly = vec![Scalar::from(1)];
        
        for j in 0..n {
            if i != j {
                let mut new_basis_poly = vec![Scalar::zero(); basis_poly.len() + 1];
                for k in 0..basis_poly.len() {
                    new_basis_poly[k] -= basis_poly[k] * roots[j];
                    new_basis_poly[k + 1] += basis_poly[k];
                }
                let denom = roots[i] - roots[j];
                for coeff in &mut new_basis_poly {
                    *coeff = (*coeff / denom).unwrap();
                }
                basis_poly = new_basis_poly;
            }
        }

        for (k, coeff) in basis_poly.iter().enumerate() {
            coefficients[k] += *coeff * values[i];
        }
    }

    coefficients
}


/// Multiplies a Toeplitz matrix with a vector using FFT
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
    domain.fft_in_place_g1(&mut padded_x);
    let y_vec = padded_x.clone();

    // Element-wise multiplication in the frequency domain
    let mut fft_result: Vec<G1Element> = v_vec.iter()
        .zip(y_vec.iter())
        .map(|(a, b)| b.mul(*a))
        .collect();
    
    // Inverse FFT to get the result in the time domain
    domain.ifft_in_place_g1(&mut fft_result);

    for i in 0..d {
        result[i] = fft_result[i];
    }

    result
}

/// Struct for KZG commitments using Tabular Domain Feist-Khovratovich technique
#[derive(Clone)]
pub struct KZGTabDFK {
    domain: BLS12381Domain,
    a: G1Element,
    g2_tau: G2Element,
    u_vec: Vec<G1Element>,
    l_vec: Vec<G1Element>,
    a_vec: Vec<G1Element>,
    tau_powers_g1: Vec<G1Element>,
}

impl KZG for KZGTabDFK {
    type G = G1Element;

    /// Creates a new KZGTabDFK instance with a random tau
    fn new(n: usize) -> FastCryptoResult<Self> {
        let domain = BLS12381Domain::new(n)?;
        let n_dom = domain.size();

        let tau = Scalar::rand(&mut thread_rng());
        let g2_tau = G2Element::generator().mul(tau);

        let g_tau_n = (0..n_dom).fold(G1Element::generator(), |acc, _| acc * tau);
        let a = g_tau_n - G1Element::generator();

        let mut a_vec = vec![G1Element::zero(); n_dom];
        let mut u_vec = vec![G1Element::zero(); n_dom];

        let tau_powers_g: Vec<Scalar> =
            iterate(Scalar::generator(), |g| g * tau).take(n_dom).collect();
        let tau_powers_g1: Vec<G1Element> = itertools::iterate(G1Element::generator(), |g| g.mul(tau))
            .take(n_dom)
            .collect();

        let g = G1Element::generator();
        let l_vec: Vec<G1Element> = domain
            .ifft(&tau_powers_g)
            .iter()
            .map(|s| g.mul(s))
            .collect();

        let mut omega_i = domain.element(0);
        for i in 0..n_dom {
            let denom = tau - omega_i;
            let a_i = a.mul(denom.inverse().unwrap());
            a_vec[i] = a_i;

            let l_i_minus_1 = l_vec[i] - G1Element::generator();
            let denom = tau - omega_i;
            let u_i = l_i_minus_1.mul(denom.inverse().unwrap());
            u_vec[i] = u_i;

            if i < n_dom - 1 {
                omega_i *= domain.element(1);
            }
        }

        Ok(Self {
            domain,
            g2_tau,
            a,
            u_vec,
            l_vec,
            a_vec,
            tau_powers_g1,
        })
    }

    /// Commits to a vector using the KZG commitment scheme
    fn commit(&self, v: &[Scalar]) -> G1Element {
        let mut padded_v = vec![Scalar::zero(); self.domain.size()];
        for i in 0..v.len() {
            padded_v[i] = v[i];
        }
        G1Element::multi_scalar_mul(&padded_v, &self.l_vec).unwrap()
    }

    /// Opens a KZG commitment at a specific index
    fn open(&self, v: &[Scalar], index: usize) -> G1Element {
        let mut open = G1Element::zero();
        for j in 0..v.len() {
            if j != index {
                let omega_i = self.domain.element(index);
                let omega_j = self.domain.element(j);

                let c_i = (omega_i - omega_j).inverse().unwrap();
                let c_j = (omega_j - omega_i).inverse().unwrap();

                let w_ij = self.a_vec[index].mul(c_i) + self.a_vec[j].mul(c_j);

                let omega_j_n = (omega_j / Scalar::from(v.len() as u128)).unwrap();
                let u_ij = w_ij.mul(omega_j_n);

                open += u_ij.mul(v[j]);
            }
        }
        open += self.u_vec[index].mul(v[index]);
        open
    }

    /// Opens a KZG commitment at multiple indices
    fn open_all(&self, v: &[Scalar], indices: Vec<usize>) -> Vec<G1Element> {
        let mut roots = vec![Scalar::zero(); v.len()];
        for i in 0..v.len() {
            roots[i] = self.domain.element(i);
        }
    
        let poly = lagrange_interpolation(&roots, &v);
    
        let degree = poly.len() - 1;
        let mut t = self.tau_powers_g1.clone();
        t.truncate(degree);
        t.reverse(); // [t^d-1], ... [1]
    
        let mut test_poly = poly.clone();
        test_poly.reverse();
        test_poly.truncate(degree);
        test_poly.reverse();
    
        let mut h = compute_matrix_vector_multiplication(&test_poly, &t);
    
        let mut result_eff = vec![G1Element::zero(); self.domain.size()];
        for i in 0..degree {
            result_eff[i] = h[i];
        }
    
        self.domain.fft_in_place_g1(&mut result_eff);
    
        result_eff
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
        index: usize,
        old_v_i: &Scalar,
        new_v_i: &Scalar,
    ) -> G1Element {
        *commitment + self.l_vec[index].mul(new_v_i - old_v_i)
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
        let omega_i = self.domain.element(index);
        let omega_j = self.domain.element(index_j);

        let c_i = (omega_i - omega_j).inverse().unwrap();
        let c_j = (omega_j - omega_i).inverse().unwrap();

        let w_ij = self.a_vec[index].mul(c_i) + self.a_vec[index_j].mul(c_j);

        let omega_j_n = (omega_j / Scalar::from(self.domain.size() as u128)).unwrap();
        let u_ij = w_ij.mul(omega_j_n);
        *open + u_ij.mul(new_v_j - old_v_j)
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
        let kzg = KZGTabDFK::new(n).unwrap();
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
        let kzg = KZGTabDFK::new(n).unwrap();
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
        let kzg = KZGTabDFK::new(n).unwrap();
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
        let n = 5;
        let kzg = KZGTabDFK::new(n).unwrap();
        let v: Vec<Scalar> = (0..n).map(|_| OtherScalar::rand(&mut rng)).collect();
        let commitment = kzg.commit(&v);
        let indices: Vec<usize> = (0..n).collect();
        let mut open_values = kzg.open_all(&v, indices.clone());

        open_values.truncate(n);

        for (i, open_value) in open_values.iter().enumerate() {
            let is_valid = kzg.verify(indices[i], &v[indices[i]], &commitment, open_value);
            assert!(is_valid, "Verification of the opening should succeed for index {}", i);
        }
    }

    #[test]
    fn test_mutliply_with_toeplitz() {
        let v_scalar = vec![
            Scalar::from(1), Scalar::from(2), Scalar::from(3), Scalar::from(4),
            Scalar::from(4), Scalar::from(4), Scalar::from(4), Scalar::from(4)
        ];

        let tau = Scalar::rand(&mut thread_rng());

        let tau_powers_g1: Vec<G1Element> = itertools::iterate(G1Element::generator(), |g| g.mul(tau))
            .take(8)
            .collect();

        let h_alin = multiply_toeplitz_with_v(&v_scalar, &tau_powers_g1);
        let h_mult = compute_matrix_vector_multiplication(&v_scalar, &tau_powers_g1);

        assert!(h_alin == h_mult, "Toeplitz multiplication mismatch.");
    }
}
