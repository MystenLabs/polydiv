use std::ops::{AddAssign, Div, MulAssign, SubAssign};
use std::ops::Mul;

use ark_bls12_381::{Bls12_381, Fr, G1Projective, G2Projective};
use ark_ec::{AffineRepr, CurveGroup, Group, pairing::Pairing};
use ark_ff::{Field, One, PrimeField, UniformRand, Zero};
use ark_poly::{DenseUVPolynomial, EvaluationDomain, GeneralEvaluationDomain, Polynomial};
use fastcrypto::error::{FastCryptoError, FastCryptoResult};
use fastcrypto::groups::bls12381::Scalar;
use fastcrypto::groups::Scalar as OtherScalar;
use fastcrypto::serde_helpers::ToFromByteArray;
use rand::rngs::OsRng;

use crate::KZG;

// Function to convert fastcrypto::Scalar to ark_bls12_381::Fr
fn scalar_to_fr(scalar: &Scalar) -> Fr {
    Fr::from_be_bytes_mod_order(&scalar.to_byte_array())
}

fn polynomial_division(dividend: &[Fr], divisor: &[Fr]) -> (Vec<Fr>, Vec<Fr>) {
    let mut quotient_coeffs = vec![Fr::zero(); dividend.len() - divisor.len() + 1];
    let mut remainder_coeffs = Vec::from(dividend);

    for i in (0..quotient_coeffs.len()).rev() {
        quotient_coeffs[i] = remainder_coeffs[i + divisor.len() - 1] / divisor.last().unwrap();
        for j in 0..divisor.len() {
            remainder_coeffs[i + j] -= quotient_coeffs[i] * divisor[j];
        }
    }

    // Remove leading zeros in the remainder
    while remainder_coeffs.len() > 1 && remainder_coeffs.last().unwrap().is_zero() {
        remainder_coeffs.pop();
    }

    (quotient_coeffs, remainder_coeffs)
}


pub struct KZGOriginal {
    domain: GeneralEvaluationDomain<Fr>,
    tau_powers_g1: Vec<G1Projective>,
    tau_powers_g2: Vec<G2Projective>,
}

impl KZGOriginal {
    pub fn new(n: usize) -> FastCryptoResult<Self> {
        let domain = GeneralEvaluationDomain::<Fr>::new(n).ok_or(FastCryptoError::InvalidInput)?;

        // Generate tau using a random scalar
        let tau = Fr::rand(&mut OsRng);

        // Compute g^tau^i for i = 0 to n-1 in G1
        let g1 = G1Projective::generator();
        let mut tau_powers_g1 = Vec::with_capacity(n);
        let mut current_power_g1 = g1;

        for _ in 0..n {
            tau_powers_g1.push(current_power_g1);
            current_power_g1 *= tau;
        }

        // Compute g^tau^i for i = 0 to n-1 in G2
        let g2 = G2Projective::generator();
        let mut tau_powers_g2 = Vec::with_capacity(n);
        let mut current_power_g2 = g2;

        for _ in 0..n {
            tau_powers_g2.push(current_power_g2);
            current_power_g2 *= tau;
        }

        Ok(Self { domain, tau_powers_g1, tau_powers_g2 })
    }
}

impl crate::KZG<Scalar, G1Projective> for KZGOriginal {
    fn commit(&self, v: &[Scalar]) -> G1Projective {
        let v_fr: Vec<Fr> = v.iter().map(scalar_to_fr).collect();
        let poly = self.domain.ifft(&v_fr);
        let mut commitment = G1Projective::zero();
        for (i, coeff) in poly.iter().enumerate() {
            commitment += self.tau_powers_g1[i].mul(*coeff);
        }
        commitment
    }

    fn open(&self, v: &[Scalar], index: usize) -> G1Projective {
        let v_fr: Vec<Fr> = v.iter().map(scalar_to_fr).collect();
        let poly = self.domain.ifft(&v_fr);
        println!("{:?}", poly);
        let v_index_fr = scalar_to_fr(&v[index]);
        let mut adjusted_poly = poly.clone();
        adjusted_poly[0] -= v_index_fr;
        println!("{:?}", adjusted_poly);
        let omega = self.domain.element(index);
        let divisor = [-omega, Fr::one()];
        let (quotient, remainder) = polynomial_division(&adjusted_poly, &divisor);

        let mut open_value = G1Projective::zero();
        for (i, coeff) in quotient.iter().enumerate() {
            open_value += self.tau_powers_g1[i].mul(*coeff);
        }

        // //check that the opening is correct without group operations
        //let rhs_poly = &divisor*&quotient;
        //println!("{}", adjusted_poly == rhs_poly);

        open_value
    }

    fn verify(&self, index: usize, v_i: &Scalar, commitment: &G1Projective, open_i: &G1Projective) -> bool {
        let v_i_fr = scalar_to_fr(v_i);
        let g = self.tau_powers_g1[0];
        let g_v_i = g.mul(v_i_fr);
        let lhs = *commitment - g_v_i;
        let g_tau = self.tau_powers_g2[1];
        let omega = self.domain.element(index);
        // let omega_inv = omega.inverse().unwrap();
        let g_omega = self.tau_powers_g2[0].mul(omega);
        let rhs = g_tau - g_omega;
        // Perform the pairing check e(lhs, g) == e(open_i, rhs)
        let lhs_pairing = Bls12_381::pairing(lhs.into_affine(), self.tau_powers_g2[0]);
        let rhs_pairing = Bls12_381::pairing(open_i.into_affine(), rhs.into_affine());

        lhs_pairing == rhs_pairing
    }

    fn update(&self, commitment: &mut G1Projective, index: usize, new_v_i: &Scalar) -> G1Projective {
        *commitment
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
        
        // Create a new KZGOriginal struct
        let n = 8; 
        let kzg = KZGOriginal::new(n).unwrap();

        // Create a random vector v
        let v: Vec<Scalar> = (0..n).map(|_| OtherScalar::rand(&mut rng)).collect();

        println!("{:?}",v);

        // Create a commitment
        let commitment = kzg.commit(&v);

        // Pick a random index to open
        let index = rng.gen_range(0..n);

        // Create an opening
        let open_value = kzg.open(&v, index);

        // Verify the opening
        let is_valid = kzg.verify(index, &v[index], &commitment, &open_value);

        // Assert that the verification passes
        assert!(is_valid, "Verification of the opening should succeed.");
    }
}
