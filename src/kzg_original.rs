use ark_bls12_381::{Fr, G1Projective, G1Affine, G2Projective, G2Affine, Bls12_381};
use ark_ec::{pairing::Pairing, AffineRepr, CurveGroup, Group};
use ark_ff::{Field, UniformRand, PrimeField, BigInteger256, Zero, One};
use ark_poly::{univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain, Polynomial, DenseUVPolynomial};
use fastcrypto::error::{FastCryptoError, FastCryptoResult};
use fastcrypto::groups::bls12381::Scalar;
use fastcrypto::groups::Scalar as OtherScalar;
use rand::rngs::OsRng;
use std::ops::{AddAssign, SubAssign, MulAssign, Div};
use fastcrypto::serde_helpers::ToFromByteArray;
use std::ops::Mul;
use crate::KZG;

// Function to convert fastcrypto::Scalar to ark_bls12_381::Fr
fn scalar_to_fr(scalar: &Scalar) -> Fr {
    Fr::from_be_bytes_mod_order(&scalar.to_byte_array())
}

fn polynomial_division(dividend: &DensePolynomial<Fr>, divisor: &DensePolynomial<Fr>) -> (DensePolynomial<Fr>, DensePolynomial<Fr>) {
    let mut quotient_coeffs = vec![Fr::zero(); dividend.coeffs.len() - divisor.coeffs.len() + 1];
    let mut remainder_coeffs = dividend.coeffs.clone();

    for i in (0..quotient_coeffs.len()).rev() {
        quotient_coeffs[i] = remainder_coeffs[i + divisor.coeffs.len() - 1] / divisor.coeffs.last().unwrap();
        for j in 0..divisor.coeffs.len() {
            remainder_coeffs[i + j] -= quotient_coeffs[i] * divisor.coeffs[j];
        }
    }

    // Remove leading zeros in the remainder
    while remainder_coeffs.len() > 1 && remainder_coeffs.last().unwrap().is_zero() {
        remainder_coeffs.pop();
    }

    (DensePolynomial::from_coefficients_vec(quotient_coeffs), DensePolynomial::from_coefficients_vec(remainder_coeffs))
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
        let poly = DensePolynomial::from_coefficients_vec(self.domain.ifft(&v_fr));
        let mut commitment = G1Projective::zero();
        for (i, coeff) in poly.coeffs.iter().enumerate() {
            commitment += self.tau_powers_g1[i].mul(*coeff);
        }
        commitment
    }

    fn open(&self, v: &[Scalar], index: usize) -> G1Projective {
        let v_fr: Vec<Fr> = v.iter().map(scalar_to_fr).collect();
        let poly = DensePolynomial::from_coefficients_vec(self.domain.ifft(&v_fr));
        println!("{:?}", poly);
        let v_index_fr = scalar_to_fr(&v[index]);
        let mut adjusted_poly_coeffs = poly.coeffs.clone();
        adjusted_poly_coeffs[0] -= v_index_fr;
        let adjusted_poly = DensePolynomial::from_coefficients_vec(adjusted_poly_coeffs);
        println!("{:?}", adjusted_poly);
        let omega = self.domain.element(index);
        let divisor = DensePolynomial::from_coefficients_slice(&[-omega, Fr::one()]);
        let (quotient, remainder) = polynomial_division(&adjusted_poly, &divisor);

        let mut open_value = G1Projective::zero();
        for (i, coeff) in quotient.coeffs.iter().enumerate() {
            open_value += self.tau_powers_g1[i].mul(*coeff);
        }

        // //check that the opening is correct without group operations
        let rhs_poly = &divisor*&quotient;
        println!("{}", adjusted_poly == rhs_poly);

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
    use super::*;
    use fastcrypto::groups::bls12381::Scalar;
    use rand::Rng;  // Use `rand` crate directly
    use ark_ff::UniformRand;

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
