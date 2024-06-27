use std::ops::{AddAssign, Div, MulAssign, SubAssign};
use std::ops::Mul;

use fastcrypto::error::FastCryptoResult;
use fastcrypto::groups::{GroupElement, Pairing, Scalar as OtherScalar};
use fastcrypto::groups::bls12381::{G1Element, G2Element, Scalar};
use fastcrypto::serde_helpers::ToFromByteArray;
use rand::thread_rng;

use crate::fft::{BLS12381Domain, FFTDomain};
use crate::KZG;

fn polynomial_division(dividend: &[Scalar], divisor: &[Scalar]) -> (Vec<Scalar>, Vec<Scalar>) {
    let mut quotient_coeffs = vec![Scalar::zero(); dividend.len() - divisor.len() + 1];
    let mut remainder_coeffs = Vec::from(dividend);

    for i in (0..quotient_coeffs.len()).rev() {
        quotient_coeffs[i] = (remainder_coeffs[i + divisor.len() - 1] / *divisor.last().unwrap()).unwrap();
        for j in 0..divisor.len() {
            remainder_coeffs[i + j] -= quotient_coeffs[i] * divisor[j];
        }
    }

    // Remove leading zeros in the remainder
    while remainder_coeffs.len() > 1 && remainder_coeffs.last().unwrap() == &Scalar::zero() {
        remainder_coeffs.pop();
    }

    (quotient_coeffs, remainder_coeffs)
}


pub struct KZGOriginal {
    domain: BLS12381Domain,
    tau_powers_g1: Vec<G1Element>,
    tau_powers_g2: Vec<G2Element>,
}

impl KZGOriginal {
    pub fn new(n: usize) -> FastCryptoResult<Self> {
        let domain = BLS12381Domain::new(n)?;

        // Generate tau using a random scalar
        let tau = Scalar::rand(&mut thread_rng());

        // Compute g^tau^i for i = 0 to n-1 in G1
        let g1 = G1Element::generator();
        let mut tau_powers_g1 = Vec::with_capacity(n);
        let mut current_power_g1 = g1;

        for _ in 0..n {
            tau_powers_g1.push(current_power_g1);
            current_power_g1 *= tau;
        }

        // Compute g^tau^i for i = 0 to n-1 in G2
        let g2 = G2Element::generator();
        let mut tau_powers_g2 = Vec::with_capacity(n);
        let mut current_power_g2 = g2;

        for _ in 0..n {
            tau_powers_g2.push(current_power_g2);
            current_power_g2 *= tau;
        }

        Ok(Self { domain, tau_powers_g1, tau_powers_g2 })
    }
}

impl KZG for KZGOriginal {

    type G = G1Element;


    fn commit(&self, v: &[Scalar]) -> G1Element {
        let poly = self.domain.ifft(&v);
        let mut commitment = G1Element::zero();
        for (i, coeff) in poly.iter().enumerate() {
            commitment += self.tau_powers_g1[i].mul(*coeff);
        }
        commitment
    }

    fn open(&self, v: &[Scalar], index: usize) -> G1Element {
        let poly = self.domain.ifft(&v);
        let mut adjusted_poly = poly.clone();
        adjusted_poly[0] -= &v[index];
        let omega = self.domain.element(index);
        let divisor = [-omega, Scalar::generator()];
        let (quotient, remainder) = polynomial_division(&adjusted_poly, &divisor);

        let mut open_value = G1Element::zero();
        for (i, coeff) in quotient.iter().enumerate() {
            open_value += self.tau_powers_g1[i].mul(*coeff);
        }

        // //check that the opening is correct without group operations
        //let rhs_poly = &divisor*&quotient;
        //println!("{}", adjusted_poly == rhs_poly);

        open_value
    }

    fn verify(&self, index: usize, v_i: &Scalar, commitment: &G1Element, open_i: &G1Element) -> bool {
        let g = self.tau_powers_g1[0];
        let g_v_i = g * v_i;
        let lhs = *commitment - g_v_i;
        let g_tau = self.tau_powers_g2[1];
        let omega = self.domain.element(index);
        // let omega_inv = omega.inverse().unwrap();
        let g_omega = self.tau_powers_g2[0].mul(omega);
        let rhs = g_tau - g_omega;
        // Perform the pairing check e(lhs, g) == e(open_i, rhs)
        let lhs_pairing = lhs.pairing(&self.tau_powers_g2[0]);
        let rhs_pairing = open_i.pairing(&rhs);

        lhs_pairing == rhs_pairing
    }

    fn update(&self, commitment: &mut G1Element, index: usize, new_v_i: &Scalar) -> G1Element {
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
