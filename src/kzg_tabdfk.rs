use fastcrypto::error::{FastCryptoError, FastCryptoResult};
use fastcrypto::groups::bls12381::{G1Element, G2Element, Scalar};
use fastcrypto::groups::{GroupElement, MultiScalarMul, Pairing, Scalar as OtherScalar};
use rand::thread_rng;
use std::ops::Mul;

use crate::fft::{BLS12381Domain, FFTDomain};
use crate::KZG;


pub struct KZGTabDFK {
    domain: BLS12381Domain,
    tau_powers_g1: Vec<G1Element>,
    tau_powers_g2: Vec<G2Element>,
    a:G1Element, 
    u_vec:Vec<G1Element>,
    l_vec:Vec<G1Element>,
    a_vec:Vec<G1Element>,
}

impl KZGTabDFK {
    pub fn new(n: usize) -> FastCryptoResult<Self> {
        let domain = BLS12381Domain::new(n)?;

        // Generate tau using a random scalar
        let tau = Scalar::rand(&mut thread_rng());

        // Compute g^tau^i for i = 0 to n-1 in G1
        let tau_powers_g1: Vec<G1Element> = itertools::iterate(G1Element::generator(), |g| g * tau)
            .take(n)
            .collect();

        // Compute g^tau^i for i = 0 to n-1 in G2
        let tau_powers_g2: Vec<G2Element> = itertools::iterate(G2Element::generator(), |g| g * tau)
            .take(n)
            .collect();
        
        //Compute a = g^{A(tau)} where A = X^n-1
        let a = tau_powers_g1[n-1] - G1Element::generator();

        // Compute a_i = g^{A(tau)/tau - omega^i}
        let mut a_vec = vec![G1Element::zero(); n];
        for i in 0..n {
            let omega_i = domain.element(i);
            let denom = tau - omega_i;
            let a_i = a.mul(denom.inverse().unwrap());
            a_vec[i] = a_i;
        }

        //Compute l_i = g^L_i(tau)

        let mut l_vec = vec![G1Element::zero(); n];
        for i in 0..n {
            let mut num = Scalar::from(1u128);
            let mut denom = Scalar::from(1u128);
            let omega_i = domain.element(i);

            for j in 0..n {
                if i != j {
                    let omega_j = domain.element(j);
                    num *= tau - omega_j;
                    denom *= omega_i - omega_j;
                }
            }

            let l_i = num * denom.inverse().unwrap();
            l_vec[i] = G1Element::generator().mul(l_i);
        }

        //Compute u_i = g^{(L_i(tau) - 1)/(tau-omega^i)}
        let mut u_vec = vec![G1Element::zero(); n];
        for i in 0..n {
            let omega_i = domain.element(i);
            let l_i = l_vec[i];
            let l_i_minus_1 = l_i - tau_powers_g1[0];
            let denom = tau - omega_i;
            let u_i = l_i_minus_1.mul(denom.inverse().unwrap());
            u_vec[i] = u_i;
        }


        Ok(Self {
            domain,
            tau_powers_g1,
            tau_powers_g2,
            a,
            u_vec,
            l_vec,
            a_vec,
        })
    }
}

impl KZG for KZGTabDFK {
    // Uses the BLS12-381 construction
    type G = G1Element;

    fn commit(&self, v: &[Scalar]) -> G1Element {
        let poly = self.domain.ifft(&v);
        G1Element::multi_scalar_mul(&poly, &self.tau_powers_g1).unwrap()
    }

    fn open(&self, v: &[Scalar], index: usize) -> G1Element {
        let mut poly = self.domain.ifft(&v);
        let mut quotient_coeffs: Vec<Scalar> = vec![Scalar::zero(); poly.len()-1];

        quotient_coeffs[poly.len() - 2] = poly[poly.len() - 1];

        for j in (0..poly.len() - 2).rev() {
            quotient_coeffs[j] = poly[j + 1] + quotient_coeffs[j + 1] * self.domain.element(index);
        }
        G1Element::multi_scalar_mul(&quotient_coeffs, &self.tau_powers_g1[..quotient_coeffs.len()]).unwrap()

    }

    fn verify(
        &self,
        index: usize,
        v_i: &Scalar,
        commitment: &G1Element,
        open_i: &G1Element,
    ) -> bool {
        let lhs = *commitment - self.tau_powers_g1[0] * v_i;

        let rhs = self.tau_powers_g2[1] - self.tau_powers_g2[0] * self.domain.element(index);

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

        // Create a new KZGTabDFK struct
        let n = 8;
        let kzg = KZGTabDFK::new(n).unwrap();

        // Create a random vector v
        let v: Vec<Scalar> = (0..n).map(|_| OtherScalar::rand(&mut rng)).collect();

        println!("{:?}", v);

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
