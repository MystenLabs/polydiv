use std::ops::Mul;

use fastcrypto::error::FastCryptoResult;
use fastcrypto::groups::bls12381::{G1Element, G2Element, Scalar};
use fastcrypto::groups::{GroupElement, MultiScalarMul, Pairing, Scalar as OtherScalar};
use fastcrypto::serde_helpers::ToFromByteArray;
use itertools::{iterate, Itertools};
use rand::thread_rng;

use crate::fft::{BLS12381Domain, FFTDomain};
use crate::KZG;

#[derive(Clone)]
pub struct KZGDeriv {
    n: usize,
    omega: Scalar,
    g2_tau: G2Element,
    w_vec: Vec<G1Element>,
    u_vec: Vec<G1Element>,
}

impl KZGDeriv {
    fn element(&self, index: usize) -> Scalar {
        if index == 0 {
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
    // Uses the BLS12-381 construction
    type G = G1Element;

    fn new(n: usize) -> FastCryptoResult<Self> {
        let domain = BLS12381Domain::new(n)?;

        // Generate tau using a random scalar
        let tau = Scalar::rand(&mut thread_rng());
        let g2_tau = G2Element::generator() * tau;

        // Compute tau^i for i = 0 to n-1
        let tau_powers_g1: Vec<Scalar> =
            iterate(Scalar::generator(), |g| g * tau).take(n).collect();

        // Compute w_vec and u_i = g^{(L_i(tau) - 1)/(tau-omega^i)}
        let w_vec: Vec<G1Element> = domain
            .ifft(&tau_powers_g1)
            .iter()
            .map(|s| G1Element::generator() * s)
            .collect();

        let mut omega_i = Scalar::generator();
        let mut u_vec = (0..n)
            .map(|i| {
                // Compute u_i
                let l_i_minus_1 = w_vec[i] - G1Element::generator();
                let denom = tau - omega_i;
                let u_i = (l_i_minus_1 / denom).unwrap();

                // Update omega_i
                if i < n - 1 {
                    omega_i *= domain.element(1);
                }
                u_i
            })
            .collect();

        let omega = domain.element(1);

        Ok(Self {
            n,
            omega,
            g2_tau,
            w_vec,
            u_vec,
        })
    }

    fn commit(&self, v: &[Scalar]) -> G1Element {
        G1Element::multi_scalar_mul(&v, &self.w_vec).unwrap()
    }

    fn open(&self, v: &[Scalar], index: usize) -> G1Element {
        let mut scalars = vec![Scalar::zero(); v.len()];
        let mut v_prime = Scalar::zero();

        // Intialize omega_j = 1 and omega_j_minus_i = omega_i^-1
        let omega_i = self.element(index);
        let mut omega_j = Scalar::generator();
        let mut omega_j_minus_i = self.element(self.n - index);

        for j in 0..v.len() {
            if j != index {
                let diff_inverse = (omega_i - omega_j).inverse().unwrap();
                v_prime += v[j] * omega_j_minus_i * diff_inverse;
                scalars[j] = (v[index] - v[j]) * diff_inverse;
            } else {
                v_prime += v[j]
                    * (Scalar::from((v.len() - 1) as u128) / (Scalar::from(2u128) * omega_i))
                        .unwrap();
            }

            if j < v.len() - 1 {
                omega_j *= self.omega;
                omega_j_minus_i *= self.omega;
            }
        }

        scalars[index] = v_prime;

        G1Element::multi_scalar_mul(&scalars, &self.w_vec[..scalars.len()]).unwrap()
    }

    fn verify(
        &self,
        index: usize,
        v_i: &Scalar,
        commitment: &G1Element,
        open_i: &G1Element,
    ) -> bool {
        let lhs = *commitment - G1Element::generator() * v_i;
        let rhs = self.g2_tau - G2Element::generator() * self.element(index);

        // Perform the pairing check
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

        // Create a new KZGDeriv struct
        let n = 8;
        let kzg = KZGDeriv::new(n).unwrap();

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

    #[test]
    fn test_kzg_commit_open_update_i_verify() {
        let mut rng = rand::thread_rng();

        // Create a new KZGDeriv struct
        let n = 8;
        let kzg = KZGDeriv::new(n).unwrap();

        // Create a random vector v
        let v: Vec<Scalar> = (0..n).map(|_| OtherScalar::rand(&mut rng)).collect();

        println!("{:?}", v);

        // Create a commitment
        let mut commitment = kzg.commit(&v);

        // Pick a random index to open
        let index = rng.gen_range(0..n);

        // Create an opening
        let mut open_value = kzg.open(&v, index);

        // Set a new valie for v_i
        let new_v_index = Scalar::rand(&mut rng);

        //Update the commitment
        let new_commitment = kzg.update(&mut commitment, index, &v[index], &new_v_index);

        //Update the opening
        let new_opening = kzg.update_open_i(&mut open_value, index, &v[index], &new_v_index);

        //Verify the updated opening
        let is_valid = kzg.verify(index, &new_v_index, &new_commitment, &new_opening);

        // Assert that the verification passes
        assert!(
            is_valid,
            "Verification of the opening after updating should succeed."
        );
    }

    #[test]
    fn test_kzg_commit_open_update_j_verify() {
        let mut rng = rand::thread_rng();

        // Create a new KZGDeriv struct
        let n = 8;
        let kzg = KZGDeriv::new(n).unwrap();

        // Create a random vector v
        let v: Vec<Scalar> = (0..n).map(|_| OtherScalar::rand(&mut rng)).collect();

        println!("{:?}", v);

        // Create a commitment
        let mut commitment = kzg.commit(&v);

        // Pick a random index to open
        let index = rng.gen_range(0..n);

        // Create an opening
        let mut open_value = kzg.open(&v, index);

        // Pick a new index to update¨
        let mut index_j;
        loop {
            index_j = rng.gen_range(0..n);
            if index_j != index {
                break;
            }
        }

        // Set a new value for v_i
        let new_v_index_j = Scalar::rand(&mut rng);

        //Update the commitment
        let new_commitment = kzg.update(&mut commitment, index_j, &v[index_j], &new_v_index_j);

        //Update the opening
        let new_opening =
            kzg.update_open_j(&mut open_value, index, index_j, &v[index_j], &new_v_index_j);

        //Verify the updated opening
        let is_valid = kzg.verify(index, &v[index], &new_commitment, &new_opening);

        // Assert that the verification passes
        assert!(
            is_valid,
            "Verification of the opening after updating j's value should succeed."
        );
    }
}
