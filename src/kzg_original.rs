use fastcrypto::error::{FastCryptoError, FastCryptoResult};
use fastcrypto::groups::bls12381::{G1Element, G2Element, Scalar};
use fastcrypto::groups::{GroupElement, MultiScalarMul, Pairing, Scalar as OtherScalar};
use fastcrypto::serde_helpers::ToFromByteArray;
use rand::thread_rng;

use crate::fft::{BLS12381Domain, FFTDomain};
use crate::KZG;

fn polynomial_division(
    dividend: &[Scalar],
    divisor: &[Scalar],
) -> FastCryptoResult<(Vec<Scalar>, Vec<Scalar>)> {
    let mut remainder = Vec::from(dividend);

    let divisor_leading_term_inverse = divisor
        .last()
        .ok_or(FastCryptoError::InvalidInput)?
        .inverse()?;

    let quotient_size = dividend.len() - divisor.len() + 1;

    let mut quotient: Vec<Scalar> = (0..quotient_size)
        .rev()
        .map(|i| {
            let q_i = remainder[i + divisor.len() - 1] * divisor_leading_term_inverse;
            for j in 0..divisor.len() {
                remainder[i + j] -= q_i * divisor[j];
            }
            q_i
        })
        .collect();
    quotient.reverse();

    // Remove leading zeros in the remainder
    while remainder.len() > 1 && remainder.last().unwrap() == &Scalar::zero() {
        remainder.pop();
    }

    Ok((quotient, remainder))
}

#[derive(Clone)]
pub struct KZGOriginal {
    domain: BLS12381Domain,
    tau_powers_g1: Vec<G1Element>,
    g2_tau: G2Element,
}

impl KZG for KZGOriginal {
    // Uses the BLS12-381 construction
    type G = G1Element;

    fn new(n: usize) -> FastCryptoResult<Self> {
        let domain = BLS12381Domain::new(n)?;

        // Generate tau using a random scalar
        let tau = fastcrypto::groups::bls12381::Scalar::rand(&mut thread_rng());

        // Compute g^tau^i for i = 0 to n-1 in G1
        let tau_powers_g1: Vec<G1Element> = itertools::iterate(G1Element::generator(), |g| g * tau)
            .take(n)
            .collect();

        let g2_tau = G2Element::generator() * tau;

        Ok(Self {
            domain,
            tau_powers_g1,
            g2_tau,
        })
    }

    fn commit(&self, v: &[Scalar]) -> G1Element {
        let poly = self.domain.ifft(&v);
        G1Element::multi_scalar_mul(poly.as_slice(), &self.tau_powers_g1).unwrap()
    }

    fn open(&self, v: &[Scalar], index: usize) -> G1Element {
        let mut poly = self.domain.ifft(&v);
        poly[0] -= &v[index];

        let divisor = [-self.domain.element(index), Scalar::generator()];
        let (quotient, _) = polynomial_division(&poly, &divisor).unwrap();

        G1Element::multi_scalar_mul(&quotient, &self.tau_powers_g1[..quotient.len()]).unwrap()
    }

    fn open_all(&self, v: &[Scalar], indices: Vec<usize>) -> Vec<G1Element> {
        indices.into_iter().map(|i| self.open(v, i)).collect()
    }

    fn verify(
        &self,
        index: usize,
        v_i: &Scalar,
        commitment: &G1Element,
        open_i: &G1Element,
    ) -> bool {
        let lhs = *commitment - G1Element::generator() * v_i;

        let rhs = self.g2_tau - G2Element::generator() * self.domain.element(index);

        // Perform the pairing check e(lhs, g) == e(open_i, rhs)
        let lhs_pairing = lhs.pairing(&G2Element::generator());
        let rhs_pairing = open_i.pairing(&rhs);

        lhs_pairing == rhs_pairing
    }

    fn update(
        &self,
        commitment: &mut G1Element,
        index: usize,
        old_v_i: &Scalar,
        new_v_i: &Scalar,
    ) -> G1Element {
        *commitment
    }

    fn update_open_i(
        &self,
        open: &mut G1Element,
        index: usize,
        old_v_i: &Scalar,
        new_v_i: &Scalar,
    ) -> G1Element {
        *open
    }

    fn update_open_j(
        &self,
        open: &mut G1Element,
        index: usize,
        index_j: usize,
        old_v_j: &Scalar,
        new_v_j: &Scalar,
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

        // Create a new KZGOriginal struct
        let n = 8;
        let kzg = KZGOriginal::new(n).unwrap();

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
    fn test_kzg_commit_open_all() {
        let mut rng = rand::thread_rng();

        // Create a new KZGOriginal struct
        let n = 8;
        let kzg = KZGOriginal::new(n).unwrap();

        // Create a random vector v
        let v: Vec<Scalar> = (0..n).map(|_| OtherScalar::rand(&mut rng)).collect();

        println!("{:?}", v);

        // Create a commitment
        let commitment = kzg.commit(&v);

        // Create array with all indices
        let indices: Vec<usize> = (0..n).collect();

        // Create an opening
        let open_values = kzg.open_all(&v, indices.clone());

        // Verify all openings
        for (i, open_value) in open_values.iter().enumerate() {
            let is_valid = kzg.verify(indices[i], &v[indices[i]], &commitment, open_value);
            assert!(is_valid, "Verification of the opening should succeed for index {}", i);
        }
    }
}
