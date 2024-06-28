use std::ops::Mul;

use fastcrypto::error::FastCryptoResult;
use fastcrypto::groups::bls12381::{G1Element, G2Element, Scalar};
use fastcrypto::groups::{GroupElement, MultiScalarMul, Pairing, Scalar as OtherScalar};
use rand::thread_rng;

use crate::fft::{BLS12381Domain, FFTDomain};
use crate::KZG;

pub struct KZGTabDFK {
    domain: BLS12381Domain,
    tau_powers_g1: Vec<G1Element>,
    tau_powers_g2: Vec<G2Element>,

    // Jonas: Is this variable necessary?
    a: G1Element,

    u_vec: Vec<G1Element>,
    l_vec: Vec<G1Element>,
    a_vec: Vec<G1Element>,
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
        let g_tau_n = tau_powers_g1[n - 1].mul(tau); // g^{tau^n}
        let g = G1Element::generator();
        let a = g_tau_n - g;

        // Compute a_i = g^{A(tau)/tau - omega^i}
        let mut a_vec = vec![G1Element::zero(); n];

        //Compute l_i = g^L_i(tau)
        let mut l_vec = vec![G1Element::zero(); n];

        //Compute u_i = g^{(L_i(tau) - 1)/(tau-omega^i)}
        let mut u_vec = vec![G1Element::zero(); n];

        let mut omega_i = domain.element(0);
        for i in 0..n {
            // Compute a_i
            let denom = tau - omega_i;
            let a_i = a.mul(denom.inverse().unwrap());
            a_vec[i] = a_i;

            // Compute l_i
            let mut num = Scalar::from(1u128);
            let mut denom = Scalar::from(1u128);
            for j in 0..n {
                if i != j {
                    let omega_j = domain.element(j);
                    num *= tau - omega_j;
                    denom *= omega_i - omega_j;
                }
            }
            let l_i = G1Element::generator().mul(num * denom.inverse().unwrap());
            l_vec[i] = l_i;

            // Compute u_i
            let l_i_minus_1 = l_i - tau_powers_g1[0];
            let denom = tau - omega_i;
            let u_i = l_i_minus_1.mul(denom.inverse().unwrap());
            u_vec[i] = u_i;

            // Update omega_i
            if i < n - 1 {
                omega_i *= domain.element(1);
            }
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
        G1Element::multi_scalar_mul(&v, &self.l_vec).unwrap()
    }

    fn open(&self, v: &[Scalar], index: usize) -> G1Element {
        let mut uij_vec = vec![G1Element::zero(); v.len()];
        for j in 0..v.len() {
            if j != index {
                let omega_i = self.domain.element(index);
                let omega_j = self.domain.element(j);

                // Compute c_i and c_j
                let c_i = (omega_i - omega_j).inverse().unwrap();
                let c_j = (omega_j - omega_i).inverse().unwrap();

                // Compute w_ij = a_i^c_i * a_j^c_j
                let w_ij = self.a_vec[index].mul(c_i) + self.a_vec[j].mul(c_j);

                // Compute u_ij = w_ij^{omega_j / n}
                let omega_j_n = (omega_j / Scalar::from(v.len() as u128)).unwrap();
                let u_ij = w_ij.mul(omega_j_n);

                uij_vec[j] = u_ij;
            }
            uij_vec[index] = self.u_vec[index];
        }
        let mut open = G1Element::zero();

        for (v_i, uij_i) in v.iter().zip(uij_vec.iter()) {
            open += uij_i.mul(*v_i);
        }

        open
        // G1Element::multi_scalar_mul(&v, &uij_vec[..v.len()]).unwrap()
        // let mut poly = self.domain.ifft(&v);
        // let mut quotient_coeffs: Vec<Scalar> = vec![Scalar::zero(); poly.len()-1];

        // quotient_coeffs[poly.len() - 2] = poly[poly.len() - 1];

        // for j in (0..poly.len() - 2).rev() {
        //     quotient_coeffs[j] = poly[j + 1] + quotient_coeffs[j + 1] * self.domain.element(index);
        // }
        // G1Element::multi_scalar_mul(&quotient_coeffs, &self.tau_powers_g1[..quotient_coeffs.len()]).unwrap()
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

        // Compute c_i and c_j
        let c_i = (omega_i - omega_j).inverse().unwrap();
        let c_j = (omega_j - omega_i).inverse().unwrap();

        // Compute w_ij = a_i^c_i * a_j^c_j
        let w_ij = self.a_vec[index].mul(c_i) + self.a_vec[index_j].mul(c_j);

        // Compute u_ij = w_ij^{omega_j / n}
        let omega_j_n = (omega_j / Scalar::from(self.tau_powers_g1.len() as u128)).unwrap();
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

    #[test]
    fn test_kzg_commit_open_update_i_verify() {
        let mut rng = rand::thread_rng();

        // Create a new KZGTabDFK struct
        let n = 8;
        let kzg = KZGTabDFK::new(n).unwrap();

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

        // Create a new KZGTabDFK struct
        let n = 8;
        let kzg = KZGTabDFK::new(n).unwrap();

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
