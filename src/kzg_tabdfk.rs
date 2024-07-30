use std::ops::Mul;

use fastcrypto::error::FastCryptoResult;
use fastcrypto::groups::bls12381::{G1Element, G2Element, Scalar};
use fastcrypto::groups::{GroupElement, MultiScalarMul, Pairing, Scalar as OtherScalar};
use itertools::iterate;
use rand::thread_rng;

use crate::fft::{BLS12381Domain, FFTDomain};
use crate::KZG;

pub fn build_circulant(polynomial: &[Scalar], size: usize) -> Vec<Scalar> {
    let mut circulant = vec![Scalar::zero(); 2 * size];
    let coeffs = polynomial;

    if size == coeffs.len() - 1 {
        circulant[0] = *coeffs.last().unwrap();
        circulant[size] = *coeffs.last().unwrap();
        circulant[size + 1..size + 1 + coeffs.len() - 2]
            .copy_from_slice(&coeffs[1..coeffs.len() - 1]);
    } else {
        circulant[size + 1..size + 1 + coeffs.len() - 1].copy_from_slice(&coeffs[1..]);
    }

    circulant
}

pub fn multiply_toeplitz_with_v(
    polynomial: &[Scalar],
    v: &[G1Element],
    size: usize,
) -> Vec<G1Element> {
    let m = polynomial.len() - 1;
    let size = std::cmp::max(size, m);

    let domain = BLS12381Domain::new(2 * size).unwrap();

    let size = domain.size() / 2;
    let mut circulant = build_circulant(polynomial, size);

    let mut tmp: Vec<G1Element> = Vec::with_capacity(domain.size());

    for _ in 0..(size - v.len()) {
        tmp.push(G1Element::zero());
    }

    for i in v.iter().rev() {
        tmp.push(*i);
    }

    tmp.resize(domain.size(), G1Element::zero());
    domain.fft_in_place_group(&mut tmp);
    let circulant_fft = domain.fft(&mut circulant);

    for (i, j) in tmp.iter_mut().zip(circulant_fft.iter()) {
        *i = i.mul(*j);
    }

    domain.ifft_in_place_group(&mut tmp);
    let mut result = vec![G1Element::zero(); size];
    for i in 0..size {
        result[i] = tmp[i];
    }
    result
}

#[derive(Clone)]
pub struct KZGTabDFK {
    domain: BLS12381Domain,
    g2_tau: G2Element,
    u_vec: Vec<G1Element>,
    l_vec: Vec<G1Element>,
    a_vec: Vec<G1Element>,
    tau_powers_g1: Vec<G1Element>,
}

impl KZG for KZGTabDFK {
    type G = G1Element;

    fn new(n: usize) -> FastCryptoResult<Self> {
        let domain = BLS12381Domain::new(n)?;
        let n_dom = domain.size();

        let tau = Scalar::rand(&mut thread_rng());
        let g2_tau = G2Element::generator().mul(tau);

        let g_tau_n = (0..n_dom).fold(G1Element::generator(), |acc, _| acc * tau);
        let a = g_tau_n - G1Element::generator();

        let mut a_vec = vec![G1Element::zero(); n_dom];
        let mut u_vec = vec![G1Element::zero(); n_dom];

        let tau_powers_g: Vec<Scalar> = iterate(Scalar::generator(), |g| g * tau)
            .take(n_dom)
            .collect();
        let tau_powers_g1: Vec<G1Element> =
            itertools::iterate(G1Element::generator(), |g| g.mul(tau))
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
            u_vec,
            l_vec,
            a_vec,
            tau_powers_g1,
        })
    }

    fn commit(&self, v: &[Scalar]) -> G1Element {
        let mut padded_v = vec![Scalar::zero(); self.domain.size()];
        for i in 0..v.len() {
            padded_v[i] = v[i];
        }
        G1Element::multi_scalar_mul(&padded_v, &self.l_vec).unwrap()
    }

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

    fn open_all(&self, v: &[Scalar], indices: Vec<usize>) -> Vec<G1Element> {
        // Jonas: Why are the indices not used here?

        let domain = &self.domain;

        let poly = domain.ifft(&v);
        let poly_degree = poly.len() - 1;
        let mut t = self.tau_powers_g1.clone();
        t.truncate(poly_degree);

        let mut h = multiply_toeplitz_with_v(&poly, &t, domain.size());
        domain.fft_in_place_group(&mut h);

        h
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
        let n = 8;
        let kzg = KZGTabDFK::new(n).unwrap();
        let v: Vec<Scalar> = (0..n).map(|_| OtherScalar::rand(&mut rng)).collect();
        let commitment = kzg.commit(&v);
        let indices: Vec<usize> = (0..n).collect();
        let mut open_values = kzg.open_all(&v, indices.clone());

        open_values.truncate(n);

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
