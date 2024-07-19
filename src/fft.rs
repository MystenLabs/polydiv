use ark_bls12_381::Fr;
use ark_ff::{BigInteger, PrimeField};
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
use fastcrypto::error::{FastCryptoError, FastCryptoResult};
use fastcrypto::groups::bls12381::{Scalar, G1Element};
use fastcrypto::serde_helpers::ToFromByteArray;
use fastcrypto::groups::GroupElement;
use std::ops::{Mul, Sub, Add, Neg};
use ark_ff::Field;

pub trait FFTDomain: Sized {
    type ScalarType;

    fn new(n: usize) -> FastCryptoResult<Self>;

    fn fft(&self, v: &[Self::ScalarType]) -> Vec<Self::ScalarType>;

    fn ifft(&self, v_hat: &[Self::ScalarType]) -> Vec<Self::ScalarType>;

    fn fft_in_place_g1(&self, v: &mut [G1Element]);

    fn ifft_in_place_g1(&self, v_hat: &mut [G1Element]);

    fn element(&self, index: usize) -> Self::ScalarType;

    fn size(&self) -> usize;
}

#[derive(Clone)]
pub struct BLS12381Domain {
    domain: GeneralEvaluationDomain<Fr>,
}

impl FFTDomain for BLS12381Domain {
    type ScalarType = Scalar;

    fn new(n: usize) -> FastCryptoResult<Self> {
        let domain = GeneralEvaluationDomain::<Fr>::new(n).ok_or(FastCryptoError::InvalidInput)?;
        Ok(Self { domain })
    }

    fn fft(&self, v: &[Scalar]) -> Vec<Scalar> {
        let scalars = v.iter().map(fastcrypto_to_arkworks).collect::<Vec<_>>();
        let result = self.domain
            .fft(&scalars)
            .iter()
            .map(arkworks_to_fastcrypto)
            .collect::<Vec<_>>();
        result
    }

    fn ifft(&self, v_hat: &[Scalar]) -> Vec<Scalar> {
        let scalars = v_hat.iter().map(fastcrypto_to_arkworks).collect::<Vec<_>>();
        let result = self.domain
            .ifft(&scalars)
            .iter()
            .map(arkworks_to_fastcrypto)
            .collect::<Vec<_>>();
        result
    }

    fn fft_in_place_g1(&self, v: &mut [G1Element]) {
        let mut padded_v = v.to_vec();
        if padded_v.len() < self.domain.size() {
            padded_v.resize(self.domain.size(), G1Element::zero());
        } else if padded_v.len() > self.domain.size() {
            padded_v.truncate(self.domain.size());
        }
        let root_of_unity = self.element(1);
        let mut result = fft_g1(self, &padded_v, &root_of_unity);
        v.copy_from_slice(&result[..v.len()]);
    }

    fn ifft_in_place_g1(&self, v: &mut [G1Element]) {
        let n = v.len();
        let root_of_unity = self.element(self.domain.size() - 1);
        let mut padded_v_hat = v.to_vec();
        if padded_v_hat.len() < self.domain.size() {
            padded_v_hat.resize(self.domain.size(), G1Element::zero());
        } else if padded_v_hat.len() > self.domain.size() {
            padded_v_hat.truncate(self.domain.size());
        }
        let mut result = fft_g1(self, &padded_v_hat, &root_of_unity);
        let n_fr = Fr::from(n as u64);
        let inv_n_fr = n_fr.inverse().unwrap();
        let inv_n_scalar = arkworks_to_fastcrypto(&inv_n_fr);

        for elem in result.iter_mut() {
            *elem = elem.mul(inv_n_scalar);
        }
        v.copy_from_slice(&result[..v.len()]);
    }

    fn element(&self, index: usize) -> Scalar {
        arkworks_to_fastcrypto(&self.domain.element(index))
    }

    fn size(&self) -> usize {
        self.domain.size()
    }
}

/// Converts a fastcrypto Scalar to an arkworks Fr
fn fastcrypto_to_arkworks(s: &Scalar) -> Fr {
    let bytes = s.to_byte_array();
    Fr::from_be_bytes_mod_order(&bytes)
}

/// Converts an arkworks Fr to a fastcrypto Scalar
fn arkworks_to_fastcrypto(f: &Fr) -> Scalar {
    let bytes: [u8; 32] = f.into_bigint().to_bytes_be().try_into().unwrap();
    Scalar::from_byte_array(&bytes).unwrap()
}

/// FFT function for G1 elements
fn fft_g1(domain: &BLS12381Domain, v: &[G1Element], root_of_unity: &Scalar) -> Vec<G1Element> {
    let n = v.len();
    if n <= 1 {
        return v.to_vec();
    }

    let half_n = n / 2;
    let even: Vec<G1Element> = v.iter().step_by(2).cloned().collect();
    let odd: Vec<G1Element> = v.iter().skip(1).step_by(2).cloned().collect();

    let even_fft = fft_g1(domain, &even, &root_of_unity.mul(root_of_unity));
    let odd_fft = fft_g1(domain, &odd, &root_of_unity.mul(root_of_unity));

    let mut omega = Scalar::from(1);
    let mut result = vec![G1Element::zero(); n];

    for i in 0..half_n {
        let t = odd_fft[i].mul(&omega);
        result[i] = even_fft[i] + t;
        result[i + half_n] = even_fft[i] - t;
        omega *= root_of_unity;
    }

    result
}

#[cfg(test)]
mod tests {
    use fastcrypto::groups::bls12381::{G1Element, Scalar};
    use fastcrypto::groups::{GroupElement, Scalar as OtherScalar};
    use rand::thread_rng;
    use crate::fft::{BLS12381Domain, FFTDomain};
    use std::ops::Mul;

    #[test]
    fn test_fft() {
        let domain = BLS12381Domain::new(4).unwrap();
        let v = vec![
            Scalar::generator(),
            Scalar::zero(),
            Scalar::generator(),
            Scalar::zero(),
        ];
        let v_hat = domain.fft(&v);
        let v_prime = domain.ifft(&v_hat);
        assert_eq!(v, v_prime);
    }

    #[test]
    fn test_fft_g1() {
        let mut rng = rand::thread_rng();
        let n = 8;
        let v_scalar: Vec<Scalar> = (0..n).map(|_| OtherScalar::rand(&mut rng)).collect();
        let domain = BLS12381Domain::new(n).unwrap();
        let g = G1Element::generator();
        let v_group_elements: Vec<G1Element> = v_scalar.iter().map(|x| g.mul(x)).collect();
        let mut v_fft_g = v_group_elements.clone();
        domain.fft_in_place_g1(&mut v_fft_g);
        let v_fft = domain.fft(&v_scalar);
        let expected_v_fft_g: Vec<G1Element> = v_fft.iter().map(|x| g.mul(x)).collect();
        assert_eq!(v_fft_g, expected_v_fft_g);
    }

    #[test]
    fn test_ifft_g1() {
        let domain = BLS12381Domain::new(8).unwrap();
        let g = G1Element::generator();
        let v_scalar = vec![
            Scalar::from(1), Scalar::from(2), Scalar::from(3), Scalar::from(4),
            Scalar::from(4), Scalar::from(4), Scalar::from(4), Scalar::from(4)
        ];
        let v_group_elements: Vec<G1Element> = v_scalar.iter().map(|x| g.mul(x)).collect();
        let mut v_ifft_g = v_group_elements.clone();
        domain.ifft_in_place_g1(&mut v_ifft_g);
        let v_ifft = domain.ifft(&v_scalar);
        let expected_v_ifft_g: Vec<G1Element> = v_ifft.iter().map(|x| g.mul(x)).collect();
        assert_eq!(v_ifft_g, expected_v_ifft_g);
    }
}
