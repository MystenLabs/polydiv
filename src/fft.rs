use ark_bls12_381::Fr;
use ark_ff::{BigInteger, PrimeField};
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
use fastcrypto::error::{FastCryptoError, FastCryptoResult};
use fastcrypto::groups::bls12381::{Scalar, G1Element};
use fastcrypto::serde_helpers::ToFromByteArray;
use ark_ff::Field;
use std::ops::Mul;
use fastcrypto::groups::GroupElement;


pub trait FFTDomain: Sized {
    type ScalarType;

    /// Create a new domain with at least `n` elements. Fails
    fn new(n: usize) -> FastCryptoResult<Self>;

    /// Compute the FFT of a vector of field elements. If the input is smaller than
    /// the domain size, it is padded with zeros. If it is larger, the input is truncated.
    fn fft(&self, v: &[Self::ScalarType]) -> Vec<Self::ScalarType>;

    /// Compute the inverse FFT of a vector of field elements. If the input is smaller than
    /// the domain size, it is padded with zeros. If it is larger, the input is truncated.
    fn ifft(&self, v_hat: &[Self::ScalarType]) -> Vec<Self::ScalarType>;

    /// Compute the FFT in place for G1 elements
    fn fft_in_place_g1(&self, v: &mut [G1Element]);

    /// Compute the inverse FFT in place for G1 elements
    fn ifft_in_place_g1(&self, v_hat: &mut [G1Element]);

    /// Get the i'th power of the root of unity.
    fn element(&self, index: usize) -> Self::ScalarType;

    /// Get the size of the domain.
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
        
        // println!("FFT result: {:?}", result);
        result
    }

    fn ifft(&self, v_hat: &[Scalar]) -> Vec<Scalar> {
        let scalars = v_hat.iter().map(fastcrypto_to_arkworks).collect::<Vec<_>>();
        let result = self.domain
            .ifft(&scalars)
            .iter()
            .map(arkworks_to_fastcrypto)
            .collect::<Vec<_>>();
        // println!("Inverse FFT result: {:?}", result); // Debug print
        result

    }

    fn fft_in_place_g1(&self, v: &mut [G1Element]) {
        let n = self.domain.size();
        println!("The domain size here is {:?}", n);
        let mut padded_v = v.to_vec();
        if padded_v.len() < n {
            padded_v.resize(n, G1Element::zero());
            println!("I was here FFT");
        } else if padded_v.len() > n {
            padded_v.truncate(n);
        }

        let root_of_unity = self.domain.element(1);
        let root_of_unity_scalar = arkworks_to_fastcrypto(&root_of_unity);
        fft_in_place_g1_recursive(&mut padded_v, &root_of_unity_scalar, n);

        v.copy_from_slice(&padded_v[..v.len()]);
    }

    fn ifft_in_place_g1(&self, v_hat: &mut [G1Element]) {
        let n = self.domain.size();
        let mut padded_v_hat = v_hat.to_vec();
        if padded_v_hat.len() < n {
            padded_v_hat.resize(n, G1Element::zero());
            println!("I was here iFFT");
        } else if padded_v_hat.len() > n {
            padded_v_hat.truncate(n);
        }

        let root_of_unity = self.domain.element(n - 1); // Inverse root of unity
        let root_of_unity_scalar = arkworks_to_fastcrypto(&root_of_unity);
        fft_in_place_g1_recursive(&mut padded_v_hat, &root_of_unity_scalar, n);

        // Compute the inverse of n using arkworks
        let n_fr = Fr::from(n as u64);
        let inv_n_fr = n_fr.inverse().unwrap();
        let inv_n_scalar = arkworks_to_fastcrypto(&inv_n_fr);

        for elem in padded_v_hat.iter_mut() {
            *elem = elem.mul(inv_n_scalar);
        }

        v_hat.copy_from_slice(&padded_v_hat[..v_hat.len()]);
    }


    fn element(&self, index: usize) -> Scalar {
        arkworks_to_fastcrypto(&self.domain.element(index))
    }

    fn size(&self) -> usize {
        self.domain.size()
    }
}

/// Convert a BLS12-381 field element from the fastcrypto group to an arkworks field element.
fn fastcrypto_to_arkworks(s: &Scalar) -> Fr {
    let bytes = s.to_byte_array();
    Fr::from_be_bytes_mod_order(&bytes)
}

/// Convert an arkworks field element to a BLS12-381 field element from the fastcrypto group.
fn arkworks_to_fastcrypto(f: &Fr) -> Scalar {
    let bytes: [u8; 32] = f.into_bigint().to_bytes_be().try_into().unwrap();
    Scalar::from_byte_array(&bytes).unwrap()
}

fn fft_in_place_g1_recursive(v: &mut [G1Element], root_of_unity: &Scalar, n: usize) {
    if n <= 1 {
        return;
    }

    let half_n = n / 2;
    let mut even = vec![G1Element::zero(); half_n];
    let mut odd = vec![G1Element::zero(); half_n];

    for i in 0..half_n {
        even[i] = v[i * 2];
        odd[i] = v[i * 2 + 1];
    }

    fft_in_place_g1_recursive(&mut even, root_of_unity, half_n);
    fft_in_place_g1_recursive(&mut odd, root_of_unity, half_n);

    let mut omega = Scalar::from(1); // Initialize omega as 1
    for i in 0..half_n {
        let t = odd[i].mul(omega); // Use exponentiation
        v[i] = even[i] + t;
        v[i + half_n] = even[i] - t;
        omega *= root_of_unity;
    }
}

#[cfg(test)]
mod tests {
    use std::ops::{Add, Mul};

    use fastcrypto::groups::bls12381::{G1Element, G2Element, GTElement, Scalar};
    use fastcrypto::groups::{GroupElement, HashToGroupElement, MultiScalarMul, Pairing};

    use crate::fft::{BLS12381Domain, FFTDomain};

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
        let domain = BLS12381Domain::new(4).unwrap();
        let g = G1Element::generator();
        let v = vec![g.mul(Scalar::from(7)), G1Element::zero(), g, G1Element::zero()];
        let mut v_fft = v.clone();
        domain.fft_in_place_g1(&mut v_fft);
        let mut v_ifft = v_fft.clone();
        domain.ifft_in_place_g1(&mut v_ifft);
        assert_eq!(v, v_ifft);
    }

    #[test]
    fn test_pairing() {
        let a = G1Element::generator().mul(Scalar::from(7));
        let b = G2Element::generator().mul(Scalar::from(5));
        let c = a.pairing(&b);

        let expected = GTElement::generator().mul(Scalar::from(35));
        assert_eq!(c, expected);
    }

    #[test]
    fn test_msm() {
        let g = G1Element::hash_to_group_element(b"hello");
        let h = G1Element::hash_to_group_element(b"world");

        let scalars = vec![Scalar::from(3), Scalar::from(5)];
        let elements = vec![g, h];
        let result = G1Element::multi_scalar_mul(&scalars, &elements).unwrap();

        let expected = g.mul(Scalar::from(3)).add(&h.mul(Scalar::from(5)));
        assert_eq!(result, expected);
    }
}
