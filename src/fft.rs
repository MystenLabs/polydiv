use ark_bls12_381::Fr;
use ark_ff::{BigInteger, PrimeField};
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
use fastcrypto::error::{FastCryptoError, FastCryptoResult};
use fastcrypto::groups::bls12381::Scalar;
use fastcrypto::serde_helpers::ToFromByteArray;

pub trait FFTDomain<G>: Sized {
    /// Create a new domain with at least `n` elements. Fails
    fn new(n: usize) -> FastCryptoResult<Self>;

    /// Compute the FFT of a vector of field elements. If the input is smaller than
    /// the domain size, it is padded with zeros. If it is larger, the input is truncated.
    fn fft(&self, v: &[G]) -> Vec<G>;

    /// Compute the inverse FFT of a vector of field elements. If the input is smaller than
    /// the domain size, it is padded with zeros. If it is larger, the input is truncated.
    fn ifft(&self, v_hat: &[G]) -> Vec<G>;

    /// Get the nth root of unity of the domain.
    fn root_of_unity(&self) -> G;

    /// Get the size of the domain.
    fn size(&self) -> usize;
}

pub struct BLS12381Domain {
    domain: GeneralEvaluationDomain<Fr>,
}

impl FFTDomain<Scalar> for BLS12381Domain {
    fn new(n: usize) -> FastCryptoResult<Self> {
        let domain = GeneralEvaluationDomain::<Fr>::new(n).ok_or(FastCryptoError::InvalidInput)?;
        Ok(Self { domain })
    }

    fn fft(&self, v: &[Scalar]) -> Vec<Scalar> {
        let scalars = v.iter().map(fastcrypto_to_arkworks).collect::<Vec<_>>();
        self.domain
            .fft(&scalars)
            .iter()
            .map(arkworks_to_fastcrypto)
            .collect::<Vec<_>>()
    }

    fn ifft(&self, v_hat: &[Scalar]) -> Vec<Scalar> {
        let scalars = v_hat.iter().map(fastcrypto_to_arkworks).collect::<Vec<_>>();
        self.domain
            .ifft(&scalars)
            .iter()
            .map(arkworks_to_fastcrypto)
            .collect::<Vec<_>>()
    }

    fn root_of_unity(&self) -> Scalar {
        arkworks_to_fastcrypto(&self.domain.group_gen())
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

#[cfg(test)]
mod tests {
    use crate::fft::{BLS12381Domain, FFTDomain};
    use fastcrypto::groups::bls12381::{G1Element, G2Element, GTElement, Scalar};
    use fastcrypto::groups::{GroupElement, HashToGroupElement, MultiScalarMul, Pairing};
    use std::ops::{Add, Mul};

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
