use std::ops::Mul;

use ark_bls12_381::Fr;
use ark_ff::{BigInteger, PrimeField};
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
use fastcrypto::error::{FastCryptoError, FastCryptoResult};
use fastcrypto::groups::bls12381::Scalar;
use fastcrypto::groups::{GroupElement, Scalar as ScalarTrait};
use fastcrypto::serde_helpers::ToFromByteArray;

pub trait FFTDomain: Sized {
    type ScalarType: ScalarTrait;

    /// Create a new domain with n elements.
    fn new(n: usize) -> FastCryptoResult<Self>;

    /// Compute the FFT of a vector of scalars.
    fn fft(&self, v: &[Self::ScalarType]) -> Vec<Self::ScalarType> {
        fft_group(v, &self.element(1))
    }

    /// Compute the IFFT of a vector of scalars.
    fn ifft(&self, v_hat: &[Self::ScalarType]) -> Vec<Self::ScalarType> {
        fft_group(v_hat, &self.element(self.size() - 1))
            .iter()
            .map(|x| x.mul(&self.size_inv()))
            .collect()
    }

    /// Perform a FFT on elements of a group using the domains scalar type.
    fn fft_in_place_group<G: GroupElement<ScalarType = Self::ScalarType>>(&self, v: &mut [G]) {
        // TODO: Implement actual in-place algorithm.
        // TODO: Do we need padding?
        let result = fft_group(&v, &self.element(1));
        v.copy_from_slice(&result);
    }

    /// Perform an IFFT on elements of a group using the domains scalar type.
    fn ifft_in_place_group<G: GroupElement<ScalarType = Self::ScalarType>>(&self, v_hat: &mut [G]) {
        let mut result = fft_group(&v_hat, &self.element(self.size() - 1));
        let n_inverse = self.size_inv();
        for elem in result.iter_mut() {
            *elem = elem.mul(&n_inverse);
        }
        v_hat.copy_from_slice(&result);
    }

    /// Get the i-th element of the domain which is the n-th root of unity to the index-th power.
    /// This may be computed on request so it is not suitable for repeated calls.
    fn element(&self, index: usize) -> Self::ScalarType;

    /// Get the size of the domain.
    fn size(&self) -> usize;

    /// Get the inverse of the size of the domain as a scalar.
    fn size_inv(&self) -> Self::ScalarType;
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
        let result = self
            .domain
            .fft(&scalars)
            .iter()
            .map(arkworks_to_fastcrypto)
            .collect::<Vec<_>>();
        result
    }

    fn ifft(&self, v_hat: &[Scalar]) -> Vec<Scalar> {
        let scalars = v_hat.iter().map(fastcrypto_to_arkworks).collect::<Vec<_>>();
        let result = self
            .domain
            .ifft(&scalars)
            .iter()
            .map(arkworks_to_fastcrypto)
            .collect::<Vec<_>>();
        result
    }

    fn element(&self, index: usize) -> Scalar {
        arkworks_to_fastcrypto(&self.domain.element(index))
    }

    fn size(&self) -> usize {
        self.domain.size()
    }

    fn size_inv(&self) -> Self::ScalarType {
        arkworks_to_fastcrypto(&self.domain.size_inv())
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

/// Helper function for FFT function for group elements.
fn fft_group<G: GroupElement>(v: &[G], root_of_unity: &<G as GroupElement>::ScalarType) -> Vec<G> {
    fft_group_with_offset(v, 0, 1, v.len(), root_of_unity)
}

/// Compute the FFT of the elements {v_i} for i = offset, offset + step, offset + 2*step, ..., offset + (n-1)*step
fn fft_group_with_offset<G: GroupElement>(
    v: &[G],
    offset: usize,
    step: usize,
    n: usize,
    root_of_unity: &<G as GroupElement>::ScalarType,
) -> Vec<G> {
    if n <= 1 {
        return vec![v[offset]];
    }

    // TODO: This implicitly assumes that n is even
    let half_n = n / 2;
    let root_of_unity_squared = root_of_unity.mul(root_of_unity);
    let even_fft = fft_group_with_offset(v, offset, 2 * step, n - half_n, &root_of_unity_squared);
    let odd_fft = fft_group_with_offset(v, offset + step, 2 * step, half_n, &root_of_unity_squared);

    let mut omega = G::ScalarType::from(1);
    let mut result = vec![G::zero(); n];

    for i in 0..half_n {
        let t = odd_fft[i].mul(&omega);
        result[i] = even_fft[i] + t;
        result[i + half_n] = even_fft[i] - t;
        if i < half_n - 1 {
            omega = root_of_unity.mul(omega);
        }
    }

    result
}

#[cfg(test)]
mod tests {
    use std::ops::Mul;

    use fastcrypto::groups::bls12381::{G1Element, Scalar};
    use fastcrypto::groups::{GroupElement, Scalar as OtherScalar};

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
        let mut rng = rand::thread_rng();
        let n = 8;
        let v_scalar: Vec<Scalar> = (0..n).map(|_| OtherScalar::rand(&mut rng)).collect();
        let domain = BLS12381Domain::new(n).unwrap();
        let g = G1Element::generator();
        let v_group_elements: Vec<G1Element> = v_scalar.iter().map(|x| g.mul(x)).collect();
        let mut v_fft_g = v_group_elements.clone();
        domain.fft_in_place_group(&mut v_fft_g);
        let v_fft = domain.fft(&v_scalar);
        let expected_v_fft_g: Vec<G1Element> = v_fft.iter().map(|x| g.mul(x)).collect();
        assert_eq!(v_fft_g, expected_v_fft_g);
    }

    #[test]
    fn test_ifft_g1() {
        let domain = BLS12381Domain::new(8).unwrap();
        let g = G1Element::generator();
        let v_scalar = vec![
            Scalar::from(1),
            Scalar::from(2),
            Scalar::from(3),
            Scalar::from(4),
            Scalar::from(4),
            Scalar::from(4),
            Scalar::from(4),
            Scalar::from(4),
        ];
        let v_group_elements: Vec<G1Element> = v_scalar.iter().map(|x| g.mul(x)).collect();
        let mut v_ifft_g = v_group_elements.clone();
        domain.ifft_in_place_group(&mut v_ifft_g);
        let v_ifft = domain.ifft(&v_scalar);
        let expected_v_ifft_g: Vec<G1Element> = v_ifft.iter().map(|x| g.mul(x)).collect();
        assert_eq!(v_ifft_g, expected_v_ifft_g);
    }
}
