use ark_bls12_381::Fr;
use ark_ff::{BigInteger, PrimeField};
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
use fastcrypto::error::{FastCryptoError, FastCryptoResult};
use fastcrypto::groups::bls12381::{Scalar, G1Element};
use fastcrypto::serde_helpers::ToFromByteArray;
use ark_ff::Field;
use fastcrypto::groups::GroupElement;
use std::ops::Mul;
use std::ops::Sub;


use std::ops::Neg;

pub trait FFTDomain: Sized {
    type ScalarType;

    fn new(n: usize) -> FastCryptoResult<Self>;

    fn fft(&self, v: &[Self::ScalarType]) -> Vec<Self::ScalarType>;

    fn ifft(&self, v_hat: &[Self::ScalarType]) -> Vec<Self::ScalarType>;

    fn fft_in_place_g1(&self, v: &mut [G1Element]);

    // fn fft_g1(&self, v: &mut [G1Element]) -> Vec<Self::G1Element>;

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
        let root_of_unity = self.element(1);
        // let root_of_unity_scalar = arkworks_to_fastcrypto(&root_of_unity);
        let result = fft_g1(self, v, &root_of_unity);
        v.copy_from_slice(&result[..v.len()]);
    }

    // fn fft_in_place_g1(&self, v: &mut [G1Element]) {
    //     let n = self.domain.size();
    //     let mut padded_v = v.to_vec();
    //     if padded_v.len() < n {
    //         padded_v.resize(n, G1Element::zero());
    //     } else if padded_v.len() > n {
    //         padded_v.truncate(n);

    //     }

    //     let root_of_unity = self.element(1);
    //     // let root_of_unity_scalar = arkworks_to_fastcrypto(&root_of_unity);
    //     fft_in_place_g1_recursive(self, &mut padded_v, &root_of_unity, n);

    //     v.copy_from_slice(&padded_v[..v.len()]);
    // }

    fn ifft_in_place_g1(&self, v_hat: &mut [G1Element]) {
        let n = self.domain.size();
        let mut padded_v_hat = v_hat.to_vec();
        if padded_v_hat.len() < n {
            padded_v_hat.resize(n, G1Element::zero());
        } else if padded_v_hat.len() > n {
            padded_v_hat.truncate(n);
        }

        let root_of_unity = self.domain.element(n - 1); // Inverse root of unity
        let root_of_unity_scalar = arkworks_to_fastcrypto(&root_of_unity);
        fft_in_place_g1_recursive(self, &mut padded_v_hat, &root_of_unity_scalar, n);

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

fn fastcrypto_to_arkworks(s: &Scalar) -> Fr {
    let bytes = s.to_byte_array();
    Fr::from_be_bytes_mod_order(&bytes)
}

fn arkworks_to_fastcrypto(f: &Fr) -> Scalar {
    let bytes: [u8; 32] = f.into_bigint().to_bytes_be().try_into().unwrap();
    Scalar::from_byte_array(&bytes).unwrap()
}

fn fft_g1(domain: &BLS12381Domain, v: &[G1Element], root_of_unity: &Scalar) -> Vec<G1Element> {
    let n = v.len();
    if n <= 1 {
        return v.to_vec();
    }

    let half_n = n / 2;
    let even: Vec<G1Element> = v.iter().step_by(2).cloned().collect();
    let odd: Vec<G1Element> = v.iter().skip(1).step_by(2).cloned().collect();

    println!("FFT even indices: {:?}", even);
    println!("FFT odd indices: {:?}", odd);

    let even_fft = fft_g1(domain, &even, root_of_unity);
    let odd_fft = fft_g1(domain, &odd, root_of_unity);

    let mut omega = Scalar::from(1);
    let mut result = vec![G1Element::zero(); n];

    for i in 0..half_n {
        println!("omega at step {} for n = {}: {:?}", i, n, omega);
        let t = odd_fft[i].mul(omega);
        result[i] = even_fft[i] + t;
        result[i + half_n] = even_fft[i] - t;
        omega *= root_of_unity;
    }

    println!("FFT result at size {}: {:?}", n, result);

    result
}




fn fft_in_place_g1_recursive(domain: &BLS12381Domain, v: &mut [G1Element], root_of_unity: &Scalar, n: usize) {
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

    fft_in_place_g1_recursive(domain, &mut even, root_of_unity, half_n);
    fft_in_place_g1_recursive(domain, &mut odd, root_of_unity, half_n);

    let mut omega = Scalar::from(1u128); // Initialize omega as 1
    let omega_increment = *root_of_unity;

    for i in 0..half_n {
        let t = odd[i].mul(&omega);
        v[i] = even[i] + t;
        v[i + half_n] = even[i].sub(t);
        omega *= &omega_increment; // Correctly update omega by multiplying with the root of unity increment
    }
}

#[cfg(test)]


mod tests {
    use fastcrypto::groups::bls12381::{G1Element, Scalar};
    use fastcrypto::groups::{GroupElement,Scalar as OtherScalar};
    use std::ops::Mul;
    use rand::thread_rng;

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

        // Create a random vector v
        let v_scalar: Vec<Scalar> = (0..n).map(|_| OtherScalar::rand(&mut rng)).collect();
        let domain = BLS12381Domain::new(n).unwrap();
        let g = G1Element::generator();
        // let v_scalar = vec![Scalar::from(1), Scalar::from(2), Scalar::from(3), Scalar::from(4), Scalar::from(4), Scalar::from(4), Scalar::from(4), Scalar::from(4)];
        //test for domain of size 4
        // let v_scalar = vec![Scalar::from(1), Scalar::from(2), Scalar::from(3), Scalar::from(4)];
        //test 16
        // let v_scalar = vec![Scalar::from(1), Scalar::from(2), Scalar::from(3), Scalar::from(4), Scalar::from(4), Scalar::from(4), Scalar::from(4), Scalar::from(4),Scalar::from(1), Scalar::from(2), Scalar::from(3), Scalar::from(4), Scalar::from(4), Scalar::from(4), Scalar::from(4), Scalar::from(4)];

        let v_group_elements: Vec<G1Element> = v_scalar.iter().map(|x| g.mul(x)).collect();
        let mut v_fft_g = v_group_elements.clone();
        domain.fft_in_place_g1(&mut v_fft_g);
        let v_fft = domain.fft(&v_scalar);
        let expected_v_fft_g: Vec<G1Element> = v_fft.iter().map(|x| g.mul(x)).collect();
        for i in 0..v_fft.len(){
            if expected_v_fft_g[i]!=v_fft_g[i]{
                println!("Value mismatch at index {} ", i);
                println!("expected {:?}", expected_v_fft_g[i]);
                println!("computed {:?}", v_fft_g[i]);
            }
        }
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
