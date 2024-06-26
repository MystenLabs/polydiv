// use crate::fft::{FFTDomain, BLS12381Domain};
// use fastcrypto::groups::bls12381::Scalar;

// pub struct KZGFK {
//     domain: BLS12381Domain,
// }

// impl KZGFK {
//     pub fn new(n: usize) -> Self {
//         Self { domain: BLS12381Domain::new(n).unwrap() }
//     }
// }

// impl crate::KZG<Scalar, Vec<Scalar>> for KZGFK {
//     fn commit(&self, v: &[Scalar]) -> Vec<Scalar> {
//         // FK implementation here
//     }

//     fn open(&self, commitment: &Vec<Scalar>, index: usize) -> Scalar {
//         // FK implementation here
//     }

//     fn verify(&self, index: usize, v_i: &Scalar, commitment: &Vec<Scalar>, open_i: &Scalar) -> bool {
//         // FK implementation here
//     }

//     fn update(&self, commitment: &mut Vec<Scalar>, index: usize, new_v_i: &Scalar) -> Vec<Scalar> {
//         // FK implementation here
//     }
// }
