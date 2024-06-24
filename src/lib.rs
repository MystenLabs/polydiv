use ark_ff::{BigInteger, PrimeField};
use ark_poly::EvaluationDomain;
use fastcrypto::serde_helpers::ToFromByteArray;

mod fft;

pub trait KZG<G, C> {
    fn commit(&self, v: &[G]) -> C;
    fn open(&self, commitment: &C, index: usize) -> G;
    fn verify(&self, index: usize, v_i: &G, commitment: &C, open_i: &G) -> bool;
    fn update(&self, commitment: &mut C, index: usize, new_v_i: &G) -> C;
}
