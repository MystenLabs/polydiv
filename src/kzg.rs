// Declare the modules
use ark_ff::{BigInteger, PrimeField};
use ark_poly::EvaluationDomain;
use fastcrypto::serde_helpers::ToFromByteArray;
use fastcrypto::groups::bls12381::Scalar;
use ark_bls12_381::{Fr, G1Projective, G1Affine, G2Affine, Bls12_381};


pub trait KZG<G, C> {
    fn commit(&self, v: &[G]) -> C;
    fn open(&self, v: &[Scalar], index: usize) -> G1Projective;
    fn verify(&self, index: usize, v_i: &G, commitment: &C, open_i: &G1Projective) -> bool;
    fn update(&self, commitment: &mut C, index: usize, new_v_i: &G) -> C;
}