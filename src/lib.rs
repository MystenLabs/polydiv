// Declare the modules

use ark_ff::{BigInteger, PrimeField};
use ark_poly::EvaluationDomain;
use fastcrypto::serde_helpers::ToFromByteArray;

pub mod kzg_original;
pub mod kzg_fk;
pub mod kzg_tabdfk;

mod fft;

pub trait KZG<G, C> {
    fn commit(&self, v: &[G]) -> C;
    fn open(&self, v: &[G], index: usize) -> C;
    fn verify(&self, index: usize, v_i: &G, commitment: &C, open_i: &C) -> bool;
    fn update(&self, commitment: &mut C, index: usize, new_v_i: &G) -> C;
}
