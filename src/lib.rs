// Declare the modules

use fastcrypto::groups::GroupElement;

pub mod kzg_fk;
pub mod kzg_original;
pub mod kzg_tabdfk;

pub mod fft;

pub trait KZG {
    type G: GroupElement;

    fn commit(&self, v: &[<Self::G as GroupElement>::ScalarType]) -> Self::G;
    fn open(&self, v: &[<Self::G as GroupElement>::ScalarType], index: usize) -> Self::G;
    fn verify(
        &self,
        index: usize,
        v_i: &<Self::G as GroupElement>::ScalarType,
        commitment: &Self::G,
        open_i: &Self::G,
    ) -> bool;
    fn update(
        &self,
        commitment: &mut Self::G,
        index: usize,
        new_v_i: &<Self::G as GroupElement>::ScalarType,
    ) -> Self::G;
}
