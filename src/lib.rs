// Declare the modules

use fastcrypto::error::FastCryptoResult;
use fastcrypto::groups::GroupElement;

pub mod kzg_deriv;
pub mod kzg_fk;
pub mod kzg_original;
pub mod kzg_tabdfk;

pub mod fft;

pub trait KZG: Sized + Clone {
    type G: GroupElement;

    fn new(n: usize) -> FastCryptoResult<Self>;

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
        old_v_i: &<Self::G as GroupElement>::ScalarType,
        new_v_i: &<Self::G as GroupElement>::ScalarType,
    ) -> Self::G;

    fn update_open_i(
        &self,
        open: &mut Self::G,
        index: usize,
        old_v_i: &<Self::G as GroupElement>::ScalarType,
        new_v_i: &<Self::G as GroupElement>::ScalarType,
    ) -> Self::G;

    fn update_open_j(
        &self,
        open: &mut Self::G,
        index: usize,
        index_j: usize,
        old_v_j: &<Self::G as GroupElement>::ScalarType,
        new_v_j: &<Self::G as GroupElement>::ScalarType,
    ) -> Self::G;
}
