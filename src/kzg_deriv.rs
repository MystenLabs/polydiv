use fastcrypto::error::{FastCryptoError, FastCryptoResult};
use fastcrypto::groups::bls12381::{G1Element, G2Element, Scalar};
use fastcrypto::groups::{GroupElement, MultiScalarMul, Pairing, Scalar as OtherScalar};
use rand::thread_rng;
use std::ops::Mul;

use crate::fft::{BLS12381Domain, FFTDomain};
use crate::KZG;


pub struct KZGDeriv {
    domain: BLS12381Domain,
    tau_powers_g1: Vec<G1Element>,
    tau_powers_g2: Vec<G2Element>,
    w_vec: Vec<G1Element>,
    u_vec: Vec<G1Element>,
}

impl KZGDeriv {
    pub fn new(n: usize) -> FastCryptoResult<Self> {
        let domain = BLS12381Domain::new(n)?;

        // Generate tau using a random scalar
        let tau = Scalar::rand(&mut thread_rng());

        // Compute g^tau^i for i = 0 to n-1 in G1
        let tau_powers_g1: Vec<G1Element> = itertools::iterate(G1Element::generator(), |g| g * tau)
            .take(n)
            .collect();

        // Compute g^tau^i for i = 0 to n-1 in G2
        let tau_powers_g2: Vec<G2Element> = itertools::iterate(G2Element::generator(), |g| g * tau)
            .take(n)
            .collect();
        
        // Compute w_vec 
        let mut w_vec = vec![G1Element::zero(); n];
        for i in 0..n {
            let mut num = Scalar::from(1u128);
            let mut denom = Scalar::from(1u128);
            let omega_i = domain.element(i);

            for j in 0..n {
                if i != j {
                    let omega_j = domain.element(j);
                    num *= tau - omega_j;
                    denom *= omega_i - omega_j;
                }
            }

            let w_i = num * denom.inverse().unwrap();
            w_vec[i] = G1Element::generator().mul(w_i);
        }

         //Compute u_i = g^{(L_i(tau) - 1)/(tau-omega^i)}
         let mut u_vec = vec![G1Element::zero(); n];
         for i in 0..n {
             let omega_i = domain.element(i);
             let l_i = w_vec[i];
             let l_i_minus_1 = l_i - tau_powers_g1[0];
             let denom = tau - omega_i;
             let u_i = l_i_minus_1.mul(denom.inverse().unwrap());
             u_vec[i] = u_i;
         }
        

        Ok(Self {
            domain,
            tau_powers_g1,
            tau_powers_g2,
            w_vec,
            u_vec,
        })
    }
}

impl KZG for KZGDeriv {
    // Uses the BLS12-381 construction
    type G = G1Element;

    fn commit(&self, v: &[Scalar]) -> G1Element {
        G1Element::multi_scalar_mul(&v, &self.w_vec).unwrap()
    }

    fn open(&self, v: &[Scalar], index: usize) -> G1Element {
        // let mut poly = self.domain.ifft(&v);
        // poly[0] -= &v[index];

        // let divisor = [-self.domain.element(index), Scalar::generator()];
        // let (quotient, _) = polynomial_division(&poly, &divisor).unwrap();

        let mut e_vec = vec![Scalar::zero(); v.len()];
        for j in (0..v.len()){
            if j!=index{
                let omega_i = self.domain.element(index);
                let omega_j = self.domain.element(j);
                let omega_j_minus_i = (omega_j/omega_i).unwrap();
                e_vec[j] = (omega_j_minus_i/(omega_i - omega_j)).unwrap();
            }
            let omega_i = self.domain.element(index);
            e_vec[index] = (Scalar::from((v.len() -1) as u128)/(Scalar::from(2u128)*omega_i)).unwrap();
        }

        let mut to_mul = vec![Scalar::zero(); v.len()];
        // to_mul[index] = e_vec[index]*v[index];
        // for j in (0..v.len()){
        //     if j!=index{
        //         let omega_i = self.domain.element(index);
        //         let omega_j = self.domain.element(j);
        //         let to_add = (v[j]-v[index])/(omega_j - omega_i);
        //         to_mul[j] = e_vec[j]*v[j]+to_add.unwrap();
        //     }
            
        // }
        let mut v_prime = Scalar::zero();
    
        for (v_i, e_i) in v.iter().zip(e_vec.iter()) {
        v_prime += v_i*e_i;
        }
    
        
        for j in (0..v.len()){
            if j!=index{
                let omega_i = self.domain.element(index);
                let omega_j = self.domain.element(j);
                let to_mul_j = (v[j]-v[index])/(omega_j - omega_i);
                to_mul[j] = to_mul_j.unwrap();
            }
        }
        to_mul[index] = v_prime;

        G1Element::multi_scalar_mul(&to_mul, &self.w_vec[..to_mul.len()]).unwrap()
    }

    fn verify(
        &self,
        index: usize,
        v_i: &Scalar,
        commitment: &G1Element,
        open_i: &G1Element,
    ) -> bool {
        let lhs = *commitment - self.tau_powers_g1[0] * v_i;

        let rhs = self.tau_powers_g2[1] - self.tau_powers_g2[0] * self.domain.element(index);

        // Perform the pairing check e(lhs, g) == e(open_i, rhs)
        let lhs_pairing = lhs.pairing(&self.tau_powers_g2[0]);
        let rhs_pairing = open_i.pairing(&rhs);

        lhs_pairing == rhs_pairing
    }

    fn update(&self, commitment: &mut G1Element, index: usize,old_v_i: &Scalar,  new_v_i: &Scalar) -> G1Element {
        *commitment + self.w_vec[index].mul(new_v_i - old_v_i)
    }

    fn update_open_i(&self, open: &mut G1Element, index: usize, old_v_i: &Scalar, new_v_i: &Scalar) -> G1Element{
        *open + self.u_vec[index].mul(new_v_i - old_v_i)
    }

    fn update_open_j(&self, open: &mut G1Element, index: usize, index_j:usize, old_v_j: &Scalar,  new_v_j: &Scalar) -> G1Element{
        let omega_i = self.domain.element(index);
        let omega_j = self.domain.element(index_j);
        let to_mul_1 = (new_v_j - old_v_j)/(omega_j - omega_i);
        let omega_j_minus_i = (omega_j/omega_i).unwrap();
        let to_mul = (omega_j_minus_i)/(omega_i - omega_j);
        let to_mul_2 = (new_v_j - old_v_j)*to_mul.unwrap();

        *open + self.w_vec[index_j].mul(to_mul_1.unwrap()) + self.w_vec[index].mul(to_mul_2)
    }
}

#[cfg(test)]
mod tests {
    use fastcrypto::groups::bls12381::Scalar;
    use rand::Rng;

    use super::*;

    #[test]
    fn test_kzg_commit_open_verify() {
        let mut rng = rand::thread_rng();

        // Create a new KZGDeriv struct
        let n = 8;
        let kzg = KZGDeriv::new(n).unwrap();

        // Create a random vector v
        let v: Vec<Scalar> = (0..n).map(|_| OtherScalar::rand(&mut rng)).collect();

        println!("{:?}", v);

        // Create a commitment
        let commitment = kzg.commit(&v);

        // Pick a random index to open
        let index = rng.gen_range(0..n);

        // Create an opening
        let open_value = kzg.open(&v, index);

        // Verify the opening
        let is_valid = kzg.verify(index, &v[index], &commitment, &open_value);

        // Assert that the verification passes
        assert!(is_valid, "Verification of the opening should succeed.");
    }

    #[test]
    fn test_kzg_commit_open_update_i_verify(){
        let mut rng = rand::thread_rng();

        // Create a new KZGDeriv struct
        let n = 8;
        let kzg = KZGDeriv::new(n).unwrap();

        // Create a random vector v
        let v: Vec<Scalar> = (0..n).map(|_| OtherScalar::rand(&mut rng)).collect();

        println!("{:?}", v);

        // Create a commitment
        let mut commitment = kzg.commit(&v);

        // Pick a random index to open
        let index = rng.gen_range(0..n);

        // Create an opening
        let mut open_value = kzg.open(&v, index);

        // Set a new valie for v_i
        let new_v_index = Scalar::rand(&mut rng);

        //Update the commitment
        let mut new_commitment = kzg.update(&mut commitment, index, &v[index], &new_v_index);

        //Update the opening
        let mut new_opening = kzg.update_open_i(&mut open_value, index, &v[index],&new_v_index);


        //Verify the updated opening
        let is_valid = kzg.verify(index, &new_v_index, &new_commitment, &new_opening);

        // Assert that the verification passes
        assert!(is_valid, "Verification of the opening after updating should succeed.");
    }

    #[test]
    fn test_kzg_commit_open_update_j_verify(){
        let mut rng = rand::thread_rng();

        // Create a new KZGDeriv struct
        let n = 8;
        let kzg = KZGDeriv::new(n).unwrap();

        // Create a random vector v
        let v: Vec<Scalar> = (0..n).map(|_| OtherScalar::rand(&mut rng)).collect();

        println!("{:?}", v);

        // Create a commitment
        let mut commitment = kzg.commit(&v);

        // Pick a random index to open
        let index = rng.gen_range(0..n);

        // Create an opening
        let mut open_value = kzg.open(&v, index);

        // Pick a new index to update

        let index_j = rng.gen_range(0..n);


        // Set a new value for v_i
        let new_v_index_j = Scalar::rand(&mut rng);

        //Update the commitment
        let mut new_commitment = kzg.update(&mut commitment, index_j, &v[index_j], &new_v_index_j);

        //Update the opening
        let mut new_opening = kzg.update_open_j(&mut open_value, index, index_j, &v[index_j],&new_v_index_j);


        //Verify the updated opening
        let is_valid = kzg.verify(index, &v[index], &new_commitment, &new_opening);

        // Assert that the verification passes
        assert!(is_valid, "Verification of the opening after updating j's value should succeed.");
    }
}
