use std::ops::Mul;

use fastcrypto::error::FastCryptoResult;
use fastcrypto::groups::bls12381::{G1Element, G2Element, Scalar};
use fastcrypto::groups::{GroupElement, MultiScalarMul, Pairing, Scalar as OtherScalar};
use itertools::iterate;
use rand::thread_rng;

use crate::fft::{BLS12381Domain, FFTDomain};
use crate::KZG;



#[derive(Clone)]
pub struct KZGTabDFK {
    domain: BLS12381Domain,

    // Jonas: Is this variable necessary?
    a: G1Element,

    g2_tau: G2Element,

    u_vec: Vec<G1Element>,
    l_vec: Vec<G1Element>,
    a_vec: Vec<G1Element>,
    tau_powers_g1: Vec<G1Element>,
}

pub fn multiply_toeplitz_with_vector_fft(
    toeplitz_matrix: &[Scalar],
    vector: &[G1Element]
) -> Vec<G1Element> {
    let m = toeplitz_matrix.len();

    // Instantiate FFT domain with size 2m 
    let domain = BLS12381Domain::new(2 * m).unwrap();

    // Compute FFT of the circulant vector
    // let fft_circulant = domain.fft(&circulant_vector);

    let mut c_vector = vec![Scalar::zero(); 2 * m];
    // println!("{:?}", toeplitz_matrix);
    for i in 0..m{
        c_vector[m+i] = toeplitz_matrix[i];
    }

    // println!("{:?}", c_vector);

    let v = domain.fft(&c_vector);

    // println!("{:?}", v);
    

    // Pad the input vector to match the size of the circulant vector
    let mut padded_vector = vec![G1Element::zero(); 2 * m];
    for i in 0..m {
        padded_vector[i] = vector[i];
    }

    // println!("{:?}", padded_vector);

    // Compute FFT of the padded vector
    domain.fft_in_place_g1(&mut padded_vector);

    let mut omega_vector = vec![Scalar::zero(); 2 * m ];
    for i in 0..2*m{
        let item = domain.element(i);
        omega_vector[i] = item;
    }


    // println!("{:?}", omega_vector);

    // Element-wise multiplication in the frequency domain
    let fft_result_1: Vec<G1Element> = v.iter()
        .zip(padded_vector.iter())
        .map(|(a, b)| b.mul(*a)) 
        .collect();

    // Second product in the exponent
    let fft_result: Vec<G1Element> = fft_result_1.iter()
        .zip(omega_vector.iter())
        .map(|(a, b)| a.mul(*b)) 
        .collect();
    // fft_result
    // Compute inverse FFT
    let mut ifft_result = fft_result.clone();
    domain.ifft_in_place_g1(&mut ifft_result);
    // ifft_result


    // Extract the result corresponding to the original vector length
    let res = ifft_result[..m].to_vec();
    let mut result = res.clone();
    domain.fft_in_place_g1(&mut result);
    result

}

impl KZG for KZGTabDFK {
    // Uses the BLS12-381 construction
    type G = G1Element;

    fn new(n: usize) -> FastCryptoResult<Self> {
        let domain = BLS12381Domain::new(n)?;

        // Generate tau using a random scalar
        let tau = Scalar::rand(&mut thread_rng());

        let g2_tau = G2Element::generator() * tau;

        //Compute a = g^{A(tau)} where A = X^n-1
        let g_tau_n = (0..n).fold(G1Element::generator(), |acc, _| acc * tau);

        //let g_tau_n = tau_powers_g1[n - 1].mul(tau); // g^{tau^n}
        let a = g_tau_n - G1Element::generator();

        // Compute a_i = g^{A(tau)/tau - omega^i}
        let mut a_vec = vec![G1Element::zero(); n];

        //Compute u_i = g^{(L_i(tau) - 1)/(tau-omega^i)}
        let mut u_vec = vec![G1Element::zero(); n];

        // Compute tau^i for i = 0 to n-1
        let tau_powers_g: Vec<Scalar> =
            iterate(Scalar::generator(), |g| g * tau).take(n).collect();
        let tau_powers_g1: Vec<G1Element> = itertools::iterate(G1Element::generator(), |g| g.mul(tau))
            .take(n)
            .collect();

        // Compute l_i = g^L_i(tau)
        let l_vec: Vec<G1Element> = domain
            .ifft(&tau_powers_g)
            .iter()
            .map(|s| G1Element::generator() * s)
            .collect();

        let mut omega_i = domain.element(0);
        for i in 0..n {
            // Compute a_i
            let denom = tau - omega_i;
            let a_i = a.mul(denom.inverse().unwrap());
            a_vec[i] = a_i;

            // Compute u_i
            let l_i_minus_1 = l_vec[i] - G1Element::generator();
            let denom = tau - omega_i;
            let u_i = l_i_minus_1.mul(denom.inverse().unwrap());
            u_vec[i] = u_i;

            // Update omega_i
            if i < n - 1 {
                omega_i *= domain.element(1);
            }
        }

        Ok(Self {
            domain,
            g2_tau,
            a,
            u_vec,
            l_vec,
            a_vec,
            tau_powers_g1,
        })
    }

    fn commit(&self, v: &[Scalar]) -> G1Element {
        G1Element::multi_scalar_mul(&v, &self.l_vec).unwrap()
    }

    fn open(&self, v: &[Scalar], index: usize) -> G1Element {
        let mut open = G1Element::zero();
        for j in 0..v.len() {
            if j != index {
                let omega_i = self.domain.element(index);
                let omega_j = self.domain.element(j);

                // Compute c_i and c_j
                let c_i = (omega_i - omega_j).inverse().unwrap();
                let c_j = (omega_j - omega_i).inverse().unwrap();

                // Compute w_ij = a_i^c_i * a_j^c_j
                let w_ij = self.a_vec[index].mul(c_i) + self.a_vec[j].mul(c_j);

                // Compute u_ij = w_ij^{omega_j / n}
                let omega_j_n = (omega_j / Scalar::from(v.len() as u128)).unwrap();
                let u_ij = w_ij.mul(omega_j_n);

                open += u_ij.mul(v[j]);
            }
        }
        open += self.u_vec[index].mul(v[index]);
        open
    }

    fn open_all(&self, v: &[Scalar], indices: Vec<usize>) -> Vec<G1Element> {
        let poly = self.domain.ifft(v);
        let mut t = self.tau_powers_g1.clone();
        t.reverse();
        let h = multiply_toeplitz_with_vector_fft(&poly, &t);
        h
    }

    fn verify(
        &self,
        index: usize,
        v_i: &Scalar,
        commitment: &G1Element,
        open_i: &G1Element,
    ) -> bool {
        let lhs = *commitment - G1Element::generator() * v_i;

        let rhs = self.g2_tau - G2Element::generator() * self.domain.element(index);

        // Perform the pairing check e(lhs, g) == e(open_i, rhs)
        let lhs_pairing = lhs.pairing(&G2Element::generator());
        let rhs_pairing = open_i.pairing(&rhs);

        lhs_pairing == rhs_pairing
    }

    fn update(
        &self,
        commitment: &mut G1Element,
        index: usize,
        old_v_i: &Scalar,
        new_v_i: &Scalar,
    ) -> G1Element {
        *commitment + self.l_vec[index].mul(new_v_i - old_v_i)
    }

    fn update_open_i(
        &self,
        open: &mut G1Element,
        index: usize,
        old_v_i: &Scalar,
        new_v_i: &Scalar,
    ) -> G1Element {
        *open + self.u_vec[index].mul(new_v_i - old_v_i)
    }

    fn update_open_j(
        &self,
        open: &mut G1Element,
        index: usize,
        index_j: usize,
        old_v_j: &Scalar,
        new_v_j: &Scalar,
    ) -> G1Element {
        let omega_i = self.domain.element(index);
        let omega_j = self.domain.element(index_j);

        // Compute c_i and c_j
        let c_i = (omega_i - omega_j).inverse().unwrap();
        let c_j = (omega_j - omega_i).inverse().unwrap();

        // Compute w_ij = a_i^c_i * a_j^c_j
        let w_ij = self.a_vec[index].mul(c_i) + self.a_vec[index_j].mul(c_j);

        // Compute u_ij = w_ij^{omega_j / n}
        let omega_j_n = (omega_j / Scalar::from(self.domain.size() as u128)).unwrap();
        let u_ij = w_ij.mul(omega_j_n);
        *open + u_ij.mul(new_v_j - old_v_j)
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

        // Create a new KZGTabDFK struct
        let n = 8;
        let kzg = KZGTabDFK::new(n).unwrap();

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
    fn test_kzg_commit_open_update_i_verify() {
        let mut rng = rand::thread_rng();

        // Create a new KZGTabDFK struct
        let n = 8;
        let kzg = KZGTabDFK::new(n).unwrap();

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
        let new_commitment = kzg.update(&mut commitment, index, &v[index], &new_v_index);

        //Update the opening
        let new_opening = kzg.update_open_i(&mut open_value, index, &v[index], &new_v_index);

        //Verify the updated opening
        let is_valid = kzg.verify(index, &new_v_index, &new_commitment, &new_opening);

        // Assert that the verification passes
        assert!(
            is_valid,
            "Verification of the opening after updating should succeed."
        );
    }

    #[test]
    fn test_kzg_commit_open_update_j_verify() {
        let mut rng = rand::thread_rng();

        // Create a new KZGTabDFK struct
        let n = 8;
        let kzg = KZGTabDFK::new(n).unwrap();

        // Create a random vector v
        let v: Vec<Scalar> = (0..n).map(|_| OtherScalar::rand(&mut rng)).collect();

        println!("{:?}", v);

        // Create a commitment
        let mut commitment = kzg.commit(&v);

        // Pick a random index to open
        let index = rng.gen_range(0..n);

        // Create an opening
        let mut open_value = kzg.open(&v, index);

        // Pick a new index to update¨
        let mut index_j;
        loop {
            index_j = rng.gen_range(0..n);
            if index_j != index {
                break;
            }
        }

        // Set a new value for v_i
        let new_v_index_j = Scalar::rand(&mut rng);

        //Update the commitment
        let new_commitment = kzg.update(&mut commitment, index_j, &v[index_j], &new_v_index_j);

        //Update the opening
        let new_opening =
            kzg.update_open_j(&mut open_value, index, index_j, &v[index_j], &new_v_index_j);

        //Verify the updated opening
        let is_valid = kzg.verify(index, &v[index], &new_commitment, &new_opening);

        // Assert that the verification passes
        assert!(
            is_valid,
            "Verification of the opening after updating j's value should succeed."
        );
    }
}
