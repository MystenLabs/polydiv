use fastcrypto::error::{FastCryptoError, FastCryptoResult};
use fastcrypto::groups::bls12381::{G1Element, G2Element, Scalar};
use fastcrypto::serde_helpers::ToFromByteArray;
use fastcrypto::groups::{GroupElement, MultiScalarMul, Pairing, Scalar as OtherScalar};
use rand::thread_rng;
use ark_ff::FftField;
use ark_poly::{Polynomial, EvaluationDomain};
use ark_poly::polynomial::univariate::DensePolynomial as Poly;
use std::ops::Mul;


use crate::fft::{BLS12381Domain, FFTDomain};
use crate::KZG;



// testing purposes only
fn compute_matrix_vector_multiplication(
    coefficients: &[Scalar],
    tau_powers: &[G1Element],
) -> Vec<G1Element> {
    let d = coefficients.len();
    let mut result = vec![G1Element::zero(); d];
    let mut matrix = vec![vec![Scalar::zero(); d]; d];

    // Iterate over each row of the result vector
    let mut coeff = coefficients.clone().to_vec();
    coeff.reverse();

    for i in 0..d{
        // println!("Row number {:?}",i);
        // println!("{:?}",coeff);
        matrix[i] = coeff.clone();
        coeff.pop();
        coeff.insert(0,Scalar::zero());
    }

    let rows = matrix.len();
    let cols = matrix[0].len();

    for i in 0..rows {
        for j in 0..cols {
            result[i] += tau_powers[j].mul(matrix[i][j]);
        }
    }

    println!("{:?}", result.len());
    result
}

pub fn multiply_poly_with_vector_fft(
    poly_coeff: &[Scalar],
    vector: &[G1Element]
) -> Vec<G1Element> {
    let m = poly_coeff.len() - 1;//degree = d
    println!("{:?}", m);


    // Instantiate FFT domain with size 2m 
    let domain = BLS12381Domain::new(2 * m).unwrap(); //2^k closest to 2m

    // Compute FFT of the circulant vector
    // let fft_circulant = domain.fft(&circulant_vector);

    let mut c_vector = vec![Scalar::zero(); 2 * m];
    // println!("{:?}", poly_coeff);
    for i in 0..m{
        c_vector[m+i] = poly_coeff[i+1];//c[m] = f_1, c[m+1] = f_2...c[2m] = f[m+1] = f_d
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
        let item = domain.element(i);//these are the n-th roots of unity when n = 2^k, n is the closest power of 2 to 2m. Do I need 2m-th root of unity instead? 
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
    ifft_result


    

}





#[derive(Clone)]
pub struct KZGFK {
    domain: BLS12381Domain,
    tau_powers_g1: Vec<G1Element>,
    g2_tau: G2Element,
}

impl KZG for KZGFK {
    // Uses the BLS12-381 construction
    type G = G1Element;

    fn new(n: usize) -> FastCryptoResult<Self> {
        let domain = BLS12381Domain::new(n)?;

        // Generate tau using a random scalar
        let tau = fastcrypto::groups::bls12381::Scalar::rand(&mut thread_rng());

        // Compute g^tau^i for i = 0 to n-1 in G1
        let tau_powers_g1: Vec<G1Element> = itertools::iterate(G1Element::generator(), |g| g.mul(tau))
            .take(n)
            .collect();

        // println!("{:?}",tau_powers_g1);

        // Compute g^tau^i for i = 0 to n-1 in G2
        let g2_tau = G2Element::generator().mul(tau);

        Ok(Self {
            domain,
            tau_powers_g1,
            g2_tau,
        })
    }

    fn commit(&self, v: &[Scalar]) -> G1Element {
        let poly = self.domain.ifft(v);
        G1Element::multi_scalar_mul(&poly, &self.tau_powers_g1).unwrap()
        // let mut result = G1Element::zero();
        // for i in 0..v.len(){
        //     result+=self.tau_powers_g1[i].mul(v[i]);
        // }
        // result
    }

    fn open(&self, v: &[Scalar], index: usize) -> G1Element {
        let poly = self.domain.ifft(v);
    
        let mut quotient_coeffs: Vec<Scalar> = vec![Scalar::zero(); poly.len() - 1];

        quotient_coeffs[poly.len() - 2] = poly[poly.len() - 1];

        let omega_i = self.domain.element(index);
        for j in (0..poly.len() - 2).rev() {
            quotient_coeffs[j] = poly[j + 1] + quotient_coeffs[j + 1] * omega_i;
        }

        G1Element::multi_scalar_mul(
            &quotient_coeffs,
            &self.tau_powers_g1[..quotient_coeffs.len()],
        )
        .unwrap()
    }

    fn open_all(&self, v: &[Scalar], indices: Vec<usize>) -> Vec<G1Element> {
        let mut poly = self.domain.ifft(v); //f_0 to f_7 
        // println!("{:?}", poly);
        // let mut omega_vector = vec![Scalar::zero(); v.len()];
        // for i in 0..v.len(){
        //     let item = self.domain.element(i);//these are the n-th roots of unity when n = 2^k, n is the closest power of 2 to 2m. Do I need 2m-th root of unity instead? 
        //     omega_vector[i] = item;
        // }
        // let mut poly = lagrange_interpolation_coefficients(&omega_vector, &v);
        println!("The very OG {:?}", poly);
        let degree = poly.len() - 1; 
        // println!("{:?}", self.tau_powers_g1);
        let mut t = self.tau_powers_g1.clone();
        t.truncate(degree);
        t.reverse(); //[t^d-1], ... [1]

        // println!("{:?}", t.len());
        // // poly.reverse();
        

        // // println!("{:?}", poly);
        // let mut h = multiply_poly_with_vector_fft(&poly, &t);
        // h.truncate(degree);
        // println!("{:?}", h);

        // let mut result_eff = vec![G1Element::zero(); poly.len()];
        // for i in 0..poly.len()-1{
        //     result_eff[i] = h_eff[i];
        // }

        // self.domain.fft_in_place_g1(&mut result_eff);
        // println!("{:?}",result_eff );

        //test
        let mut test_poly = poly.clone();
        test_poly.reverse();
        test_poly.truncate(degree);
        test_poly.reverse();
        // println!("Original polymnomial coefficients are {:?}",poly);
        let h_test = compute_matrix_vector_multiplication(&test_poly, &t);
        println!("{:?}", h_test);
        let mut result = vec![G1Element::zero(); poly.len()];
        // println!("{:?}",result );
        
        for i in 0..poly.len()-1{
            result[i] = h_test[i];
        }
        
        let mut fin_result = vec![G1Element::zero(); poly.len()];
        for i in 0..degree+1{
            for j in 0..degree{
                let temp = result[j].mul(self.domain.element(i*j));
                fin_result[i] = fin_result[i] + temp;
            }
        }

        // self.domain.fft_in_place_g1(&mut result);
        // println!("{:?}",result);
        fin_result
        

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

        // Perform the pairing check
        lhs.pairing(&G2Element::generator()) == open_i.pairing(&rhs)
    }

    fn update(
        &self,
        commitment: &mut G1Element,
        index: usize,
        old_v_i: &Scalar,
        new_v_i: &Scalar,
    ) -> G1Element {
        *commitment
    }

    fn update_open_i(
        &self,
        open: &mut G1Element,
        index: usize,
        old_v_i: &Scalar,
        new_v_i: &Scalar,
    ) -> G1Element {
        *open
    }

    fn update_open_j(
        &self,
        open: &mut G1Element,
        index: usize,
        index_j: usize,
        old_v_j: &Scalar,
        new_v_j: &Scalar,
    ) -> G1Element {
        *open
    }
}

#[cfg(test)]
mod tests {
    use fastcrypto::groups::bls12381::Scalar;
    use rand::Rng;

    use super::*;

    // #[test]
    // fn test_kzg_commit_open_verify() {
    //     let mut rng = rand::thread_rng();

    //     // Create a new KZGFK struct
    //     let n = 9;
    //     let kzg = KZGFK::new(n).unwrap();

    //     // Create a random vector v
    //     let v: Vec<Scalar> = (0..n).map(|_| OtherScalar::rand(&mut rng)).collect();


    //     // Create a commitment
    //     let commitment = kzg.commit(&v);

    //     // Pick a random index to open
    //     let index = rng.gen_range(0..n);

    //     // Create an opening
    //     let open_value = kzg.open(&v, index);

    //     // Verify the opening
    //     let is_valid = kzg.verify(index, &v[index], &commitment, &open_value);

    //     // Assert that the verification passes
    //     assert!(is_valid, "Verification of the opening should succeed.");
    // }

    #[test]
    fn test_kzg_commit_open_all() {
        let mut rng = rand::thread_rng();

        // Create a new KZGFK struct
        let n = 32;
        let kzg = KZGFK::new(n).unwrap();

        // Create a random vector v
        let v: Vec<Scalar> = (0..n).map(|_| OtherScalar::rand(&mut rng)).collect();

        // // println!("{:?}", v);

        // Create a commitment
        let commitment = kzg.commit(&v);

        // Create array with all indices
        let indices: Vec<usize> = (0..n).collect();

        // Create an opening
        let open_values = kzg.open_all(&v, indices.clone());

        // println!("{:?}", open_values);

        // Verify all openings
        for (i, open_value) in open_values.iter().enumerate() {
            let is_valid = kzg.verify(indices[i], &v[indices[i]], &commitment, open_value);
            if is_valid{
                println!("Verification success for {:?}", i);
            }
            assert!(is_valid, "Verification of the opening should succeed for index {}", i);
        }
    }
}
