use std::ops::Mul;
use fastcrypto::error::{FastCryptoResult,FastCryptoError};
use fastcrypto::groups::bls12381::{G1Element, G2Element, Scalar};
use fastcrypto::groups::{GroupElement, MultiScalarMul, Pairing, Scalar as OtherScalar};
use fastcrypto::serde_helpers::ToFromByteArray;
use itertools::{iterate, Itertools};
use rand::thread_rng;
use rayon::prelude::*;
use std::sync::{Mutex, Arc};
extern crate nalgebra as na;
use na::{DMatrix};




use crate::fft::{BLS12381Domain, FFTDomain};
use crate::KZG;

// fn inverse_matrix(matrix: Vec<Vec<Scalar>>) -> Vec<Vec<Scalar>> {
//     let n = matrix.len();
//     let mut augmented = vec![vec![Scalar::zero(); 2 * n]; n];

//     // Create the augmented matrix [matrix | identity]
//     for i in 0..n {
//         for j in 0..n {
//             augmented[i][j] = matrix[i][j].clone();
//         }
//         augmented[i][n + i] = Scalar::from((1) as u128);
//     }

//     // Perform Gaussian elimination
//     for i in 0..n {
//         // Find the pivot
//         let mut pivot = i;
//         for j in (i + 1)..n {
//             if augmented[j][i] != Scalar::zero() {
//                 pivot = j;
//                 break;
//             }
//         }
//         if augmented[pivot][i] == Scalar::zero() {
//             return vec![];; // Matrix is singular
//         }

//         // Swap rows i and pivot
//         augmented.swap(i, pivot);

//         // Normalize the pivot row
//         let pivot_inv = augmented[i][i].inverse().unwrap(); // Ensure invertibility
//         for j in 0..2 * n {
//             augmented[i][j] *= pivot_inv;
//         }

//         // Eliminate the current column in all other rows
//         for j in 0..n {
//             if j != i {
//                 let factor = augmented[j][i].clone();
//                 for k in 0..2 * n {
//                     augmented[j][k] -= factor * augmented[i][k].clone();
//                 }
//             }
//         }
//     }

//     // Extract the inverse matrix from the augmented matrix
//     let mut inverse = vec![vec![Scalar::zero(); n]; n];
//     for i in 0..n {
//         for j in 0..n {
//             inverse[i][j] = augmented[i][n + j].clone();
//         }
//     }

//     inverse
// }


// fn convert_to_dmatrix(matrix_data: Vec<Vec<Scalar>>) -> DMatrix<Scalar> {
//     let rows = matrix_data.len();
//     let cols = matrix_data[0].len();
//     DMatrix::from_vec(rows, cols, matrix_data.into_iter().flatten().collect())
// }

fn add_vectors(v1: Vec<G1Element>, v2: Vec<G1Element>, v3: Vec<G1Element>) -> Vec<G1Element> {
    v1.into_iter().zip(v2.into_iter()).zip(v3.into_iter())
        .map(|((a, b), c)| a+b+c)
        .collect()
}

fn matrix_multiply(a: &Vec<Vec<Scalar>>, b: &Vec<Vec<Scalar>>) -> Vec<Vec<Scalar>> {
    // Check if the matrices can be multiplied
    // if a[0].len() != b.len() {
    //     return; // The number of columns in A must be equal to the number of rows in B
    // }

    let rows_a = a.len();
    let cols_a = a[0].len();
    let cols_b = b[0].len();

    // Initialize the result matrix with zeros
    let mut result = vec![vec![Scalar::zero(); cols_b]; rows_a];

    // Perform matrix multiplication
    for i in 0..rows_a {
        for j in 0..cols_b {
            for k in 0..cols_a {
                result[i][j] += a[i][k] * b[k][j];
            }
        }
    }

    result
}

fn multiply_matrix_scalar_vector(matrix: &[Vec<Scalar>], v: &[Scalar]) -> Vec<Scalar> {
    let n = matrix.len();
    let mut result = vec![Scalar::zero(); n];

    for i in 0..n {
        for j in 0..v.len() {
            if matrix[i][j] != Scalar::zero() {
                result[i] = result[i] + matrix[i][j].mul(v[j]);
            }
        }
    }

    result
}

fn multiply_d_matrix_by_vector(vector: &[Scalar]) -> Vec<Scalar> {
    let n = vector.len();
    let mut result = vec![Scalar::zero(); n];

    for i in 0..n-1 {
        result[i] = Scalar::from((i + 1) as u128) * vector[i + 1];
    }

    result
}
// this currently does n MSMs, can improve this to an O(n) algorithm as we did for d_matrix multiplicaiton
fn multiply_matrix_vector(matrix: &[Vec<Scalar>], v: &[G1Element]) -> Vec<G1Element> {
    let n = v.len();
    let mut result = vec![G1Element::zero(); n];

    for i in 0..n {
        let scalars = &matrix[i];
        result[i] = G1Element::multi_scalar_mul(&scalars, &v).unwrap();
    }

    result
}




#[derive(Clone)]
pub struct KZGDeriv {
    domain: BLS12381Domain,
    n: usize,
    omega: Scalar,
    g2_tau: G2Element,
    w_vec: Vec<G1Element>,
    check_w_vec: Vec<G1Element>,
    u_vec: Vec<G1Element>,
    tau_powers_g1: Vec<G1Element>,
}

impl KZGDeriv {
    fn element(&self, index: usize) -> Scalar {
        // Compute omega^index recursively

        if index > self.n {
            return self.element(index & self.n);
        } else if index == 0 {
            return Scalar::generator();
        } else if index == 1 {
            return self.omega;
        }

        let half = self.element(index / 2);
        if index % 2 == 0 {
            return half * half;
        }
        half * half * self.omega
    }
}

impl KZG for KZGDeriv {
    // Uses the BLS12-381 construction
    type G = G1Element;

    fn new(n: usize) -> FastCryptoResult<Self> {
        let domain = BLS12381Domain::new(n)?;

        // Generate tau using a random scalar
        let tau = Scalar::rand(&mut thread_rng());
        let g2_tau = G2Element::generator() * tau;

        // Compute tau^i for i = 0 to n-1
        let tau_powers_g: Vec<Scalar> =
            iterate(Scalar::generator(), |g| g * tau).take(n).collect();

        let g = G1Element::generator();
        // Compute w_vec and u_i = g^{(L_i(tau) - 1)/(tau-omega^i)}
        let w_vec: Vec<G1Element> = domain
            .ifft(&tau_powers_g)
            .iter()
            .map(|s| g.mul(s))
            .collect();

        let mut omega_i = Scalar::generator();
        let mut u_vec = (0..n)
            .map(|i| {
                // Compute u_i
                let l_i_minus_1 = w_vec[i] - G1Element::generator();
                let denom = tau - omega_i;
                let u_i = (l_i_minus_1 / denom).unwrap();

                // Update omega_i
                if i < n - 1 {
                    omega_i *= domain.element(1);
                }
                u_i
            })
            .collect();

        let omega = domain.element(1);

        let tau_powers_g1: Vec<G1Element> = tau_powers_g.iter().map(|x| g.mul(x)).collect();

        let mut check_w_vec = tau_powers_g1.clone();
        domain.ifft_in_place_g1(&mut check_w_vec);
        println!("w_vec {:?}", w_vec);
        println!("check_w_vec{:?}", check_w_vec);

        println!("This test is {}", check_w_vec == w_vec);

        Ok(Self {
            domain,
            n,
            omega,
            g2_tau,
            w_vec,
            check_w_vec,
            u_vec,
            tau_powers_g1
        })
    }

    fn commit(&self, v: &[Scalar]) -> G1Element {
        G1Element::multi_scalar_mul(&v, &self.w_vec).unwrap()
    }

    fn open(&self, v: &[Scalar], index: usize) -> G1Element {
        let scalars = Arc::new(Mutex::new(vec![Scalar::zero(); v.len()]));
        let v_prime = Arc::new(Mutex::new(Scalar::zero()));

        // Initialize omega_i and omega_j_minus_i
        let omega_i = self.element(index);
        let omega_base = Scalar::generator();
        let omega_j_minus_i_base = self.element(self.n - index);

        (0..v.len()).into_par_iter().for_each(|j| {
            let mut local_omega_j = omega_base;
            let mut local_omega_j_minus_i = omega_j_minus_i_base;

            // Compute local omegas
            for _ in 0..j {
                local_omega_j *= self.omega;
                local_omega_j_minus_i *= self.omega;
            }

            if j != index {
                let diff_inverse = (omega_i - local_omega_j).inverse().unwrap();
                let mut v_prime_guard = v_prime.lock().unwrap();
                *v_prime_guard += v[j] * local_omega_j_minus_i * diff_inverse;

                let mut scalars_guard = scalars.lock().unwrap();
                scalars_guard[j] = (v[index] - v[j]) * diff_inverse;
            } else {
                let mut v_prime_guard = v_prime.lock().unwrap();
                *v_prime_guard += v[j]
                    * (Scalar::from((v.len() - 1) as u128) / (Scalar::from(2u128) * omega_i))
                        .unwrap();
            }
        });

        let v_prime = v_prime.lock().unwrap().clone();
        let mut scalars_guard = scalars.lock().unwrap();
        scalars_guard[index] = v_prime;

        G1Element::multi_scalar_mul(&scalars_guard, &self.w_vec[..scalars_guard.len()]).unwrap()


        // implement open with edivk

    }

    fn open_all(&self, v: &[Scalar], indices: Vec<usize>) -> Vec<G1Element> {

        let n = v.len();

        //compute D - dont need this
        let mut d_matrix = vec![vec![Scalar::zero(); n]; n];
        for i in 0..n-1 { // n-1 to avoid out of bounds access for j = i+1
            d_matrix[i][i + 1] = Scalar::from((i + 1) as u128);
        }

        

        // println!("{:?}",d_matrix);

        //compute ColEDivHat
        let mut c_matrix = vec![vec![Scalar::zero(); n]; n];
        for i in 0..n {
            for j in 0..n {
                if i == 0 && j == n - 1 {
                    c_matrix[i][j] = -((Scalar::from((n - 1) as u128)) / (Scalar::from(2 as u128))).unwrap();
                } else if j + 1 == i && i >= 1 && i < n {
                    c_matrix[i][j] = -Scalar::from(i as u128) + (Scalar::from((n + 1) as u128) / Scalar::from(2u128)).unwrap();
                } else {
                    c_matrix[i][j] = Scalar::zero();
                }
            }
        }

        let mut d_matrix = vec![vec![Scalar::zero(); n]; n];
        for i in 0..n {
            for j in 0..n {
                if i == 0 && j == n - 1 {
                    d_matrix[i][j] = ((Scalar::from((n - 1) as u128)) / (Scalar::from(2 as u128))).unwrap();
                } else if j + 1 == i && i >= 1 && i < n {
                    d_matrix[i][j] = Scalar::from(i as u128) - (Scalar::from((n + 1) as u128) / Scalar::from(2u128)).unwrap();
                } else {
                    d_matrix[i][j] = Scalar::zero();
                }
            }
        }

        // println!("{:?}",c_matrix);


        


        //first compute tau*Dhatv - result1
        let idftv = self.domain.ifft(&v);
        let d_msm_idftv:Vec<Scalar> = multiply_d_matrix_by_vector(&idftv);
        let dhatv = self.domain.fft(&d_msm_idftv);
        let result1: Vec<G1Element> = self.w_vec.iter()
            .zip(dhatv.iter())
            .map(|(a, b)| a.mul(*b)) 
            .collect();
        
        //next compute ColEDiv.tau*v - result2

        let mut powtau = self.w_vec.clone();
        self.domain.fft_in_place_g1(&mut powtau);
        let mut col_hat_dft_tau: Vec<G1Element> = multiply_matrix_vector(&c_matrix, &powtau);
        self.domain.ifft_in_place_g1(&mut col_hat_dft_tau);
        let result2: Vec<G1Element> = col_hat_dft_tau.iter()
            .zip(v.iter())
            .map(|(a, b)| a.mul(*b)) 
            .collect();
        
        // finally compute diadiv.powtau*v - result3

        let mut mult:Vec<G1Element> = self.w_vec.iter()
        .zip(v.iter())
        .map(|(a, b)| a.mul(*b)) 
        .collect();

        self.domain.fft_in_place_g1(&mut mult);
        let mut diadiv_idft_tau_v : Vec<G1Element> = multiply_matrix_vector(&d_matrix, &mult);
        self.domain.ifft_in_place_g1(&mut diadiv_idft_tau_v);

        let  result3 = diadiv_idft_tau_v.clone();

        let result = add_vectors(result1.clone(),result2.clone(),result3.clone());

        


        

        //test code to check the equation is correct

        //compute result1 
        
        let mut d_hat = vec![vec![Scalar::zero(); n]; n];
        for i in 0..n{
            for j in 0..n{
                let omega_i = self.element(i);
                let omega_j = self.element(j);
                if i == j{
                    d_hat[i][j] = (Scalar::from((n - 1) as u128) / (Scalar::from(2u128) * omega_i)).unwrap();
                }
                else{
                    d_hat[i][j] =  self.element(j) * self.element(n - i) * (omega_i - omega_j).inverse().unwrap();
                }
            }
        }
        // println!("{:?}", d_hat);

        let d_hat_v = multiply_matrix_scalar_vector(&d_hat, v);

        let result1_test: Vec<G1Element> = self.w_vec.iter()
            .zip(d_hat_v.iter())
            .map(|(a, b)| a.mul(*b)) 
            .collect();

        println!("Check result1 {}", result1 == result1_test);


        //compute result2 

        let mut c_div = vec![vec![Scalar::zero(); n]; n];
        for i in 0..n{
            for j in 0..n{
                let omega_i = self.element(i);
                let omega_j = self.element(j);
                if i == j{
                    c_div[i][j] = Scalar::zero();
                }
                else{
                    c_div[i][j] = (omega_i - omega_j).inverse().unwrap();
                }
            }
        }

        let mut d_div = vec![vec![Scalar::zero(); n]; n];
        for i in 0..n{
            for j in 0..n{
                let omega_i = self.element(i);
                let omega_j = self.element(j);
                if i == j{
                    d_div[i][j] = Scalar::zero();
                }
                else{
                    d_div[i][j] = (omega_j - omega_i).inverse().unwrap();
                }
            }
        }

        // println!("{:?}", c_div);

        let mut col_pow_tau = multiply_matrix_vector(&c_div, &self.w_vec);

        // println!("{:?}", col_pow_tau);

        let result2_test: Vec<G1Element> = col_pow_tau.iter()
            .zip(v.iter())
            .map(|(a, b)| a.mul(*b)) 
            .collect();
            // println!("{:?}", result2);
        
        println!("Check result2 {}", result2 == result2_test);
        
        //compute result3 

        let mut mult_test:Vec<G1Element> = self.w_vec.iter()
        .zip(v.iter())
        .map(|(a, b)| a.mul(*b)) 
        .collect();

        let result3_test = multiply_matrix_vector(&d_div, &mult_test);

        let result_test = add_vectors(result1_test,result2_test,result3_test.clone());

        println!("Check result3 {}", result3 == result3_test);
        println!("{}", result == result_test);
        // result_test
        result

        // // test with edivk 

        // let mut result_ediv = vec![G1Element::zero();v.len()];
        // let mut edivk = vec![vec![Scalar::zero(); n]; n];
        // for k in 0..n{
        //     for i in 0..n{
        //         for j in 0..n{
        //             let omega_i = self.element(i);
        //             let omega_j = self.element(j);
        //             let omega_k = self.element(k);
        //             if i == j && i != k{
        //                 edivk[i][j] = (omega_i - omega_k).inverse().unwrap();
        //             }
        //             else if i == k && j != k{
        //                 edivk[i][j] = omega_j * self.element(n-k)*(omega_k - omega_j).inverse().unwrap();
        //             }

        //             else if j == k && i != k{
        //                 edivk[i][j] = (omega_k - omega_i).inverse().unwrap();

        //             }
        //             else if i == j && i == k{
        //                 edivk[i][j] = self.element(n-k)*((Scalar::from((n-1) as u128)) / (Scalar::from(2 as u128))).unwrap();
        //             }
        //             else{
        //                 edivk[i][j] = Scalar::zero();
        //             }
        //         }
        //     }
        //     let edivk_v = multiply_matrix_scalar_vector(&edivk, &v);

        //     result_ediv[k]=G1Element::multi_scalar_mul(&edivk_v, &self.tau_powers_g1).unwrap();
        // }

        

        // println!("{:?}", edivk);

        // let mut m_hat = vec![vec![Scalar::zero(); n]; n];
        // for i in 0..n{
        //     for j in 0..n{
        //         let omega_i = self.element(i);
        //         let omega_j = self.element(j);
        //         if i != j{
        //             m_hat[i][j] = - omega_j*((Scalar::from((1) as u128)) / (Scalar::from(n as u128))).unwrap();
        //         }
        //         else{
        //             m_hat[i][j] = omega_j*((Scalar::from((n-1) as u128)) / (Scalar::from(n as u128))).unwrap();
        //         }
        //     }
        // }

        // for k in 0..n{

        //     let mut e_k = vec![vec![Scalar::zero(); n]; n];
        //     for i in 0..n{
        //         for j in 0..n{
        //             if i == 0 && j == k{
        //                 e_k[i][j] = 1.into();
        //             }
        //             else{
        //                 e_k[i][j] = 0.into();
        //             }
        //         }
        //     }

        //     let mut i_matrix = vec![vec![Scalar::zero(); n]; n];
        //     for i in 0..n{
        //         for j in 0..n{
        //             if i == j{
        //                 i_matrix[i][j] = 1.into();
        //             }
        //             else{
        //                 i_matrix[i][j] = 0.into();
        //             }
        //         }
        //     }

        //     let omega_k = self.element(k);

        //     let mut omega_k_i = vec![vec![Scalar::zero(); n]; n];

        //     for i in 0..n{
        //         for j in  0..n{
        //             omega_k_i[i][j] = omega_k*i_matrix[i][j];
        //         }
        //     }

        //     let mut m_hat_minus_omega_k_i = vec![vec![Scalar::zero(); n]; n];
        //     for i in 0..n{
        //         for j in 0..n{
        //             m_hat_minus_omega_k_i[i][j] = m_hat[i][j]-omega_k*i_matrix[i][j];
        //         }
        //     }

        //     let matrix = convert_to_dmatrix(m_hat_minus_omega_k_i);
        //     let inverse_vec = inverse_matrix(m_hat_minus_omega_k_i);
        //     // match matrix.try_inverse() {
        //     //     Some(inverse) => {
        //     //         // Convert the inverse DMatrix back to a vector of vectors
        //     //         inverse_vec = (0..inverse.nrows())
        //     //             .map(|i| (0..inverse.ncols()).map(|j| inverse[(i, j)].clone()).collect())
        //     //             .collect();
        
        //     //         println!("The inverse of the matrix is:");
        //     //         for row in &inverse_vec {
        //     //             println!("{:?}", row);
        //     //         }
        //     //     }
        //     //     None => {
        //     //         println!("The matrix is not invertible.");
        //     //     }
        //     // }

        //     let mut dft_ek = self.domain.fft(&e_k);
        //     let mut right_matrix = vec![vec![Scalar::zero(); n]; n];
        //     for i in 0..n{
        //         for j in 0..n{
        //             right_matrix[i][j] = i_matrix[i][j] - dft_ek[i][j];
        //         }
        //     }

        //     let mut edivk = matrix_multiply(&inverse_vec, &right_matrix);

        //     let ediv_vecval = multiply_matrix_scalar_vector(&edivk, &v);

        //     result_ediv[k] = G1Element::multi_scalar_mul(&ediv_vecval, &self.tau_powers_g1).unwrap();








        // }

        // result_ediv
        


        

        // indices.into_iter().map(|i| self.open(v, i)).collect()
    }

    fn verify(
        &self,
        index: usize,
        v_i: &Scalar,
        commitment: &G1Element,
        open_i: &G1Element,
    ) -> bool {
        let lhs = *commitment - G1Element::generator() * v_i;
        let rhs = self.g2_tau - G2Element::generator() * self.element(index);

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
        *commitment + self.w_vec[index].mul(new_v_i - old_v_i)
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
        assert_ne!(index, index_j, "index and index_j should be different.");

        let omega_i = self.element(index);
        let omega_j = self.element(index_j);
        let omega_i_inverse = self.element(self.n - index);

        let to_mul_1 = (old_v_j - new_v_j) * (omega_i - omega_j).inverse().unwrap();
        let to_mul_2 = -to_mul_1 * omega_j * omega_i_inverse;

        *open
            + G1Element::multi_scalar_mul(
                &[to_mul_1, to_mul_2],
                &[self.w_vec[index_j], self.w_vec[index]],
            )
            .unwrap()
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
        let n = 4;
        let kzg = KZGDeriv::new(n).unwrap();

        // Create a random vector v
        let v: Vec<Scalar> = (0..n).map(|_| OtherScalar::rand(&mut rng)).collect();

        // println!("{:?}", v);

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

        // Create a new KZGDeriv struct
        let n = 8;
        let kzg = KZGDeriv::new(n).unwrap();

        // Create a random vector v
        let v: Vec<Scalar> = (0..n).map(|_| OtherScalar::rand(&mut rng)).collect();

        // println!("{:?}", v);

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

        // Create a new KZGDeriv struct
        let n = 8;
        let kzg = KZGDeriv::new(n).unwrap();

        // Create a random vector v
        let v: Vec<Scalar> = (0..n).map(|_| OtherScalar::rand(&mut rng)).collect();

        // println!("{:?}", v);

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
    #[test]
    fn test_kzg_commit_open_all() {
        let mut rng = rand::thread_rng();

        // Create a new KZGDeriv struct
        let n = 8;
        let kzg = KZGDeriv::new(n).unwrap();

        // Create a random vector v
        let v: Vec<Scalar> = (0..n).map(|_| OtherScalar::rand(&mut rng)).collect();

        // println!("{:?}", v);

        // Create a commitment
        let commitment = kzg.commit(&v);

        // Create array with all indices
        let indices: Vec<usize> = (0..n).collect();

        // Create an opening
        let open_values = kzg.open_all(&v, indices.clone());

        println!("{:?}", open_values);

        // Verify all openings
        for (i, open_value) in open_values.iter().enumerate() {
            let is_valid = kzg.verify(indices[i], &v[indices[i]], &commitment, open_value);
            assert!(is_valid, "Verification of the opening should succeed for index {}", i);
        }
    }

}
