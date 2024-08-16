# polydiv
Faster Polynomial Division

Polydiv is using L'Hôpital's rule for faster FFT-based polynomial division, applicable to Zero-Knowledge Proofs and Accumulator schemes. It's surprising that a concept many learn in high school has such significant applications in advanced cryptography. This approach optimizes polynomial divisions, offering notable speedups and enhanced parallelization. The technique supports various proof systems, making it a versatile tool for improving performance in cryptographic protocols like Groth16 and KZG commitments. The innovation lies in efficiently handling 0/0 cases in polynomial division, crucial for cryptographic proofs.
