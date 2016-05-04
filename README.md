# FV-NFLlib
A header-only library implementing the Fan-Vercauteren homomorphic encryption scheme

## What is it?

FV-NFLlib is a software library implementing an [homomorphic encryption](https://en.wikipedia.org/wiki/Homomorphic_encryption) scheme (HE).
FV-NFLlib implements the Fan-Vercauteren ([FV](http://eprint.iacr.org/2012/144)) scheme, and is based on the [NFLlib](https://github.com/quarkslab/NFLlib) C++ library dedicated to ideal lattice cryptography.

## License

GPLv3

# Install Steps

To use the FV-NFLlib library, you need GMP, Mpfr, a C+11 compiler and the NFLlib library. The library only consists of a single header file **FV.hpp** that has to be included during the compilation of your program.

# Tests

Three tests are provided for the library

- **tests/Test_binary_tree**: this program benchmarks the key generation, the encryption/decryption, homomorphic addition and homomorphic multiplication procedures. A binary tree of multiplications is evaluated to ensure correctness and the final noise bound is given.
- **tests/Test_ec_additions**: this program computes in the clear and homomorphically an elliptic curve addition over the NIST P-256 curve. The code of the elliptic curve addition is templated, and called twice (once with FV::mess_t and once with FV::ciphertext_t)
- **tests/Test_encrypt_poly**: this program encrypts two polynomials and computes their product homomorphically. It outputs a small SAGE program to test the correctness.

## Contributors

This research-oriented library has been done by members of [CryptoExperts](https://www.cryptoexperts.com), and is part of the [HEAT](https://heat-project.eu/) project, and the CryptoComp project.
