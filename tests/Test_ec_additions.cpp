/**
 *
 * This file is part of FV-NFLlib
 *
 * Copyright (C) 2016  CryptoExperts
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
 *
 */

#include <gmpxx.h>
#include <chrono>
#include <iostream>
#include <nfl.hpp>
#include <thread>
#include <vector>

/// include the FV homomorphic encryption library
namespace FV {
namespace params {
using poly_t = nfl::poly_from_modulus<uint64_t, 1 << 10, 2480>;
mpz_class plaintextModulus(
    "11579208921035624876269744694940757353008614341529031419553363130886709785"
    "3951");
using gauss_struct = nfl::gaussian<uint16_t, uint64_t, 2>;
using gauss_t = FastGaussianNoise<uint16_t, uint64_t, 2>;
gauss_t fg_prng_sk(8.0, 128, 1 << 14);
gauss_t fg_prng_evk(8.0, 128, 1 << 14);
gauss_t fg_prng_pk(8.0, 128, 1 << 14);
gauss_t fg_prng_enc(8.0, 128, 1 << 14);
}  // namespace params
}  // namespace FV
#include <FV.hpp>

/**
 * Parameters for the NIST P-256 Curve
 */
mpz_class coeff_a = FV::params::plaintextModulus - 3;
mpz_class coeff_b(
    "41058363725152142129326129780047268409114441015993725554835256314039467401"
    "291");
mpz_class order_G(
    "11579208921035624876269744694940757352999695522413576034242225906106851204"
    "4369");
mpz_class G_x(
    "48439561293906451759052585252797914202762949526041747995844080717082404635"
    "286");
mpz_class G_y(
    "36134250956749795798585127919587881956611106672985015071877198253568414405"
    "109");
mpz_class G_z("1");
mpz_class G_t(
    "69187469364232031836548821531971153808731075654725806004116076052366752432"
    "012");

/// Include Elliptic Curve Multiplication
#include "elliptic_curve.hpp"

int main() {
  // Timer
  auto start = std::chrono::steady_clock::now();
  auto end = std::chrono::steady_clock::now();

  std::cout << "Timings:" << std::endl;

  // Precomputation in the clear
  ECPrecomputation<FV::mess_t> precomputation;

  // Point G and variable for G+G
  ECPoint<FV::mess_t> *G = new ECPoint<FV::mess_t>(G_x, G_y, G_z, G_t);
  ECPoint<FV::mess_t> *GpG = new ECPoint<FV::mess_t>();

  // Compute 4G = (G+G)+(G+G) with 2 EC additions
  start = std::chrono::steady_clock::now();
  ec_addition<FV::mess_t>(*GpG, *G, *G, precomputation);
  ec_addition<FV::mess_t>(*GpG, *GpG, *GpG, precomputation);
  end = std::chrono::steady_clock::now();
  std::cout << "\tEC Addition in clear: \t\t"
            << FV::util::get_time_us(start, end, 2) << " us" << std::endl;

  FV::mess_t result_x = *(GpG->X) + FV::params::plaintextModulus;
  FV::mess_t result_y = *(GpG->Y) + FV::params::plaintextModulus;

  // Multiply by the inverse of Z
  FV::mess_t iZ = (*(*GpG).Z).invert();
  *(*GpG).X *= iZ;
  *(*GpG).Y *= iZ;

  /**
   * Let's do it homomorphically
   */

  // Keygen
  start = std::chrono::steady_clock::now();
  FV::sk_t *secret_key = new FV::sk_t();
  end = std::chrono::steady_clock::now();
  std::cout << "\tSecret key generation: \t\t"
            << FV::util::get_time_us(start, end, 1) << " us" << std::endl;

  start = std::chrono::steady_clock::now();
  FV::evk_t *evaluation_key = new FV::evk_t(*secret_key, 32);
  end = std::chrono::steady_clock::now();
  std::cout << "\tEvaluation key generation: \t"
            << FV::util::get_time_us(start, end, 1) << " us" << std::endl;

  start = std::chrono::steady_clock::now();
  FV::pk_t *public_key = new FV::pk_t(*secret_key, *evaluation_key);
  end = std::chrono::steady_clock::now();
  std::cout << "\tPublic key generation: \t\t"
            << FV::util::get_time_us(start, end, 1) << " us" << std::endl;

  // Precomputation
  ECPrecomputation<FV::ciphertext_t> precomputation_encrypted;

  // Point G and variable for G+G
  FV::ciphertext_t *eG_x = new FV::ciphertext_t();
  FV::ciphertext_t *eG_y = new FV::ciphertext_t();
  FV::ciphertext_t *eG_z = new FV::ciphertext_t();
  FV::ciphertext_t *eG_t = new FV::ciphertext_t();

  FV::mess_t *mG_x = new FV::mess_t(G_x);
  FV::mess_t *mG_y = new FV::mess_t(G_y);
  FV::mess_t *mG_z = new FV::mess_t(G_z);
  FV::mess_t *mG_t = new FV::mess_t(G_t);

  start = std::chrono::steady_clock::now();
  FV::encrypt(*eG_x, *public_key, *mG_x);
  FV::encrypt(*eG_y, *public_key, *mG_y);
  FV::encrypt(*eG_z, *public_key, *mG_z);
  FV::encrypt(*eG_t, *public_key, *mG_t);
  end = std::chrono::steady_clock::now();
  std::cout << "\tPoint Encryption: \t\t"
            << FV::util::get_time_us(start, end, 1) << " us" << std::endl;

  ECPoint<FV::ciphertext_t> *eG =
      new ECPoint<FV::ciphertext_t>(*eG_x, *eG_y, *eG_z, *eG_t);
  ECPoint<FV::ciphertext_t> *eGpG = new ECPoint<FV::ciphertext_t>();

  // Compute 4G = (G+G)+(G+G) with 2 EC additions
  start = std::chrono::steady_clock::now();
  ec_addition<FV::ciphertext_t>(*eGpG, *eG, *eG, precomputation_encrypted);
  ec_addition<FV::ciphertext_t>(*eGpG, *eGpG, *eGpG, precomputation_encrypted);
  end = std::chrono::steady_clock::now();
  std::cout << "\tHomom. EC Addition: \t\t"
            << FV::util::get_time_us(start, end, 2) / 1000 << " ms"
            << std::endl;

  start = std::chrono::steady_clock::now();
  FV::decrypt(*mG_x, *secret_key, *public_key, *(*eGpG).X);
  FV::decrypt(*mG_y, *secret_key, *public_key, *(*eGpG).Y);
  FV::decrypt(*mG_z, *secret_key, *public_key, *(*eGpG).Z);
  FV::decrypt(*mG_t, *secret_key, *public_key, *(*eGpG).T);
  end = std::chrono::steady_clock::now();
  std::cout << "\tEC Point Decryption: \t\t"
            << FV::util::get_time_us(start, end, 1) << " us" << std::endl;

  // Noise
  unsigned noise_x = noise(*mG_x, *secret_key, *public_key, *(*eGpG).X);
  std::cout << "noise in ciphertext: \t" << noise_x << "/"
            << public_key->noise_max << std::endl;

  // Multiply by the inverse of Z
  FV::mess_t invZ = (*mG_z).invert();
  *mG_x *= invZ;
  *mG_y *= invZ;

  // Results
  std::cout << "4*G (clear): \t\t(" << *(*GpG).X << "," << *(*GpG).Y << ")"
            << std::endl;
  std::cout << "4*G (enc.): \t\t(" << *mG_x << "," << *mG_y << ")" << std::endl;

  // Clean
  delete secret_key;
  delete evaluation_key;
  delete public_key;

  delete mG_x;
  delete mG_y;
  delete mG_z;
  delete mG_t;

  delete eG_x;
  delete eG_y;
  delete eG_z;
  delete eG_t;

  delete eG;
  delete eGpG;
  delete G;
  delete GpG;

  return 0;
}