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

#pragma once

template <class E>
class ECPrecomputation;

template <>
class ECPrecomputation<FV::mess_t> {
 public:
  FV::mess_t *two, *two_a, *three_a, *two_b, *twelve_b, *four_b, *a_sq,
      *three_a_sq, *a, *two_a_sq, *four_a_b, *a_3_8_b_2;

  ECPrecomputation() {
    two = new FV::mess_t(2);
    two_a = new FV::mess_t(2 * coeff_a);
    three_a = new FV::mess_t(3 * coeff_a);
    two_b = new FV::mess_t(2 * coeff_b);
    twelve_b = new FV::mess_t(12 * coeff_b);
    four_b = new FV::mess_t(4 * coeff_b);
    a_sq = new FV::mess_t(coeff_a * coeff_a);
    three_a_sq = new FV::mess_t(3 * coeff_a * coeff_a);
    a = new FV::mess_t(coeff_a);
    two_a_sq = new FV::mess_t(2 * coeff_a * coeff_a);
    four_a_b = new FV::mess_t(4 * coeff_a * coeff_b);
    a_3_8_b_2 =
        new FV::mess_t(coeff_a * coeff_a * coeff_a + 8 * coeff_b * coeff_b);
  }

  ~ECPrecomputation() {
    delete two;
    delete two_a;
    delete three_a;
    delete two_b;
    delete twelve_b;
    delete four_b;
    delete a_sq;
    delete a;
    delete two_a_sq;
    delete four_a_b;
    delete three_a_sq;
    delete a_3_8_b_2;
  }
};

template <>
class ECPrecomputation<FV::ciphertext_t> {
  using P = FV::ciphertext_t::type;

 public:
  P *two, *two_a, *three_a, *two_b, *twelve_b, *four_b, *a_sq, *three_a_sq, *a,
      *two_a_sq, *four_a_b, *a_3_8_b_2;

  ECPrecomputation() {
    ECPrecomputation<FV::mess_t> pc;

    two = FV::util::alloc_aligned<P, 32>(1, pc.two->getValue());
    two_a = FV::util::alloc_aligned<P, 32>(1, pc.two_a->getValue());
    three_a = FV::util::alloc_aligned<P, 32>(1, pc.three_a->getValue());
    two_b = FV::util::alloc_aligned<P, 32>(1, pc.two_b->getValue());
    twelve_b = FV::util::alloc_aligned<P, 32>(1, pc.twelve_b->getValue());
    four_b = FV::util::alloc_aligned<P, 32>(1, pc.four_b->getValue());
    a_sq = FV::util::alloc_aligned<P, 32>(1, pc.a_sq->getValue());
    three_a_sq = FV::util::alloc_aligned<P, 32>(1, pc.three_a_sq->getValue());
    a = FV::util::alloc_aligned<P, 32>(1, pc.a->getValue());
    two_a_sq = FV::util::alloc_aligned<P, 32>(1, pc.two_a_sq->getValue());
    four_a_b = FV::util::alloc_aligned<P, 32>(1, pc.four_a_b->getValue());
    a_3_8_b_2 = FV::util::alloc_aligned<P, 32>(1, pc.a_3_8_b_2->getValue());

    // Store in Double-CRT form
    two->ntt_pow_phi();
    two_a->ntt_pow_phi();
    three_a->ntt_pow_phi();
    two_b->ntt_pow_phi();
    twelve_b->ntt_pow_phi();
    four_b->ntt_pow_phi();
    a_sq->ntt_pow_phi();
    three_a_sq->ntt_pow_phi();
    a->ntt_pow_phi();
    two_a_sq->ntt_pow_phi();
    four_a_b->ntt_pow_phi();
    a_3_8_b_2->ntt_pow_phi();
  }

  ~ECPrecomputation() {
    FV::util::free_aligned(1, two);
    FV::util::free_aligned(1, two_a);
    FV::util::free_aligned(1, three_a);
    FV::util::free_aligned(1, two_b);
    FV::util::free_aligned(1, twelve_b);
    FV::util::free_aligned(1, four_b);
    FV::util::free_aligned(1, a_sq);
    FV::util::free_aligned(1, a);
    FV::util::free_aligned(1, two_a_sq);
    FV::util::free_aligned(1, four_a_b);
    FV::util::free_aligned(1, three_a_sq);
    FV::util::free_aligned(1, a_3_8_b_2);
  }
};

template <class E>
struct ECPoint {
  E *X, *Y, *Z, *T;
  ECPoint(const E &x, const E &y, const E &z, const E &t) {
    X = new E(x);
    Y = new E(y);
    Z = new E(z);
    T = new E(t);
  }
  ECPoint() {
    X = new E();
    Y = new E();
    Z = new E();
    T = new E();
  }
  ~ECPoint() {
    delete X;
    delete Y;
    delete Z;
    delete T;
  }
};

// C = A+B
template <class E>
void ec_addition(ECPoint<E> &C, const ECPoint<E> &A, const ECPoint<E> &B,
                 const ECPrecomputation<E> &precomputation) {
  E *X1X2 = new E(*A.X), *X1Z2 = new E(*A.X), *X1T2 = new E(*A.X);
  E *Y1Y2 = new E(*A.Y), *Y1Z2 = new E(*A.Y), *Y1T2 = new E(*A.Y);
  E *Z1X2 = new E(*A.Z), *Z1Y2 = new E(*A.Z), *Z1Z2 = new E(*A.Z);
  E *T1X2 = new E(*A.T), *T1Y2 = new E(*A.T), *T1T2 = new E(*A.T);

  *X1X2 *= *B.X;
  *X1Z2 *= *B.Z;
  *X1T2 *= *B.T;

  *Y1Y2 *= *B.Y;
  *Y1Z2 *= *B.Z;
  *Y1T2 *= *B.T;

  *Z1X2 *= *B.X;
  *Z1Y2 *= *B.Y;
  *Z1Z2 *= *B.Z;

  *T1X2 *= *B.X;
  *T1Y2 *= *B.Y;
  *T1T2 *= *B.T;

  E *tmp = new E();

  E *F = new E(*T1T2);
  *tmp = *X1X2;
  *tmp *= *precomputation.two_a;
  *F -= *tmp;
  *tmp = *X1Z2;
  *tmp += *Z1X2;
  *tmp *= *precomputation.four_b;
  *F -= *tmp;
  *tmp = *Z1Z2;
  *tmp *= *precomputation.a_sq;
  *F += *tmp;

  E *H = new E(*X1T2);
  *H += *T1X2;
  *tmp = *Y1Y2;
  *tmp *= *precomputation.two;
  *H += *tmp;
  *tmp = *X1Z2;
  *tmp += *Z1X2;
  *tmp *= *precomputation.a;
  *H += *tmp;
  *tmp = *Z1Z2;
  *tmp *= *precomputation.two_b;
  *H += *tmp;

  E *GA1 = new E(*T1T2);
  *tmp = *X1X2;
  *tmp *= *precomputation.two_a;
  *GA1 += *tmp;
  *tmp = *Z1X2;
  *tmp *= *precomputation.four_b;
  *GA1 += *tmp;
  *tmp = *X1Z2;
  *tmp *= *precomputation.twelve_b;
  *GA1 += *tmp;
  *tmp = *Z1Z2;
  *tmp *= *precomputation.three_a_sq;
  *GA1 -= *tmp;

  E *GA2 = new E(*T1T2);
  *tmp = *X1X2;
  *tmp *= *precomputation.two_a;
  *GA2 += *tmp;
  *tmp = *Z1X2;
  *tmp *= *precomputation.four_b;
  *GA2 += *tmp;
  *tmp = *X1Z2;
  *tmp *= *precomputation.twelve_b;
  *GA2 += *tmp;
  *tmp = *Z1Z2;
  *tmp *= *precomputation.three_a_sq;
  *GA2 -= *tmp;

  E *GB1 = new E(0);
  *tmp = *X1T2;
  *tmp *= *precomputation.four_b;
  *GB1 += *tmp;
  *tmp = *X1X2;
  *tmp *= *precomputation.two_a_sq;
  *GB1 -= *tmp;
  *tmp = *X1Z2;
  *tmp *= *precomputation.four_a_b;
  *GB1 -= *tmp;
  *tmp = *T1T2;
  *tmp *= *precomputation.three_a;
  *GB1 += *tmp;
  *tmp = *Z1Z2;
  *tmp *= *precomputation.a_3_8_b_2;
  *GB1 -= *tmp;

  E *GB2 = new E(0);
  *tmp = *T1X2;
  *tmp *= *precomputation.four_b;
  *GB2 += *tmp;
  *tmp = *X1X2;
  *tmp *= *precomputation.two_a_sq;
  *GB2 -= *tmp;
  *tmp = *Z1X2;
  *tmp *= *precomputation.four_a_b;
  *GB2 -= *tmp;
  *tmp = *T1T2;
  *tmp *= *precomputation.three_a;
  *GB2 += *tmp;
  *tmp = *Z1Z2;
  *tmp *= *precomputation.a_3_8_b_2;
  *GB2 -= *tmp;

  E *G1 = new E(*T1Y2);
  *G1 *= *GA1;

  E *G2 = new E(*Z1Y2);
  *G2 *= *GB1;

  E *G3 = new E(*Y1T2);
  *G3 *= *GA2;

  E *G4 = new E(*Y1Z2);
  *G4 *= *GB2;

  *C.X = *F;
  *C.X *= *H;

  *C.Y = *G1;
  *C.Y += *G2;
  *C.Y += *G3;
  *C.Y += *G4;

  *C.Z = *H;
  *C.Z *= *H;

  *C.T = *F;
  *C.T *= *F;

  // Clean
  delete tmp;
  delete F;
  delete H;
  delete GA1;
  delete GA2;
  delete GB1;
  delete GB2;
  delete G1;
  delete G2;
  delete G3;
  delete G4;
}

#if 0

#Sage verification code

p = 2^256-2^224+2^192+2^96-1
proof.arithmetic(False) # turn off primality checking
F = GF(p)
n = 115792089210356248762697446949407573529996955224135760342422259061068512044369
a = F(-3)
b = F(41058363725152142129326129780047268409114441015993725554835256314039467401291)
E = EllipticCurve([a, b])
G = E(48439561293906451759052585252797914202762949526041747995844080717082404635286,36134250956749795798585127919587881956611106672985015071877198253568414405109)
G.set_order(115792089210356248762697446949407573529996955224135760342422259061068512044369)

#endif
