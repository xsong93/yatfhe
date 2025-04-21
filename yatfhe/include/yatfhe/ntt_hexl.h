//
// Created by Xintong on 25-4-16.
//

#ifndef HLS_YATFHE_NTT_HEXL_H
#define HLS_YATFHE_NTT_HEXL_H

#include <hexl/hexl.hpp>
#include "yatfhe/polynomial.h"


using namespace intel::hexl;

namespace NttHexl {

    NTT &getNttHexl();

    void initNttHexl(uint64_t degree, uint64_t q);

    void printHexlParams();

    NttPolynomial& getNttRoterPoly(int32_t r);

    void initNttRotMap(int32_t degree);

    void applyNtt(NttPolynomial &out, const TorusPolynomial &in);

    void applyIntt(TorusPolynomial &out, const NttPolynomial &in);

    void calModularInnerProductNtt(NttPolynomial& res, const vector<NttPolynomial>& in1, const vector<NttPolynomial>& in2);

    void calModularInnerProductNtt(NttPolynomial &acc, const NttPolynomial &in1, const NttPolynomial &in2);

    void EltwiseSignedMultMod(uint64_t* result, const uint64_t* operand1,
                    const uint64_t* operand2, int32_t sign, uint64_t n, uint64_t modulus,
                    uint64_t input_mod_factor);

    template<typename T, typename R>
    void applyNttForAB(T &out, R &in) {
        for (auto row = 0; row < in.a.size(); row++) {
            applyNtt(out.a[row], in.a[row]);
        }
        applyNtt(out.b, in.b);
    }

    template<typename T, typename R>
    void applyInttForAB(T &out, R &in) {
        for (auto row = 0; row < in.a.size(); row++) {
            applyIntt(out.a[row], in.a[row]);
        }
        applyIntt(out.b, in.b);
    }

    template <int InputModFactor>
    void EltwiseSignedMultModNative(uint64_t* result, const uint64_t* operand1,
                              const uint64_t* operand2, const int32_t sign, uint64_t n,
                              uint64_t modulus) {
      HEXL_CHECK(InputModFactor == 1 || InputModFactor == 2 || InputModFactor == 4,
                 "Require InputModFactor = 1, 2, or 4")
      HEXL_CHECK(result != nullptr, "Require result != nullptr");
      HEXL_CHECK(operand1 != nullptr, "Require operand1 != nullptr");
      HEXL_CHECK(operand2 != nullptr, "Require operand2 != nullptr");
      HEXL_CHECK(n != 0, "Require n != 0");
      HEXL_CHECK(modulus > 1, "Require modulus > 1");
      HEXL_CHECK(modulus < (1ULL << 62), "Require modulus < (1ULL << 62)");
      HEXL_CHECK_BOUNDS(operand1, n, InputModFactor * modulus,
                        "operand1 exceeds bound " << (InputModFactor * modulus));
      HEXL_CHECK_BOUNDS(operand2, n, InputModFactor * modulus,
                        "operand2 exceeds bound " << (InputModFactor * modulus));

      constexpr int64_t beta = -2;
      HEXL_CHECK(beta <= -2, "beta must be <= -2 for correctness");

      constexpr int64_t alpha = 62;  // ensures alpha - beta = 64

      uint64_t gamma = Log2(InputModFactor);
      HEXL_UNUSED(gamma);
      HEXL_CHECK(alpha >= gamma + 1, "alpha must be >= gamma + 1 for correctness");

      const uint64_t ceil_log_mod = Log2(modulus) + 1;  // "n" from Algorithm 2
      uint64_t prod_right_shift = ceil_log_mod + beta;

      // Barrett factor "mu"
      // TODO(fboemer): Allow MultiplyFactor to take bit shifts != 64
      HEXL_CHECK(ceil_log_mod + alpha >= 64, "ceil_log_mod + alpha < 64");
      uint64_t barr_lo =
          MultiplyFactor(uint64_t(1) << (ceil_log_mod + alpha - 64), 64, modulus)
              .BarrettFactor();

      const uint64_t twice_modulus = 2 * modulus;

      HEXL_LOOP_UNROLL_4
      for (size_t i = 0; i < n; ++i) {
        uint64_t prod_hi, prod_lo, c2_hi, c2_lo, Z;

        uint64_t x = ReduceMod<InputModFactor>(*operand1, modulus, &twice_modulus);
        uint64_t y = ReduceMod<InputModFactor>(*operand2, modulus, &twice_modulus);

        // Multiply inputs
        MultiplyUInt64(x, y, &prod_hi, &prod_lo);

        // floor(U / 2^{n + beta})
        uint64_t c1 = (prod_lo >> (prod_right_shift)) +
                      (prod_hi << (64 - (prod_right_shift)));

        // c2 = floor(U / 2^{n + beta}) * mu
        MultiplyUInt64(c1, barr_lo, &c2_hi, &c2_lo);

        // alpha - beta == 64, so we only need high 64 bits
        uint64_t q_hat = c2_hi;

        // only compute low bits, since we know high bits will be 0
        Z = prod_lo - q_hat * modulus;

        // Conditional subtraction
        *result = (Z >= modulus) ? (Z - modulus) : Z;
        if (sign == -1) {
            *result = modulus - *result;
        }

        ++operand1;
        ++operand2;
        ++result;
      }
    }
}

#endif //HLS_YATFHE_NTT_HEXL_H
