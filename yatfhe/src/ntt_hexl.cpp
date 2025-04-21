//
// Created by Xintong on 25-4-16.
//

#include "yatfhe/ntt_hexl.h"
#include "yatfhe/polynomial.h"
#include "yautil/tool.h"

using namespace intel::hexl;

namespace NttHexl {

    NTT &getNttHexl() {
        static NTT nttHexl;
        return nttHexl;
    }

    void initNttHexl(uint64_t degree, uint64_t q) {
        static bool initialized = false;
        if (!initialized) {
            getNttHexl() = NTT(degree, q);
            initialized = true;
        }
    }

    void printHexlParams() {
        cout << "q: " << getNttHexl().GetModulus() <<
             ", N: " << getNttHexl().GetDegree() <<
             ", ROU: " << getNttHexl().GetMinimalRootOfUnity() << endl;
    }

    std::unordered_map<int32_t, NttPolynomial>& getNttRotMap() {
        static std::unordered_map<int32_t, NttPolynomial> nttRotMap;
        return nttRotMap;
    }

    NttPolynomial& getNttRoterPoly(int32_t r) {
        return getNttRotMap().find(r)->second;
    }

    void initNttRotMap(int32_t degree) {
        static bool initialized = false;
        if (!initialized) {
            for (int32_t i = 0; i < degree; i++) {
                NttPolynomial tmp{degree};
                TorusPolynomial tmpT{degree};
                tmpT.coeffs[i] = 1;
                applyNtt(tmp, tmpT);
                getNttRotMap().insert({i, std::move(tmp)});
            }
            initialized = true;
        }
    }

    void applyNtt(NttPolynomial &out, const TorusPolynomial &in) {
        auto N = in.N;
        auto q = getNttHexl().GetModulus();
        for (size_t i = 0; i < N; i++) {
            if (in.coeffs[i] >= 0) {
                out.coeffs[i] = Ntt64(in.coeffs[i]);
            } else {
                out.coeffs[i] = Ntt64(in.coeffs[i] + q);
            }
        }
        getNttHexl().ComputeForward(out.coeffs.data(), out.coeffs.data(), 1, 1);
    }

    void applyIntt(TorusPolynomial &out, const NttPolynomial &in) {
        auto N = in.N;
        auto q = getNttHexl().GetModulus();
        auto halfQ = (q + 1) >> 1;
        std::vector<uint64_t> tmp(N);
        getNttHexl().ComputeInverse(tmp.data(), in.coeffs.data(), 1, 1);
        for (int i = 0; i < N; i++) {
            if (tmp[i] >= halfQ) {
                tmp[i] -= q;
            }
            if (tmp[i] > TORUS_MAX) {
                out.coeffs[i] = Torus(tmp[i] - TORUS_Q);
            } else {
                out.coeffs[i] = Torus(tmp[i]);
            }
        }
    }

    void calModularInnerProductNtt(NttPolynomial& res, const vector<NttPolynomial>& in1, const vector<NttPolynomial>& in2) {
        for (auto i = 0; i < in1.size(); i++) {
            calModularInnerProductNtt(res, in1[i], in2[i]);
        }
    }

    void calModularInnerProductNtt(NttPolynomial &acc, const NttPolynomial &in1, const NttPolynomial &in2) {
        auto N = in2.N;
        auto q = getNttHexl().GetModulus();
        NttPolynomial tmp{N};
        NttPolynomial tmp1 = acc;
        EltwiseMultMod(tmp.coeffs.data(), in1.coeffs.data(), in2.coeffs.data(), N, q, 1);
        EltwiseAddMod(acc.coeffs.data(), tmp.coeffs.data(), tmp1.coeffs.data(), N, q);
    }

    void EltwiseSignedMultMod(uint64_t* result, const uint64_t* operand1,
                        const uint64_t* operand2, const int32_t sign, uint64_t n, uint64_t modulus,
                        uint64_t input_mod_factor) {
      HEXL_CHECK(result != nullptr, "Require result != nullptr");
      HEXL_CHECK(operand1 != nullptr, "Require operand1 != nullptr");
      HEXL_CHECK(operand2 != nullptr, "Require operand2 != nullptr");
      HEXL_CHECK(n != 0, "Require n != 0");
      HEXL_CHECK(modulus > 1, "Require modulus > 1");
      HEXL_CHECK(input_mod_factor * modulus < (1ULL << 63),
                 "Require input_mod_factor * modulus < (1ULL << 63)");
      HEXL_CHECK(
          input_mod_factor == 1 || input_mod_factor == 2 || input_mod_factor == 4,
          "Require input_mod_factor = 1, 2, or 4")
      HEXL_CHECK_BOUNDS(operand1, n, input_mod_factor * modulus,
                        "operand1 exceeds bound " << (input_mod_factor * modulus))
      HEXL_CHECK_BOUNDS(operand2, n, input_mod_factor * modulus,
                        "operand2 exceeds bound " << (input_mod_factor * modulus))

      HEXL_VLOG(3, "Calling EltwiseMultModNative");
      switch (input_mod_factor) {
        case 1:
          EltwiseSignedMultModNative<1>(result, operand1, operand2, sign, n, modulus);
          break;
        case 2:
          EltwiseSignedMultModNative<2>(result, operand1, operand2, sign, n, modulus);
          break;
        case 4:
          EltwiseSignedMultModNative<4>(result, operand1, operand2, sign, n, modulus);
          break;
      }
      return;
    }
}