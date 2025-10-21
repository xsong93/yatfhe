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

    NttPolynomial &getNttRoterPoly(int32_t rTrue, int32_t isWrap);

    NttPolynomial &getNttRoterPolyMinusOne(int32_t rTrue, int32_t isWrap);

    NttPolynomial &getNttGadgetRecomper(const int32_t currL);

    void initNttRotMap(int32_t degree);

    void initNttRotMinusOneMap(int32_t degree);

    void initNttGadgetRecompMap(int32_t bitLength, int32_t radixBit, int32_t level, int32_t degree);

    void applyNtt(NttPolynomial &out, const TorusPolynomial &in);

    void applyIntt(TorusPolynomial &out, const NttPolynomial &in);

    void calModularInnerProductNtt(NttPolynomial& res, const vector<NttPolynomial>& in1, const vector<NttPolynomial>& in2);

    void calModularInnerProductNtt(NttPolynomial &acc, const NttPolynomial &in1, const NttPolynomial &in2);

    template<typename... NttPolyArgs>
    void addNttPolynomial(NttPolynomial& output, const NttPolyArgs&... inputs) {
        const auto q = getNttHexl().GetModulus();
        const int N = output.N;
        for (int i = 0; i < N; i++) {
            uint64_t tmp = (inputs.coeffs[i] + ...);
            output.coeffs[i] = tmp % q;
        }
    }

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

}

#endif //HLS_YATFHE_NTT_HEXL_H
