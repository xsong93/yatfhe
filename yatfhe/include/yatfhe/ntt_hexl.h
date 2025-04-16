//
// Created by Xintong on 25-4-16.
//

#ifndef HLS_YATFHE_NTT_HEXL_H
#define HLS_YATFHE_NTT_HEXL_H

#include <hexl/hexl.hpp>
#include "yatfhe/polynomial.h"


using namespace intel::hexl;

namespace NttHexl {

    NTT &nttHexl();

    void initNttHexl(uint64_t degree, uint64_t q);

    void printHexlParams();

    void applyNtt(LagrangePolynomial &out, const TorusPolynomial &in);

    void applyIntt(TorusPolynomial &out, const LagrangePolynomial &in);

    void calModularInnerProductNtt(LagrangePolynomial& res, const vector<LagrangePolynomial>& in1, const vector<LagrangePolynomial>& in2);

    void calModularInnerProductNtt(LagrangePolynomial &acc, const LagrangePolynomial &in1, const LagrangePolynomial &in2);

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
