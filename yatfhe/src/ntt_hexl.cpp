//
// Created by Xintong on 25-4-16.
//

#include "yatfhe/ntt_hexl.h"
#include "yatfhe/polynomial.h"
#include "yautil/tool.h"

using namespace intel::hexl;

namespace NttHexl {

    NTT &nttHexl() {
        static NTT nttHexl;
        return nttHexl;
    }

    void initNttHexl(uint64_t degree, uint64_t q) {
        static bool initialized = false;
        if (!initialized) {
            nttHexl() = NTT(degree, q);
            initialized = true;
        }
    }

    void printHexlParams() {
        cout << "q: " << nttHexl().GetModulus() <<
             ", N: " << nttHexl().GetDegree() <<
             ", ROU: " << nttHexl().GetMinimalRootOfUnity() << endl;
    }

    std::unordered_map<int32_t, LagrangePolynomial>& getNttRotMap() {
        static std::unordered_map<int32_t, LagrangePolynomial> nttRotMap;
        return nttRotMap;
    }

    LagrangePolynomial& getNttRoterPoly(int32_t r) {
        return getNttRotMap().find(r)->second;
    }

    void initNttRotMap(int32_t degree) {
        static bool initialized = false;
        if (!initialized) {
            for (int32_t i = 0; i < degree; i++) {
                LagrangePolynomial tmp{degree};
                TorusPolynomial tmpT{degree};
                tmpT.coeffs[i] = 1;
                applyNtt(tmp, tmpT);
                getNttRotMap().insert({i, std::move(tmp)});
            }
            initialized = true;
        }
    }

    void applyNtt(LagrangePolynomial &out, const TorusPolynomial &in) {
        auto N = in.N;
        auto q = nttHexl().GetModulus();
        for (size_t i = 0; i < N; i++) {
            if (in.coeffs[i] >= 0) {
                out.coeffs[i] = Ntt64(in.coeffs[i]);
            } else {
                out.coeffs[i] = Ntt64(in.coeffs[i] + q);
            }
        }
        nttHexl().ComputeForward(out.coeffs.data(), out.coeffs.data(), 1, 1);
    }

    void applyIntt(TorusPolynomial &out, const LagrangePolynomial &in) {
        auto N = in.N;
        auto q = nttHexl().GetModulus();
        auto halfQ = (q + 1) >> 1;
        std::vector<uint64_t> tmp(N);
        nttHexl().ComputeInverse(tmp.data(), in.coeffs.data(), 1, 1);
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

    void lagrangePolynomialRotate(LagrangePolynomial& res, const LagrangePolynomial& in, int r) {
        auto N = res.N;
        auto q = nttHexl().GetModulus();
        auto roter = NttHexl::getNttRoterPoly(r).coeffs.data();
        EltwiseMultMod(res.coeffs.data(), in.coeffs.data(), roter, N, q, 1);
    }

    void calModularInnerProductNtt(LagrangePolynomial& res, const vector<LagrangePolynomial>& in1, const vector<LagrangePolynomial>& in2) {
        for (auto i = 0; i < in1.size(); i++) {
            calModularInnerProductNtt(res, in1[i], in2[i]);
        }
    }

    void calModularInnerProductNtt(LagrangePolynomial &acc, const LagrangePolynomial &in1, const LagrangePolynomial &in2) {
        auto N = in2.N;
        auto q = nttHexl().GetModulus();
        LagrangePolynomial tmp{N};
        LagrangePolynomial tmp1 = acc;
        EltwiseMultMod(tmp.coeffs.data(), in1.coeffs.data(), in2.coeffs.data(), N, q, 1);
        EltwiseAddMod(acc.coeffs.data(), tmp.coeffs.data(), tmp1.coeffs.data(), N, q);
    }
}