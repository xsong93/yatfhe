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

    void rotateNttPolynomial(NttPolynomial& res, const NttPolynomial& in, int r) {
        auto N = res.N;
        auto q = getNttHexl().GetModulus();
        auto roter = NttHexl::getNttRoterPoly(r).coeffs.data();
        EltwiseMultMod(res.coeffs.data(), in.coeffs.data(), roter, N, q, 1);
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
}