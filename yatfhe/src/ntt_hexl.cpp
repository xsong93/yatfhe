//
// Created for anonymous review.
//

#include <unordered_map>
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

    std::unordered_map<int32_t, NttPolynomial> &getNttRotMap() {
        static std::unordered_map<int32_t, NttPolynomial> nttRotMap;
        return nttRotMap;
    }

    std::unordered_map<int32_t, NttPolynomial> &getNttRotMinusOneMap() {
        static std::unordered_map<int32_t, NttPolynomial> nttRotMapMinusOneMap;
        return nttRotMapMinusOneMap;
    }

    std::unordered_map<int32_t, NttPolynomial> &getNttRotInverseMap() {
        static std::unordered_map<int32_t, NttPolynomial> nttRotInverseMap;
        return nttRotInverseMap;
    }

    std::unordered_map<int32_t, NttPolynomial> &getNttRotMinusOneInverseMap() {
        static std::unordered_map<int32_t, NttPolynomial> nttRotMinusOneInverseMap;
        return nttRotMinusOneInverseMap;
    }

    std::unordered_map<int32_t, NttPolynomial> &getNttGadgetRecompMap() {
        static std::unordered_map<int32_t, NttPolynomial> nttGadgetRecompMap;
        return nttGadgetRecompMap;
    }

    NttPolynomial &getNttRoterPoly(const int32_t rTrue, const int32_t isWrap) {
        return isWrap == 1 ? getNttRotMap().find(rTrue)->second : getNttRotInverseMap().find(rTrue)->second;
    }

    NttPolynomial &getNttRoterPolyMinusOne(const int32_t rTrue, const int32_t isWrap) {
        return isWrap == 1 ? getNttRotMinusOneMap().find(rTrue)->second : getNttRotMinusOneInverseMap().find(rTrue)->second;
    }

    NttPolynomial &getNttGadgetRecomper(const int32_t currL) {
        return getNttGadgetRecompMap().find(currL)->second;
    }

    void initNttRotMap(int32_t degree) {
        static bool NttRotInitialized = false;
        if (!NttRotInitialized) {
            for (int32_t i = 0; i < degree; i++) {
                NttPolynomial tmp{degree};
                TorusPolynomial tmpT{degree};
                tmpT.coeffs[i] = 1;
                applyNtt(tmp, tmpT);
                getNttRotMap().insert({i, tmp});
                tmpT.coeffs[i] = -1;
                applyNtt(tmp, tmpT);
                getNttRotInverseMap().insert({i, tmp});
            }
            NttRotInitialized = true;
        }
    }

    void initNttRotMinusOneMap(int32_t degree) {
        static bool NttRotMinusOneInitialized = false;
        if (!NttRotMinusOneInitialized) {
            for (int32_t i = 0; i < degree; i++) {
                NttPolynomial tmp{degree};
                TorusPolynomial tmpT{degree};
                tmpT.coeffs[i] = 1;
                tmpT.coeffs[0] -= 1;
                applyNtt(tmp, tmpT);
                getNttRotMinusOneMap().insert({i, tmp});
                tmpT = TorusPolynomial{degree};
                tmpT.coeffs[i] = -1;
                tmpT.coeffs[0] -= 1;
                applyNtt(tmp, tmpT);
                getNttRotMinusOneInverseMap().insert({i, tmp});
            }
            NttRotMinusOneInitialized = true;
        }
    }

    void initNttGadgetRecompMap(int32_t bitLength, int32_t radixBit, int32_t level, int32_t degree) {
        for (auto l = 0; l < level; l++) {
            TorusPolynomial tmp{degree};
            NttPolynomial ntt{degree};
            tmp.coeffs[0] = 1 << (bitLength - (l + 1) * radixBit);
            applyNtt(ntt, tmp);
            getNttGadgetRecompMap().insert({l, ntt});
        }
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
                out.coeffs[i] = Torus(longModP(tmp[i], TORUS_Q));
                continue;
            }
            out.coeffs[i] = Torus(tmp[i]);
        }
    }

    void calModularInnerProductNtt(NttPolynomial &res, const vector<NttPolynomial> &in1, const vector<NttPolynomial> &in2) {
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