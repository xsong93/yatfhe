//
// Created by Xintong on 25-4-16.
//

#include <unordered_map>
#include <stdexcept>
#include <string>
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
        static uint64_t curDegree = 0, curQ = 0;
        if (degree != curDegree || q != curQ) {
            getNttHexl() = NTT(degree, q);
            curDegree = degree;
            curQ = q;
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

    namespace {
        bool rotTablesStale(int32_t degree, int32_t& curDegree, uint64_t& curQ) {
            const uint64_t q = getNttHexl().GetModulus();
            if (degree == curDegree && q == curQ) return false;
            curDegree = degree;
            curQ = q;
            return true;
        }
    }

    void initNttRotMap(int32_t degree) {
        static int32_t curDegree = -1;
        static uint64_t curQ = 0;
        if (!rotTablesStale(degree, curDegree, curQ)) return;
        getNttRotMap().clear();
        getNttRotInverseMap().clear();
        for (int32_t i = 0; i < degree; i++) {
            NttPolynomial tmp{degree};
            TorusPolynomial tmpT{degree};
            tmpT.coeffs[i] = 1;
            applyNtt(tmp, tmpT);
            getNttRotMap().insert_or_assign(i, tmp);
            tmpT.coeffs[i] = -1;
            applyNtt(tmp, tmpT);
            getNttRotInverseMap().insert_or_assign(i, tmp);
        }
    }

    void initNttRotMinusOneMap(int32_t degree) {
        static int32_t curDegree = -1;
        static uint64_t curQ = 0;
        if (!rotTablesStale(degree, curDegree, curQ)) return;
        getNttRotMinusOneMap().clear();
        getNttRotMinusOneInverseMap().clear();
        for (int32_t i = 0; i < degree; i++) {
            NttPolynomial tmp{degree};
            TorusPolynomial tmpT{degree};
            tmpT.coeffs[i] = 1;
            tmpT.coeffs[0] -= 1;
            applyNtt(tmp, tmpT);
            getNttRotMinusOneMap().insert_or_assign(i, tmp);
            tmpT = TorusPolynomial{degree};
            tmpT.coeffs[i] = -1;
            tmpT.coeffs[0] -= 1;
            applyNtt(tmp, tmpT);
            getNttRotMinusOneInverseMap().insert_or_assign(i, tmp);
        }
    }

    void initNttGadgetRecompMap(int32_t bitLength, int32_t radixBit, int32_t level, int32_t degree) {
        static int32_t curBitLength = -1, curRadixBit = -1, curDegree = -1;
        static uint64_t curQ = 0;
        const uint64_t q = getNttHexl().GetModulus();
        auto& map = getNttGadgetRecompMap();
        if (bitLength != curBitLength || radixBit != curRadixBit
                || degree != curDegree || q != curQ) {
            map.clear();
            curBitLength = bitLength;
            curRadixBit = radixBit;
            curDegree = degree;
            curQ = q;
        }
        for (auto l = 0; l < level; l++) {
            const int shift = bitLength - (l + 1) * radixBit;
            if (shift < 0) {
                // 1 << negative is undefined, and a gadget wider than the torus has no
                // meaning.
                throw std::invalid_argument(
                    "initNttGadgetRecompMap: level " + std::to_string(l + 1) +
                    " of radix 2^" + std::to_string(radixBit) +
                    " exceeds the torus width " + std::to_string(bitLength) +
                    "; require l * radixBits <= torusBits");
            }
            TorusPolynomial tmp{degree};
            NttPolynomial ntt{degree};
            tmp.coeffs[0] = static_cast<Torus>(1) << shift;
            applyNtt(ntt, tmp);
            map.insert_or_assign(l, ntt);
        }
    }

    void applyIntt(TorusPolynomial &out, const NttPolynomial &in) {
        auto N = in.N;
        auto q = getNttHexl().GetModulus();
        auto halfQ = (q + 1) >> 1;
        // Per-thread scratch: this runs 2x per external product.
        thread_local std::vector<uint64_t> tmp;
        if (tmp.size() < static_cast<size_t>(N)) tmp.resize(N);
        getNttHexl().ComputeInverse(tmp.data(), in.coeffs.data(), 1, 1);
        // Branch outside the loop: the power-of-two path vectorises, the general
        // path stays correct for a non-power-of-two TORUS_Q (Q_CRT).
        if (TORUS_IS_POW2) {
            const int shift = TORUS_SHIFT;
            for (int i = 0; i < N; i++) {
                // Centre mod qNtt, then reduce to the torus.
                const uint64_t centred = tmp[i] - (q & static_cast<uint64_t>(-static_cast<int64_t>(tmp[i] >= halfQ)));
                out.coeffs[i] = static_cast<Torus>(longModPow2(static_cast<int64_t>(centred), shift));
            }
        } else {
            const int64_t torusQ = TORUS_Q;
            for (int i = 0; i < N; i++) {
                const uint64_t centred = tmp[i] - (q & static_cast<uint64_t>(-static_cast<int64_t>(tmp[i] >= halfQ)));
                out.coeffs[i] = static_cast<Torus>(longModP(static_cast<int64_t>(centred), torusQ));
            }
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
        // Per-thread scratch, called 4x per level.
        thread_local std::vector<NttType> tmp;
        if (tmp.size() < static_cast<size_t>(N)) tmp.resize(N);
        EltwiseMultMod(tmp.data(), in1.coeffs.data(), in2.coeffs.data(), N, q, 1);
        // EltwiseAddMod supports output == input1, so acc can be updated in-place
        EltwiseAddMod(acc.coeffs.data(), acc.coeffs.data(), tmp.data(), N, q);
    }
}