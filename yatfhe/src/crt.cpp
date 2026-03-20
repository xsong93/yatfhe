//
// Created for anonymous review.
//
#include "yatfhe/crt.h"
#include "yatfhe/numeric.h"
#include "yatfhe/trlwe.h"
#include "yautil/tool.h"

int calApproxCRTError(const std::vector<int>& f_tilde, const std::vector<int>& coeffs) {
    int infNorm = 0;
    for (auto i = 0; i < f_tilde.size(); i++) {
        infNorm  = std::max(infNorm, std::abs(f_tilde[i] - coeffs[i]));
    }
    return infNorm;
}

void calGadgetVector(std::vector<long>& w, const int Qlow, const std::vector<int>& highModuli) {
    int Q = highModuli[0];
    for (auto i = 1; i < highModuli.size(); i++) {
        Q *= highModuli[i];
    }
    auto Q1 = Q / highModuli[0];
    auto Q2 = Q / highModuli[1];
    w[0] = Qlow * Q1 * (modInverse(Qlow * Q1, highModuli[0]));
    w[1] = Qlow * Q2 * (modInverse(Qlow * Q2, highModuli[1]));
//    printf("w1:%ld, w2:%ld\n", w[0], w[1]);
}

void decompCrtExact(std::vector<std::vector<int8_t>>& f, const std::vector<int>& coeffs, const YatfheParameters& param) {
    for (size_t d = 0; d < param.d; d++) {
        auto& qd = param.qd[d];
        for (size_t j = 0; j < coeffs.size(); j++) {
            f[d][j] = static_cast<int8_t>(longModP(coeffs[j], qd));
        }
    }
}

void decompCrtExactIO(std::vector<std::vector<int8_t>>& f, const std::vector<int>& coeffs, const YatfheParameters& param) {
    for (size_t j = 0; j < coeffs.size(); j++) {
        for (size_t d = 0; d < param.d; d++) {
            auto qd = param.qd[d];
            f[j][d] = static_cast<int8_t>(longModP(coeffs[j], qd));
        }
    }
}

void reconstructCrtExact(std::vector<int32_t>& f_tilde, const std::vector<std::vector<int8_t>>& f, const YatfheParameters& param) {
    for (int i = 0; i < f_tilde.size(); i++) {
        long acc = 0;
        for (size_t j = 0; j < param.d; j++) {
            acc += f[j][i] * param.z[j];
        }
        f_tilde[i] = static_cast<int32_t>(longModP(acc, param.qCRT));
    }
}

void reconstructCrtExactIO(std::vector<int32_t>& f_tilde, const std::vector<std::vector<int8_t>>& f, const YatfheParameters& param) {
    for (int i = 0; i < f_tilde.size(); i++) {
        long acc = 0;
        for (size_t j = 0; j < param.d; j++) {
            acc += f[i][j] * param.z[j];
        }
        f_tilde[i] = static_cast<int32_t>(longModP(acc, param.qCRT));
    }
}

void decompCrtApprox(std::vector<std::vector<int8_t>>& f, const std::vector<int>& coeffs, const int Qlow, const std::vector<int>& lowModuli, const std::vector<int>& highModuli) {
    for (int i = 0; i < coeffs.size(); ++i) {
        int fi = coeffs[i];

        // Calculate the low part contribution
        int lowSum = 0;
        for (int q_u : lowModuli) {
            long inv = modInverse(Qlow / q_u, q_u);
            lowSum += Qlow / q_u * static_cast<int32_t>(longModP(inv * fi, q_u));
        }
        // Calculate f for high moduli
        for (auto j = 0; j < f.size(); j++) {
            f[j][i] =  static_cast<int8_t >(intModP(fi - lowSum, highModuli[j]));
        }
    }
}

void reconstructCrtApproxVec(std::vector<int>& f_tilde, const std::vector<std::vector<int8_t>>& f, const std::vector<long>& w, long q) {
    for (int i = 0; i < f_tilde.size(); i++) {
        long acc = 0;
        for (size_t j = 0; j < w.size(); j++) {
            acc += f[j][i] * w[j];
        }
        f_tilde[i] = static_cast<int32_t>(longModP(acc, q));
    }
}

int32_t reconstructCrtApproxSingleEle(const std::vector<int8_t>& f, const YatfheParameters& param) {
    long acc = 0;
    for (size_t j = 0; j < param.dh; j++) {
        acc += f[j] * param.w[j];
    }
    return static_cast<int32_t>(longModP(acc, param.qCRT));
}

void decompTrlweApproxCrt(std::vector<Trlwe8>& out, const std::vector<Trlwe8>& in, const YatfheParameters& param) {
    int32_t lowSumsA[param.k][param.N];
    int32_t lowSumsB[param.N];
    auto dh = param.dh;
    for (size_t k = 0; k < param.k; k++) {
        for (size_t j = 0; j < param.N; j++) {
            lowSumsA[k][j] = 0;
            for (size_t u = 0; u < param.dl; u++) {
                auto aLo = in[u + dh].a[k].coeffs[j];
                lowSumsA[k][j] += param.qLowDivQl[u] * aLo;
            }
        }
    }
    for (size_t j = 0; j < param.N; j++) {
        lowSumsB[j] = 0;
        for (size_t u = 0; u < param.dl; u++) {
            auto bLo = in[u + dh].b.coeffs[j];
            lowSumsB[j] += param.qLowDivQl[u] * bLo;
        }
    }
    for (size_t i = 0; i < param.dh; i++) {
        auto qHi = param.qh[i];
        auto& outA = out[i].a;
        auto& inA = in[i].a;
        for (size_t k = 0; k < param.k; k++) {
            auto& coeffAOut = outA[k].coeffs;
            auto& coeffAHi = inA[k].coeffs;
            for (size_t j = 0; j < param.N; j++) {
                coeffAOut[j] = static_cast<int8_t>(intModP(static_cast<int32_t>(coeffAHi[j]) - lowSumsA[k][j], qHi));
            }
        }
        auto& coeffBHi = in[i].b.coeffs;
        auto& coeffBOut = out[i].b.coeffs;
        for (size_t j = 0; j < param.N; j++) {
            coeffBOut[j] = static_cast<int8_t>(intModP(static_cast<int32_t>(coeffBHi[j]) - lowSumsB[j], qHi));
        }
    }
}

// l -> l*d
void broadcastTrlweApproxCrt(std::vector<std::vector<Trlwe8>>& out, const std::vector<Trlwe8>& in, const YatfheParameters& param) {
    for (size_t d = 0; d < param.d; d++) {
        auto qd = param.qd[d];
        auto qdHalf = param.qdHalf[d];
        for (size_t dh = 0; dh < param.dh; dh++) {
            auto& outA = out[d][dh].a;
            auto& inA = in[dh].a;
            for (size_t k = 0; k < param.k; k++) {
                auto& coeffInA = inA[k].coeffs;
                auto& coeffOutA = outA[k].coeffs;
                for (size_t j = 0; j < param.N; j++) {
                    auto valA = coeffInA[j];
                    coeffOutA[j] = (valA > qdHalf || valA < -qdHalf) ? static_cast<int8_t>(intModP(valA, qd)) : valA;
                }
            }
            auto& inB = in[dh].b;
            auto& outB = out[d][dh].b;
            for (size_t j = 0; j < param.N; j++) {
                auto valB = inB.coeffs[j];
                outB.coeffs[j] = (valB > qdHalf || valB < -qdHalf) ? static_cast<int8_t>(intModP(valB, qd)) : valB;
            }
        }
    }
}