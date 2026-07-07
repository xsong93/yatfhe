//
// Created by Xintong Song on 2023/12/25.
//
#include <iostream>
#include "yatfhe/polynomial.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/numeric.h"
#include "yatfhe/ntt_hexl.h"
#include "yatfhe/gadget_decomposition.h"

/**
 * For a random rotator input, this method converts the rotator to a value within the range of polynomial length.
 * Besides, in order to correctly show the negacyclic property, use a bit indicator to keep track of the negative signs after rotation.
 * For example. N = 5
 * if a = 3, a0 + x^1 * a1 + x^2 * a2 + x^3* a3 + x^4 * a4 ---> x^3 * a0 + x^4 * a1 - a2 - x * a3 - x^2 * a4 -> (-a2, -a3, -a4, a0, a1)
 * if a = -1, a0 + x^1 * a1 + x^2 * a2 + x^3* a3 + x^4 * a4 ---> -x^4 * a0 + a1 + x * a2 + x^2 * a3 + x^3 * a4 -> (a1, a2, a3, a4, -a0) <=> a = 9 <=> -1 * (a = 4)
 * @param aTrue Minimized rotator a.
 * @param isWrap 1: no wrap. -1: wrap around.
 * @param a Original rotator.
 * @param N Polynomial length,
 */
void validateRotator(int& aTrue, int& isWrap, const int a, const int N) {
    aTrue = a % (2 * N);
    if (aTrue < 0) {
        aTrue += 2 * N;
    }
    isWrap = (aTrue < N) ? 1 : -1;
    aTrue = (aTrue < N) ? aTrue : aTrue - N;
}

void intPolyToDoublePoly(DoublePolynomial& output, const IntPolynomial& input) {
    for (auto i = 0; i < output.N; i++) {
        output.coeffs[i] = static_cast<double>(input.coeffs[i]);
    }
}

void torusPolyToDoublePoly(DoublePolynomial& output, const TorusPolynomial& input) {
    for (auto i = 0; i < output.N; i++) {
        output.coeffs[i] = torus32ToDouble(input.coeffs[i]);
    }
}

void doublePolyToTorusPoly(TorusPolynomial& output, const DoublePolynomial& input) {
    for (auto i = 0; i < output.N; i++) {
        output.coeffs[i] = doubleToTorus32(input.coeffs[i]);
    }
}

void torusPolyToIntPoly(IntPolynomial& output, const TorusPolynomial& input, const int mSize) {
    for (auto i = 0; i < input.N; i++) {
        output.coeffs[i] = modSwitchFromTorus32(input.coeffs[i], mSize);
    }
}

void intPolyToTorusPoly(TorusPolynomial& output, const IntPolynomial& input, const int mSize) {
    for (auto i = 0; i < input.N; i++) {
        output.coeffs[i] = modSwitchToTorus32(input.coeffs[i], mSize);
    }
}

void roundErrorTorusPoly(TorusPolynomial& target, const int torusBase) {
    for (auto i = 0 ; i < target.N; i++) {
        target.coeffs[i] = roundTorus32Error(target.coeffs[i], torusBase);
    }
}

void roundErrorDoublePoly(DoublePolynomial& target, const int torusBase) {
    for (auto i = 0 ; i < target.N; i++) {
        target.coeffs[i] = roundError(target.coeffs[i], torusBase);
    }
}

// vj = ((pj / q) mod p) / p, with every slot (including slot 0) centered on
// its sample point rather than starting there. Slot 0's left half wraps
// around to the end of the array, negated, since it represents the
// negacyclic reflection of the small-negative-index side of index 0.
void generateTestPolynomial(TorusPolynomial& v, const int modP, const int modQ) {
    const int N = v.N;
    const int boxSize = modQ / modP;
    const int halfBoxSize = boxSize / 2;

    for (auto i = 0; i < N; i++) {
        int tmp = intModP(i / boxSize, modP);
        v.coeffs[i] = modSwitchToTorus32(tmp, modP);
    }

    for (auto i = 0; i < halfBoxSize; i++) {
        v.coeffs[i] = -v.coeffs[i];
    }
    std::rotate(v.coeffs.begin(), v.coeffs.begin() + halfBoxSize, v.coeffs.begin() + N);
}

// 1: l
void generateTestPolynomialLt1(TorusPolynomial& v, const int t) {
    const int modP = MESSAGE_P;
    const int modQ = 2 * v.N;
    const int threshold = intModP(t, modP);
    const int boxSize = modQ / modP;
    const int halfBoxSize = boxSize / 2;

    for (auto i = 0; i < v.N; i++) {
        // Value in message space corresponding to coefficient index i.
        int msg = intModP(i / boxSize, modP);
        int tmp = msg < threshold ? 1 : 0;
        v.coeffs[i] = modSwitchToTorus32(tmp, MESSAGE_P);
    }
    for (auto i = 0; i < halfBoxSize; i++) {
        v.coeffs[i] = -v.coeffs[i];
    }
    std::rotate(v.coeffs.begin(), v.coeffs.begin() + halfBoxSize, v.coeffs.begin() + v.N);
}

// 0: leq
void generateTestPolynomialCompLeq0(TorusPolynomial& v, const int t) {
    const int modP = MESSAGE_P;
    const int modQ = 2 * v.N;
    const int threshold = intModP(t, modP);
    const int boxSize = modQ / modP;
    const int halfBoxSize = boxSize / 2;

    for (auto i = 0; i < v.N; i++) {
        // Value in message space corresponding to coefficient index i.
        int msg = intModP(i / boxSize, modP);
        int tmp = msg <= threshold ? 0 : 1;
        v.coeffs[i] = modSwitchToTorus32(tmp, MESSAGE_P);
    }
    for (auto i = 0; i < halfBoxSize; i++) {
        v.coeffs[i] = -v.coeffs[i];
    }
    std::rotate(v.coeffs.begin(), v.coeffs.begin() + halfBoxSize, v.coeffs.begin() + v.N);
}

void generateTestPolynomialCompWithValue(TorusPolynomial& tv, const int t, const Integer v) {
    const int modP = MESSAGE_P;
    const int modQ = 2 * tv.N;
    const int threshold = intModP(t, modP);
    const int boxSize = modQ / modP;
    const int halfBoxSize = boxSize / 2;

    for (auto i = 0; i < tv.N; i++) {
        auto msg = intModP(i / boxSize, modP);
        auto tmp = msg == threshold ? v : 0;
        tv.coeffs[i] = modSwitchToTorus32(tmp, MESSAGE_P);
    }
    for (auto i = 0; i < halfBoxSize; i++) {
        tv.coeffs[i] = -tv.coeffs[i];
    }
    std::rotate(tv.coeffs.begin(), tv.coeffs.begin() + halfBoxSize, tv.coeffs.begin() + tv.N);
}

void generateTestPolynomialOne(TorusPolynomial& v) {
    const int modP = MESSAGE_P;
    v.coeffs[0] = modSwitchToTorus32(1, modP);
}

void generateTestPolynomialValue(TorusPolynomial& tv, const Integer v) {
    const int modP = MESSAGE_P;
    const int modQ = 2 * tv.N;
    const int boundary = modQ / modP / 2;

    for (auto i = 0; i < tv.N; i++) {
        const auto tmp = i < boundary ? v : i >= tv.N - boundary ? -v : 0;
        tv.coeffs[i] = modSwitchToTorus32(tmp, MESSAGE_P);
    }
}

// output = (X^{a}) * input
void rotateTorusPolynomial(TorusPolynomial& out, const int a, const TorusPolynomial& input) {
    const auto N = input.N;
    int aTrue, isWrap;
    validateRotator(aTrue, isWrap, a, N);
    Torus tmp = 0;
    for (auto i = 0; i < N; i++) {
        tmp = (i < aTrue) ? (-input.coeffs[i - aTrue + N] * isWrap) : (input.coeffs[i - aTrue] * isWrap);
        out.coeffs[i] = static_cast<Torus>(longModP(tmp, TORUS_Q));
    }
}

// Fused: accum += X^a * input, no temporary buffer needed.
// Split into two contiguous loops to enable auto-vectorization (no branch on i).
void rotateAccumulateTorusPolynomial(TorusPolynomial& accum, const int aTrue, const int isWrap,
                                      const TorusPolynomial& input) {
    const int N = input.N;
    // i in [0, aTrue): rotated[i] = -input[i - aTrue + N] * isWrap
    for (int i = 0; i < aTrue; i++) {
        const int64_t rotated = -static_cast<int64_t>(input.coeffs[i - aTrue + N]) * isWrap;
        accum.coeffs[i] = static_cast<Torus>(longModP(accum.coeffs[i] + rotated, TORUS_Q));
    }
    // i in [aTrue, N): rotated[i] = input[i - aTrue] * isWrap
    for (int i = aTrue; i < N; i++) {
        const int64_t rotated = static_cast<int64_t>(input.coeffs[i - aTrue]) * isWrap;
        accum.coeffs[i] = static_cast<Torus>(longModP(accum.coeffs[i] + rotated, TORUS_Q));
    }
}

// output = (X^{a} - 1) * input = x^a * input - input
void rotateTorusPolynomialMinusOne(TorusPolynomial& out, const int a, const TorusPolynomial& input) {
    const auto N = input.N;
    int aTrue, isWrap;
    validateRotator(aTrue, isWrap, a, N);
    Torus tmp = 0;
    for (auto i = 0; i < N; i++) {
        tmp = (i < aTrue) ? (-input.coeffs[i - aTrue + N] * isWrap) : (input.coeffs[i - aTrue] * isWrap);
        out.coeffs[i] = subTorus(TORUS_Q, tmp, input.coeffs[i]);
    }
}

void rotateIntPolynomial(IntPolynomial& out, const int a, const IntPolynomial& input, const int64_t p) {
    const auto N = input.N;
    int aTrue, isWrap;
    validateRotator(aTrue, isWrap, a, N);
    Integer tmp;
    for (auto i = 0; i < N; i++) {
        tmp = (i < aTrue) ? (-input.coeffs[i - aTrue + N] * isWrap) : (input.coeffs[i - aTrue] * isWrap);
        out.coeffs[i] = static_cast<Integer>(longModP(tmp, p));
    }
}

void rotateInt8Polynomial(Int8Polynomial& out, const int a, const Int8Polynomial& input, int modP) {
    const auto N = input.N;
    int aTrue, isWrap;
    validateRotator(aTrue, isWrap, a, N);
    int tmp = 0;
    for (size_t i = 0; i < N; i++) {
        tmp = ((i < aTrue) ? (-input.coeffs[i - aTrue + N]) : (input.coeffs[i - aTrue])) * isWrap;
        out.coeffs[i] = static_cast<int8_t>(intModP(tmp, modP));
    }
}

void rotateInt8PolynomialMinusOne(Int8Polynomial& out, const int a, const Int8Polynomial& input, int modP) {
    const auto N = input.N;
    int aTrue, isWrap;
    validateRotator(aTrue, isWrap, a, N);
    int tmp = 0;
    for (size_t i = 0; i < N; i++) {
        tmp = ((i < aTrue) ? (-input.coeffs[i - aTrue + N]) : (input.coeffs[i - aTrue])) * isWrap;
        out.coeffs[i] = static_cast<int8_t>(intModP(tmp - input.coeffs[i], modP));
    }
}

void multTorusPolynomial(TorusPolynomial& res, const TorusPolynomial& poly1, const TorusPolynomial& poly2) {
    const int N = res.N;
    for (auto i = 0; i < N; i++) {
        int64_t tmp = 0;
        for (auto j = 0; j < N; j++) {
            if (j <= i) {
                tmp += multTorus(TORUS_Q, poly1.coeffs[j], poly2.coeffs[i - j]);
            } else {
                tmp -= multTorus(TORUS_Q, poly1.coeffs[j], poly2.coeffs[N + i - j]);
            }
        }
        res.coeffs[i] = static_cast<Torus>(longModP(tmp, TORUS_Q));
    }
}

void multIntPolynomialModQ(IntPolynomial& res, const IntPolynomial& poly1, const IntPolynomial& poly2, const int64_t q) {
    const int N = res.N;
    for (auto i = 0; i < N; i++) {
        int64_t tmp = 0;
        for (auto j = 0; j < N; j++) {
            if (j <= i) {
                tmp += modMulQ(poly1.coeffs[j], poly2.coeffs[i - j], q);
            } else {
                tmp -= modMulQ(poly1.coeffs[j], poly2.coeffs[N + i - j], q);
            }
        }
        res.coeffs[i] = static_cast<Torus>(longModP(tmp, q));
    }
}

void multInt8Polynomial(Int8Polynomial& res, const Int8Polynomial& poly1, const Int8Polynomial& poly2, const int q) {
    const int N = res.N;
    for (auto i = 0; i < N; i++) {
        int tmp = 0;
        for (auto j = 0; j < N; j++) {
            if (j <= i) {
                tmp += modMulQ(poly1.coeffs[j], poly2.coeffs[i - j], q);
            } else {
                tmp -= modMulQ(poly1.coeffs[j], poly2.coeffs[N + i - j], q);
            }
        }
        res.coeffs[i] = static_cast<int8_t>(intModP(tmp, q));
    }
}

void multInt8PolynomialAcc(Int8Polynomial& res, const Int8Polynomial& poly1, const Int8Polynomial& poly2, const int q) {
    const int N = res.N;
    for (auto i = 0; i < N; i++) {
        int tmp = 0;
        for (auto j = 0; j < N; j++) {
            if (j <= i) {
                tmp += modMulQ(poly1.coeffs[j], poly2.coeffs[i - j], q);
            } else {
                tmp -= modMulQ(poly1.coeffs[j], poly2.coeffs[N + i - j], q);
            }
        }
        res.coeffs[i] = static_cast<int8_t>(intModP(tmp + res.coeffs[i], q));
    }
}

void multIntPolynomialAcc(IntPolynomial& res, const IntPolynomial& poly1, const IntPolynomial& poly2) {
    const int N = res.N;
    for (auto i = 0; i < N; i++) {
        Torus tmp = 0;
        for (auto j = 0; j < N; j++) {
            tmp = (j <= i) ? (tmp + poly1.coeffs[j] * poly2.coeffs[i - j]) : (tmp - poly1.coeffs[j] * poly2.coeffs[N + i - j]);
        }
        res.coeffs[i] += tmp;
    }
}

void multTorusPolynomialAcc(TorusPolynomial& res, const TorusPolynomial& poly1, const TorusPolynomial& poly2) {
    const int N = res.N;
    for (auto i = 0; i < N; i++) {
        int64_t tmp = 0;
        for (auto j = 0; j < N; j++) {
            if (j <= i) {
                tmp += multTorus(TORUS_Q, poly1.coeffs[j], poly2.coeffs[i - j]);
            } else {
                tmp -= multTorus(TORUS_Q, poly1.coeffs[j], poly2.coeffs[N + i - j]);
            }
        }
        res.coeffs[i] = addTorus(TORUS_Q, res.coeffs[i], static_cast<Torus>(longModP(tmp, TORUS_Q)));
    }
}

// res = poly1 + poly2
void addIntPolynomial(IntPolynomial& res, const IntPolynomial& poly1, const IntPolynomial& poly2) {
    const int N = res.N;
    for (int i = 0; i < N; i++) {
        res.coeffs[i] = poly1.coeffs[i] + poly2.coeffs[i];
    }
}

// void addTorusPolynomial(TorusPolynomial& res, const TorusPolynomial& poly1, const TorusPolynomial& poly2) {
//     const int N = res.N;
//     for (int i = 0; i < N; i++) {
//         res.coeffs[i] = addTorus(poly1.coeffs[i], poly2.coeffs[i]);
//     }
// }

/**
 * Add or sub a value to every coefficients of the target polynomial.
 * @param poly Target polynomial.
 * @param offset Offset value.
 * @param isAdd True: add offset. Otherwise, subtract offset.
 */
void addSubIntPolynomialWithOffset(IntPolynomial& poly, const int offset, const bool isAdd) {
    for (auto i = 0; i < poly.N; i++) {
        poly.coeffs[i] = isAdd ? (poly.coeffs[i] + offset) : (poly.coeffs[i] - offset);
    }
}

// res = poly1 - poly2
void subIntPolynomial(IntPolynomial& res, const IntPolynomial& poly1, const IntPolynomial& poly2) {
    const int N = res.N;
    for (int i = 0; i < N; i++) {
        res.coeffs[i] = poly1.coeffs[i] - poly2.coeffs[i];
    }
}

void subTorusPolynomial(TorusPolynomial& res, const TorusPolynomial& poly1, const TorusPolynomial& poly2) {
    const int N = res.N;
    for (int i = 0; i < N; i++) {
        res.coeffs[i] = subTorus(TORUS_Q, poly1.coeffs[i], poly2.coeffs[i]);
    }
}

void rotateNttPolynomial(NttPolynomial& res, const NttPolynomial& in, int r) {
    auto N = res.N;
    int rTrue, isWrap;
    validateRotator(rTrue, isWrap, r, N);
    auto q = NttHexl::getNttHexl().GetModulus();
    auto roter = NttHexl::getNttRoterPoly(rTrue, isWrap).coeffs.data();
    EltwiseMultMod(res.coeffs.data(), in.coeffs.data(), roter, N, q, 1);
}

void rotateNttPolynomialMinusOne(NttPolynomial& res, const NttPolynomial& in, int r) {
    auto N = res.N;
    int rTrue, isWrap;
    validateRotator(rTrue, isWrap, r, N);
    auto q = NttHexl::getNttHexl().GetModulus();
    auto roter = NttHexl::getNttRoterPolyMinusOne(rTrue, isWrap).coeffs.data();
    EltwiseMultMod(res.coeffs.data(), in.coeffs.data(), roter, N, q, 1);
}

void genNttPolynomialWithValueAt(NttPolynomial& lagrangePolynomial, const int value, const int position) {
    TorusPolynomial tmp{lagrangePolynomial.N};
    tmp.coeffs[position] = value;
    NttHexl::applyNtt(lagrangePolynomial, tmp);
}

// accum += poly
void accumulateNttPolynomial(NttPolynomial& accum, NttPolynomial& poly) {
    const auto N = accum.N;
    const auto q = NttHexl::getNttHexl().GetModulus();
    for (auto i = 0; i < N; i++) {
        accum.coeffs[i] = AddUIntMod(accum.coeffs[i], poly.coeffs[i], q);
    }
}

// void addNttPolynomial(NttPolynomial& output, const NttPolynomial& input1, const NttPolynomial& input2) {
//     const auto q = NttHexl::getNttHexl().GetModulus();
//     for (auto i = 0; i < input1.N; i++) {
//         output.coeffs[i] = AddUIntMod(input1.coeffs[i], input2.coeffs[i], q);
//     }
// }

void subNttPolynomial(NttPolynomial& output, const NttPolynomial& input1, const NttPolynomial& input2) {
    const auto q = NttHexl::getNttHexl().GetModulus();
    for (auto i = 0; i < input1.N; i++) {
        output.coeffs[i] = SubUIntMod(input1.coeffs[i], input2.coeffs[i], q);
    }
}

void inverseGadgetDecomposePolynomial(vector<TorusPolynomial>& output, const IntPolynomial& input, const YatfheParameters& param) {
    const auto N = input.N;
    const int l = output.size();

    for (auto lvl = 0; lvl < l; lvl++) {
        for (auto j = 0; j < N; j++) {
            output[lvl].coeffs[j] = (input.coeffs[j]) << (param.torusBits - (lvl + 1) * param.radixBits);
        }
    }
}