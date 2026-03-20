//
// Created by Xintong Song on 2024/1/11.
//

#include "yatfhe/ntt_old.h"
#include "yautil/ntt_constants.h"
#include "yautil/tool.h"
using namespace std;

// Function to perform Number Theoretic Transform (NTT)
void applyNttOld(NttPolynomial& out, const IntPolynomial& in) {
    auto& input = in.coeffs;
    auto& output = out.coeffs;
    const auto N = out.N;

    for (int i = 0; i < N; i++) {
        uint64_t inputValue = input[i] < 0 ? input[i] + MODULUS : input[i];
        output[i] = modMulOld(inputValue, phi_normal_2[i]);
    }
    bitRevShuffleOld(output, N);
    int32_t wbarr = 0;

    // Loop for the NTT algorithm
    for (int transSize = 2; transSize <= N; transSize *= 2) {
        uint64_t wb = 1;
        for (int t = 0; t < (transSize >> 1); t++) {
            for (int trans = 0; trans < (N / transSize); trans++) {
                int i = trans * transSize + t;
                int j = i + (transSize >> 1);

                // Perform butterfly operations
                uint64_t a = output[i];
                uint64_t b = (wb == 1) ? output[j] : modMulOld(output[j], wb);
                output[i] = modAddOld(a, b);
                output[j] = modSubOld(a, b);
            }
            wb = wb_normal_2[wbarr++];
        }
    }
}

void applyNttTorusOld(NttPolynomial& out, const TorusPolynomial & in, const int mSize) {
    IntPolynomial intPolynomial(in.N);
    torusPolyToIntPoly(intPolynomial, in, mSize);
    printArray(intPolynomial.coeffs, "intPolynomial@applyNttTorus");
    applyNttOld(out, intPolynomial);
}

void applyInttOld(IntPolynomial& out, const NttPolynomial& in) {
    vector<NttType> input(in.coeffs.size());
    copy(in.coeffs.begin(), in.coeffs.end(), input.begin());
    int32_t N = in.N;
    NttPolynomial temp {N};
    vector<NttType>& tmp = temp.coeffs;
    vector<Integer>& output = out.coeffs;
    int inv = 0;
    bitRevShuffleOld(input, N);
    for (int transSize = 2; transSize <= N; transSize = transSize * 2) {
        uint64_t wb = 1;
        for (int t = 0; t < (transSize >> 1); t++) {
            for (int trans = 0; trans < (N / transSize); trans++) {
                int i = trans * transSize + t;
                int j = i + (transSize >> 1);
                uint64_t a = input[i];
                uint64_t b = (wb == 1) ? input[j] : modMulOld(input[j], wb);
                input[i] = modAddOld(a, b);
                input[j] = modSubOld(a, b);
            }
            wb = wb_inverse_2[inv++];
        }
    }
    for (int i = 0; i < N; i++) {
        tmp[i] = modMulOld(input[i], scale_2);    //scale_2*phi inversev, modulus  (phi inverse sclaed)
        tmp[i] = modMulOld(tmp[i], phi_inverse_2[i]);
    }

    uint64_t med = MODULUS / 2;
    for (int i = 0; i < N; i++) {
        output[i] = (int32_t)((tmp[i] & 0xffffffff) - (tmp[i] > med));
    }
}

void bitRevShuffleOld(std::vector<NttType>& x, int N) {
    int j = 0;
    int b = 0;

    for (int i = 1; i < N; i++) {
        b = N >> 1;  // Initialize b to half of N
        while (j >= b) {
            j -= b;  // Perform bit-reversal
            b >>= 1;
        }
        j += b;  // Move to the next position

        // Swap elements if the bit-reversed index is greater than the current index
        if (j > i) {
            NttType temp = x[j];
            x[j] = x[i];
            x[i] = temp;
        }
    }
}

uint64_t modAddOld(uint64_t x, uint64_t y) {
    return ((MODULUS - x) > y) ? (x + y) : (x + y - MODULUS);
}

uint64_t modSubOld(uint64_t x, uint64_t y) {
    return (x >= y) ? (x - y) : (MODULUS - y + x);
}

uint64_t modMulOld(uint64_t x, uint64_t y) {
    // Break down x and y into 32-bit components
    auto x0 = (uint32_t)x;
    auto x1 = (uint32_t)(x >> 32);
    auto y0 = (uint32_t)y;
    auto y1 = (uint32_t)(y >> 32);

    // Perform 64-bit multiplication
    uint64_t x0y0 = (uint64_t)x0 * (uint64_t)y0;
    uint64_t x0y1 = (uint64_t)x0 * (uint64_t)y1;
    uint64_t x1y0 = (uint64_t)x1 * (uint64_t)y0;
    uint64_t x1y1 = (uint64_t)x1 * (uint64_t)y1;

    // Compute partial products and handle carry
    auto d = (uint32_t)x0y0;
    uint64_t pp1 = (x0y0 >> 32) + (uint32_t)(x1y0) + (uint32_t)(x0y1);
    auto c = (uint32_t)pp1;
    uint64_t pp2 = (x1y0 >> 32) + (x0y1 >> 32) + (uint32_t)(x1y1);
    uint64_t pp3 = (pp1 >> 32) + (uint32_t)(pp2);

    // Handle overflow and underflow
    uint32_t a = (pp2 >> 32) + (x1y1 >> 32);
    uint64_t bpc = (uint32_t)pp3 + (uint64_t)c;
    bpc = ((bpc + (bpc >> 32)) << 32) - (bpc >> 32);
    uint64_t minus = ((uint64_t)a + ((uint64_t)(uint32_t)pp3));
    uint64_t plus = bpc + (uint64_t)d;

    // Return the result modulo MODULUS
    if (plus >= minus) {
        return (plus - minus);
    }
    return MODULUS - minus + plus;
}

// output_j = aj * bj mod p
void modularMultOld(std::vector<uint64_t>& output, const std::vector<uint64_t>& coeffsA, const std::vector<uint64_t>& coeffsB) {
    const auto N = output.size();
    for (auto j = 0; j < N; j++) {
        output[j] = modMulOld(coeffsA[j], coeffsB[j]);
    }
}

// b += a * s mod p
void modularAccumulateOld(std::vector<uint64_t>& coeffsB, const std::vector<uint64_t>& coeffsA, const std::vector<uint64_t>& coeffsS) {
    const auto N = coeffsB.size();
    for (auto j = 0; j < N; j++) {
        auto tmp = modMulOld(coeffsA[j], coeffsS[j]);
        coeffsB[j] = modAddOld(coeffsB[j], tmp);
    }
}

// b = aN * sN
void calModularInnerProductNttOld(NttPolynomial& b, const vector<NttPolynomial>& a, const vector<NttPolynomial>& s) {
    for (auto i = 0; i < a.size(); i++) {
        modularAccumulateOld(b.coeffs, a[i].coeffs, s[i].coeffs);
    }
}