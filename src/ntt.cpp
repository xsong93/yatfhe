//
// Created by Xintong Song on 2024/1/11.
//
#include "ntt.h"
#include "ntt_constants.h"

using namespace std;

// Function to perform Number Theoretic Transform (NTT)
void ntt(const TorusPolynomial& torusPolynomial, LagrangePolynomial& lagrangePolynomial) {
    const vector<Torus>& input = torusPolynomial.coeffs;
    vector<uint64_t>& output = lagrangePolynomial.coeffs;
    const int32_t N = lagrangePolynomial.N;

    for (int i = 0; i < N; i++) {
        uint64_t inputValue = input[i] < 0 ? input[i] + MODULUS : input[i];
        output[i] = modMul(inputValue, phi_normal_2[i]);
    }
    bitRevShuffle(output, N);
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
                uint64_t b = (wb == 1) ? output[j] : modMul(output[j], wb);
                output[i] = modAdd(a, b);
                output[j] = modSub(a, b);
            }
            wb = wb_normal_2[wbarr++];
        }
    }
}

void bitRevShuffle(std::vector<uint64_t>& x, int N) {
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
            uint64_t temp = x[j];
            x[j] = x[i];
            x[i] = temp;
        }
    }
}

uint64_t modAdd(uint64_t x, uint64_t y) {
    return ((MODULUS - x) > y) ? (x + y) : (x + y - MODULUS);
}

uint64_t modSub(uint64_t x, uint64_t y) {
    return (x >= y) ? (x - y) : (MODULUS - x + y);
}

uint64_t modMul(uint64_t x, uint64_t y) {
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
