//
// Created for anonymous review.
//
#include "yatfhe/ntt.h"
#include "yatfhe/torus.h"

using namespace std;

namespace NttNative {

    void bitRevShuffle(std::vector<NttType> &x) {
        int j = 0;
        int b = 0;
        int N = int(x.size());
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

}