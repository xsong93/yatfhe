//
// Created by Xintong Song on 2024/1/10.
//

#ifndef HLS_YATFHE_NTT_H
#define HLS_YATFHE_NTT_H

#include <vector>
#include "yatfhe/polynomial.h"
#include "yatfhe/numeric.h"

namespace NttNative {

    const string STR_NTT = "NWC-DIT-NR-NTT";
    const string STR_INTT = "NWC-DIF-RN-INTT";

    template<typename T>
    void bitRev(std::vector<T> &x) {
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
                T temp = x[j];
                x[j] = x[i];
                x[i] = temp;
            }
        }
    }

    void bitRevShuffle(std::vector<NttType> &x);
}

#endif //HLS_YATFHE_NTT_H
