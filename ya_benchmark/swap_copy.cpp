//
// Created for anonymous review.
//

#include "yatfhe/trlwe.h"
#include "yautil/time_counter.h"
#include "yautil/tool.h"

int main(int argc, char **argv) {
    COUNT_TIME("init timer", std::cout<<std::endl;)
    const int N = 1024;
    const int k = 3;
    Trlwe t1(k, 1024);
    Trlwe t2(k, 1024);
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < k; i++) {
            t1.a[i].coeffs[j] = j + i;
        }
        t1.b.coeffs[j] = j;
    }
    COUNT_TIME("swap", std::swap(t2, t1);)
    COUNT_TIME("copyTrlwe", copyTrlwe(t2, t1, true, true);)
}