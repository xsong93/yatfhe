//
// Created by Xintong Song on 2024/5/28.
//
#include "yatfhe/ntt.h"
#include "yatfhe/ntt64.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yautil/initializer.h"
#include "yautil/time_counter.h"

int main() {
    YatfheParameters param {};
    yatfheInit(param);
    printf("n:%d, k:%d, N:%d, b:%d, l:%d\n", param.n, param.k, param.N, param.radixBits, param.l);

    NttType a = 1234567890123456789ULL;
    NttType b = 9876543210987654321ULL;
    COUNT_TIME("modMul", {
        for (auto i = 0; i < 1000; i++) { modMul(a, b); }
    })
    COUNT_TIME("modmul", {
        for (auto i = 0; i < 1000; i++) { modmul(a, b); }
    })
    std::cout << "Result1: " << modMul(a, b) << std::endl;
    std::cout << "Result2: " << modmul(a, b) << std::endl;

}