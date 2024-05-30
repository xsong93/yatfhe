//
// Created by Xintong Song on 2024/5/28.
//
#include "yatfhe/ntt.h"
#include "yatfhe/ntt14.h"
#include "yatfhe/ntt16.h"
#include "yatfhe/ntt64.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yautil/initializer.h"
#include "yautil/time_counter.h"

int main() {
    YatfheParameters param {};
    yatfheInit(param);
    printf("n:%d, k:%d, N:%d, b:%d, l:%d\n\n", param.n, param.k, param.N, param.radixBits, param.l);

    Ntt64 a = 1234567890123456789ULL;
    Ntt64 b = 9876543210987654321ULL;
    Ntt16 a1 = 59923;
    Ntt16 b1 = 65535;
    Ntt14 a14 = 59923;
    Ntt14 b14 = 65535;
    COUNT_TIME("modMul", {
        for (auto i = 0; i < 10000; i++) { modMULT64(a, b); }
    })
    COUNT_TIME("modmul", {
        for (auto i = 0; i < 10000; i++) { modmul64(a, b); }
    })
    COUNT_TIME("16", {
        for (auto i = 0; i < 10000; i++) { modMULT16(a1, b1); }
    })
    COUNT_TIME("16-8", {
        for (auto i = 0; i < 10000; i++) { modMult16(a1, b1); }
    })
    COUNT_TIME("14", {
        for (auto i = 0; i < 10000; i++) { modMULT14(a14, b14); }
    })
}