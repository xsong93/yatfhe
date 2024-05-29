//
// Created by Xintong Song on 2024/5/23.
//
#include "yautil/initializer.h"
#include "yatfhe/ntt16.h"
#include "yatfhe/ntt64.h"
#include "yautil/time_counter.h"

void yatfheInit(YatfheParameters& param) {
    initGlobalParamsNtt64(param.N);
    initGlobalParamsNtt16(param.N);
    COUNT_TIME("init timer", std::cout << std::endl;)
}