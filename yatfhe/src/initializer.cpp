//
// Created by Xintong Song on 2024/5/23.
//
#include "yautil/initializer.h"
#include "yautil/time_counter.h"

void yatfheInit(YatfheParameters& param) {
    initGlobalParamsNtt64(param.N);
    COUNT_TIME("init timer", std::cout << std::endl;)
}