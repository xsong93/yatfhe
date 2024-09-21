//
// Created by Xintong Song on 2024/5/23.
//
#include "yautil/initializer.h"
#include "yatfhe/ntt14.h"
#include "yatfhe/ntt16.h"
#include "yatfhe/ntt24.h"
#include "yatfhe/ntt64.h"
#include "yautil/time_counter.h"
#include "yatfhe/numeric_functions.h"

void yatfheInit(YatfheParameters& param) {
    initGlobalParamsNtt64(param.N);
    initGlobalParamsNtt24(param.N);
    initGlobalParamsNtt16(param.N);
    initGlobalParamsNtt14(param.N);
    for (size_t d = 0; d < param.dl; d++) {
        auto& dh = param.dh;
        param.taoU[dh + d] = static_cast<int>((modInverse(param.qLow / param.ql[d], param.ql[d])));
        param.taoUInv[dh + d] = static_cast<int>((longModP(param.qLow / param.ql[d], param.ql[d])));
    }
    COUNT_TIME("init timer", std::cout << std::endl;)
}