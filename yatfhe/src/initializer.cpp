//
// Created by Xintong Song on 2024/5/23.
//
#include <sstream>
#include "yautil/initializer.h"

void yatfheInit(YatfheParameters& param) {
    initGlobalParamsNtt64(param.N);
}