//
// Created by Xintong Song on 2024/5/9.
//

#ifndef HLS_YATFHE_TRLGSW_H
#define HLS_YATFHE_TRLGSW_H

#include "yatfhe/trgsw.h"

struct Trlgsw {
    std::vector<Trgsw> trgsws;
    int l ;

    explicit Trlgsw(const YatfheParameters& param) :
        trgsws(param.l2, Trgsw(param)),
        l(param.l2) {};
};

#endif //HLS_YATFHE_TRLGSW_H
