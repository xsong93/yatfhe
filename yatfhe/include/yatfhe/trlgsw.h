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

struct TrlgswDft14 {
    std::vector<TrgswDft14> trgswDfts;
    int l;

    explicit TrlgswDft14(const YatfheParameters& param) :
            trgswDfts(param.l2, TrgswDft14(param)),
            l(param.l2) {};
};

void trgswEncryptNtt14(Trgsw& trgsw, TrlgswDft14& trgswDft14, const YatfheParameters& param, const TrgswKey& trgswKey, Integer mu);


#endif //HLS_YATFHE_TRLGSW_H
