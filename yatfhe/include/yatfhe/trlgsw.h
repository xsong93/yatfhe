//
// Created by Xintong Song on 2024/5/9.
//

#ifndef YATFHE_TRLGSW_H
#define YATFHE_TRLGSW_H

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

void encryptTrlgswNtt(Trlgsw& trlgsw, TrlgswDft14& trgswDft14, const YatfheParameters& param, const TrgswKey& trgswKey, Integer mu);

void externalProductTrlgsw(Trlwe& output, const Trlgsw& trlgsw, Trlwe& trlweInput, const YatfheParameters& param);

void externalProductTrlgswNtt(Trlwe& output, const TrlgswDft14& trlgswDft14Input, Trlwe& trlweInput, const YatfheParameters& param);


#endif //YATFHE_TRLGSW_H
