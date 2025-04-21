//
// Created by xintong on 4/21/25.
//

#include "yatfhe/cmux.h"

// res = bsk * (c1 - c0) + c0 = bski * [ X^aBarI * input - input] + input
void controlMux(Trlwe& res, const Trlwe& input, const int aBarI, const Trgsw& bskI, const YatfheParameters& param) {
    Trlwe tmp {param.k, param.N};
    rotateTrlweMinusOne(tmp, input, aBarI); // res = c1 - c0 = X^aBarI * input - input
    externalProductTrgsw(res, bskI, tmp, param); // res *= bskI
    accumulateTrlwe(res, input); // res += input
}

// res = bsk * (c1 - c0) + c0 = bski * [ X^aBarI * input - input] + input
void controlMuxNtt(Trlwe& res, const Trlwe& input, const int aBarI, const TrgswDft& bskI, const YatfheParameters& param) {
    Trlwe tmp {param.k, param.N};
    rotateTrlweMinusOne(tmp, input, aBarI); // res = c1 - c0 = X^aBarI * input - input
    externalProductTrgswNtt(res, bskI, tmp, param); // res *= bskI
    accumulateTrlwe(res, input); // res += input
}

void controlMuxApproxCRT(std::vector<Trlwe8>& res, const std::vector<Trlwe8>& inputs, const int aBarI, const std::vector<Trgsw8>& bskCRT, const YatfheParameters& param) {
    std::vector<Trlwe8> tmp(param.d, Trlwe8{param.k, param.N});

    for (size_t i = 0; i < param.d; i++) {
        rotateTrlwe8MinusOne(tmp[i], inputs[i], aBarI, param.qd[i]); // res = c1 - c0 = X^aBarI * input - input
    }

    externalProductTrgswApproxCrt(res, bskCRT, tmp, param); // res *= bskI

    for (size_t i = 0; i < param.d; i++) {
        accumulateTrlweModP(res[i], inputs[i], param.qd[i]); // res += input
    }
}

void controlMuxApproxCRTNtt(std::vector<Trlwe8>& res, const std::vector<Trlwe8>& inputs, const int aBarI, const std::vector<TrgswDft24>& bskCRT, const YatfheParameters& param) {
    std::vector<Trlwe8> tmp(param.d, Trlwe8{param.k, param.N});

    for (size_t i = 0; i < param.d; i++) {
        rotateTrlwe8MinusOne(tmp[i], inputs[i], aBarI, param.qd[i]); // res = c1 - c0 = X^aBarI * input - input
    }

    externalProductTrgswApproxCrtNtt(res, bskCRT, tmp, param); // res *= bskI

    for (size_t i = 0; i < param.d; i++) {
        accumulateTrlweModP(res[i], inputs[i], param.qd[i]); // res += input
    }
}
