//
// Created by xintong on 4/21/25.
//

#ifndef HLS_YATFHE_BLIND_ROTATE_H
#define HLS_YATFHE_BLIND_ROTATE_H

#include "yatfhe/trgsw.h"

void blindRotate(Trlwe& accum, const vector<Trgsw>& bsk, const ScaledTlwe& input, const YatfheParameters& param);

void blindRotateNtt(Trlwe& accum, const vector<TrgswDft>& bskDft, const ScaledTlwe& input, const YatfheParameters& param);

void blindRotateApproxCRT(std::vector<Trlwe8>& accum, const vector<vector<Trgsw8>>& bskCRT, const ScaledTlwe& input, const YatfheParameters& param);

void blindRotateApproxCRTNtt(std::vector<Trlwe8>& accum, const vector<vector<TrgswDft24>>& bskCRT, const ScaledTlwe& input, const YatfheParameters& param);

#endif //HLS_YATFHE_BLIND_ROTATE_H
