//
// Created by xintong on 4/21/25.
//

#ifndef HLS_YATFHE_CMUX_H
#define HLS_YATFHE_CMUX_H

#include "yatfhe/trgsw.h"

void controlMux(Trlwe& res, const Trlwe& input, int aBarI, const Trgsw& bskI, const YatfheParameters& param);

void controlMuxNtt(Trlwe& res, const Trlwe& input, int aBarI, const TrgswDft& bskI, const YatfheParameters& param);

void controlMuxNttMP21(Trlwe& res, const Trlwe& in0, const Trlwe& in1, const TrgswMPDft& sel, const YatfheParameters& param);

void controlMuxApproxCRT(std::vector<Trlwe8>& res, const std::vector<Trlwe8>& inputs, int aBarI, const std::vector<Trgsw8>& bskCRT, const YatfheParameters& param);

void controlMuxApproxCRTNtt(std::vector<Trlwe8>& res, const std::vector<Trlwe8>& inputs, int aBarI, const std::vector<TrgswDft24>& bskCRT, const YatfheParameters& param);


#endif //HLS_YATFHE_CMUX_H
