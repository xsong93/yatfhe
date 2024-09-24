//
// Created by Xintong Song on 2024/8/27.
//

#ifndef HLS_YATFHE_CRT_H
#define HLS_YATFHE_CRT_H

#include <iostream>
#include <vector>
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/trlwe.h"

int calApproxCRTError(const std::vector<int>& f_tilde, const std::vector<int>& coeffs);

void calGadgetVector(std::vector<long>& w, int Qlow, const std::vector<int>& highModuli);

void approxCRTDecomp(std::vector<std::vector<int8_t>>& f, const std::vector<int>& coeffs, int Qlow, const std::vector<int>& lowModuli, const std::vector<int>& highModuli);

int32_t approxCRTReconstructSingle(const std::vector<int8_t>& f, const YatfheParameters& param);

void approxCRTReconstructPoly(std::vector<int>& f_tilde, const std::vector<std::vector<int8_t>>& f, const std::vector<long>& w, long q);

void syncGadgetDecomp(std::vector<Trlwe8>& out, const std::vector<Trlwe8>& aux, const YatfheParameters& param);

void broadcastCRT(std::vector<std::vector<Trlwe8>>& out, const std::vector<Trlwe8>& in, const YatfheParameters& param);

#endif //HLS_YATFHE_CRT_H
