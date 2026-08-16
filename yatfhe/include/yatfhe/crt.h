//
// Created by Xintong Song on 2024/8/27.
//

#ifndef YATFHE_CRT_H
#define YATFHE_CRT_H

#include <iostream>
#include <vector>
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/trlwe.h"

int calApproxCRTError(const std::vector<int>& f_tilde, const std::vector<int>& coeffs);

void calGadgetVector(std::vector<long>& w, int Qlow, const std::vector<int>& highModuli);

void decompCrtExact(std::vector<std::vector<int8_t>>& f, const std::vector<int>& coeffs, const YatfheParameters& param);

void decompCrtExactIO(std::vector<std::vector<int8_t>>& f, const std::vector<int>& coeffs, const YatfheParameters& param);

void reconstructCrtExact(std::vector<int32_t>& f_tilde, const std::vector<std::vector<int8_t>>& f, const YatfheParameters& param);

void reconstructCrtExactIO(std::vector<int32_t>& f_tilde, const std::vector<std::vector<int8_t>>& f, const YatfheParameters& param);

void decompCrtApprox(std::vector<std::vector<int8_t>>& f, const std::vector<int>& coeffs, int Qlow, const std::vector<int>& lowModuli, const std::vector<int>& highModuli);

int32_t reconstructCrtApproxSingleEle(const std::vector<int8_t>& f, const YatfheParameters& param);

void reconstructCrtApproxVec(std::vector<int>& f_tilde, const std::vector<std::vector<int8_t>>& f, const std::vector<long>& w, long q);

void decompTrlweApproxCrt(std::vector<Trlwe8>& out, const std::vector<Trlwe8>& in, const YatfheParameters& param);

void broadcastTrlweApproxCrt(std::vector<std::vector<Trlwe8>>& out, const std::vector<Trlwe8>& in, const YatfheParameters& param);

#endif //YATFHE_CRT_H
