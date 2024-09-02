//
// Created by Xintong Song on 2024/8/27.
//

#ifndef HLS_YATFHE_CRT_H
#define HLS_YATFHE_CRT_H

#include <iostream>
#include <vector>

int calApproxCRTError(const std::vector<int>& f_tilde, const std::vector<int>& coeffs);

void calGadgetVector(std::vector<long>& w, int Qlow, const std::vector<int>& highModuli);

void approximateCRTDecomposition(std::vector<std::vector<int8_t>>& f, const std::vector<int>& coeffs, const int Qlow, const std::vector<int>& lowModuli, const std::vector<int>& highModuli);

void approximatePolynomialReconstruction(std::vector<int>& f_tilde, const std::vector<std::vector<int8_t>>& f, const std::vector<long>& w, long q);

#endif //HLS_YATFHE_CRT_H
