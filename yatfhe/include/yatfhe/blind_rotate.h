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

void blindRotateWithPreRotNtt(Trlwe& accum, const vector<Trlwe>& trlwe,
                              const vector<vector<TrgswMPDft>>& trgsws, const ScaledTlwe& input, int batchSize,
                              const YatfheParameters& param);

void blindRotateExternalGeneralNtt(Trlev& accum, const vector<TrgswMPDft>& trgsws, const ScaledTlwe& input, const YatfheParameters& param);

void blindRotateMPNtt(Trlwe& accum, const vector<TrgswMPDft>& bskDft, const ScaledTlwe& input, const YatfheParameters& param);

void blindRotateMPNtt(Trlwe& accum, const vector<vector<TrgswMPDft>>& bskDft, const ScaledTlwe& input, const YatfheParameters& param);

void blindRotateInternalPairWiseNtt(Trlwe& accum, vector<vector<TrgswMP>>& trgsws, vector<vector<TrgswMPDft>>& trgswDfts,
    const ScaledTlwe& input, int batchSize, const YatfheParameters& param);

void blindRotateInternalPairWiseAsymNtt(Trlwe& accum, vector<vector<Trlev>>& trlevs, vector<vector<TrgswMPDft>>& trgswDfts,
    const ScaledTlwe& input, const TrlevDft& sSquare, int batchSize, const YatfheParameters& param);

void blindRotateInternalPairWiseAsymOptNtt(Trlwe& out, vector<vector<Trlev>>& trlevs, vector<Trlwe>& trlwes,
                                           vector<vector<TrgswMPDft>>& trgswDfts, const ScaledTlwe& input,
                                           const TrlevDft& sSquare, int batchSize, const YatfheParameters& param);

#endif //HLS_YATFHE_BLIND_ROTATE_H
