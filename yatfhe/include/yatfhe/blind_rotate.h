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

void blindRotateWithPreRotNtt(Trlwe& accum, const vector<Trlwe>& trlwe, const vector<vector<TrgswMPDft>>& trgsws,
                              const ScaledTlwe& input, const YatfheParameters& param);

void blindRotateWithPreRotNttMT(Trlwe& accum, const vector<Trlwe>& trlwe, const vector<vector<TrgswMPDft>>& trgsws, const ScaledTlwe& input, const YatfheParameters& param);

void blindRotateOptNtt(Trlwe& accum, const vector<Trlwe>& bskFirst,  const vector<vector<TrgswMPDft>>& bskDft,
                       const ScaledTlwe& input, const TorusPolynomial& v, const YatfheParameters& param);

void blindRotateLazyNtt(Trlwe& accum, const vector<Trlwe>& bskFirst, vector<vector<TrgswMPDft>>& bsk, bool& initialized,
                        const ScaledTlwe& input, const TorusPolynomial& v, const TrlevDft& s2, const YatfheParameters& param);

void blindRotateLazyMTNtt(Trlwe& accum, const vector<Trlwe>& bskFirst, vector<vector<TrgswMPDft>>& bsk,
                           const vector<vector<vector<vector<vector<DecompPolynomial>>>>>& bskDecompA,
                           const ScaledTlwe& input, const TorusPolynomial& v, const TrlevDft& s2,
                           const YatfheParameters& param);

void blindRotateLazyPipeNtt(Trlwe& accum, const vector<Trlwe>& bskFirst, vector<vector<TrgswMPDft>>& bsk,
                           const vector<vector<vector<vector<vector<DecompPolynomial>>>>>& bskDecompA,
                           const ScaledTlwe& input, const TorusPolynomial& v, const TrlevDft& s2,
                           const TrgswMPDft& one, const YatfheParameters& param);

void blindRotateLazyPipeSerializationNtt(Trlwe& accum, const vector<Trlwe>& bskFirst, vector<vector<TrgswMPDft>>& bsk,
                                         const vector<vector<vector<vector<vector<DecompPolynomial>>>>>& bskDecompA,
                                         const ScaledTlwe& input, const TorusPolynomial& v, const TrlevDft& s2,
                                         const TrgswMPDft& one, const string& filename, const YatfheParameters& param);

void blindRotateExternalGeneralNtt(Trlev& accum, const vector<TrgswMPDft>& trgsws, const ScaledTlwe& input, const YatfheParameters& param);

void blindRotateJP22Ntt(Trlwe& accum, const vector<vector<TrgswMPDft>>& bskDft, const ScaledTlwe& input, const YatfheParameters& param);

void blindRotateJP22NttMT(Trlwe& accum, const vector<vector<TrgswMPDft>>& bskDft, const ScaledTlwe& input, const YatfheParameters& param);

void blindRotateMP21Ntt(Trlwe& accum, const vector<vector<TrgswMPDft>>& bskDft, const ScaledTlwe& input, const YatfheParameters& param);

#endif //HLS_YATFHE_BLIND_ROTATE_H
