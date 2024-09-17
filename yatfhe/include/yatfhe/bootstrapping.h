//
// Created by Xintong Song on 2023/12/25.
//

#ifndef HLS_YATFHE_BOOTSTRAPPING_H
#define HLS_YATFHE_BOOTSTRAPPING_H

#include <vector>
#include "yatfhe/tlwe.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/keyswitching.h"

using namespace std;

struct BootstrappingKey {
    vector<TrgswDft> bskDft {}; // n
    vector<Trgsw> bsk {}; // n
    int n {};
    int k {};
    int N {};
    int l {};
    int unfolding {};

    explicit BootstrappingKey(const YatfheParameters& parameters) :
        n(parameters.n),
        unfolding(parameters.unfolding),
        bsk(parameters.n, Trgsw(parameters)),
        bskDft(parameters.n,TrgswDft(parameters)) {};
};

struct BootstrappingKeyCRT {
    std::vector<std::vector<TrgswDft24>> bskCRT {}; // n * d
    int n {};
    int d {};

    explicit BootstrappingKeyCRT(const YatfheParameters& param) :
            n(param.n),
            bskCRT(param.n, std::vector<TrgswDft24>(param.d, TrgswDft24(param))) {};
};

void trgswFunctionalBootstrapping(Tlwe& out, const Tlwe& input, const BootstrappingKey& bsk, const TlweKeySwitchingKey& ksk, const TorusPolynomial& v, const YatfheParameters& param);

void trgswFunctionalBootstrappingNtt(Tlwe& out, const Tlwe& input, const BootstrappingKey& bsk, const TlweKeySwitchingKey& ksk, const TorusPolynomial& v, const YatfheParameters& param);

void blindRotate(Trlwe& accum, const BootstrappingKey& bsk, const ScaledTlwe& input, const YatfheParameters& param);

void blindRotateNtt(Trlwe& accum, const BootstrappingKey& bsk, const ScaledTlwe& input, const YatfheParameters& param);

void controlMux(Trlwe& res, const Trlwe& input, int aBarI, const Trgsw& bskI, const YatfheParameters& param);

void controlMuxNtt(Trlwe& res, const Trlwe& input, int aBarI, const TrgswDft& bskI, const YatfheParameters& param);

void controlMuxCRT(Trlwe8& res, const std::vector<Trlwe8>& inputs, int aBarI, const std::vector<TrgswDft24>& bskCRT, const YatfheParameters& param);

void bootstrappingKeyGenWoUnfolding(BootstrappingKey& bsk, const YatfheParameters& param, TrgswKey& trgswKey, const TlweKey& tlweKey);

void bootstrappingKeyGen(BootstrappingKey& bsk, const YatfheParameters& param, TrgswKey& trgswKey, const TlweKey& tlweKey);

void bootstrappingKeyCRTDecomp(BootstrappingKeyCRT& bskCRT, const BootstrappingKey& bsk, const YatfheParameters& param);

#endif //HLS_YATFHE_BOOTSTRAPPING_H
