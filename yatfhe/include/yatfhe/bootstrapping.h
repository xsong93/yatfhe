//
// Created by Xintong Song on 2023/12/25.
//

#ifndef HLS_YATFHE_BOOTSTRAPPING_H
#define HLS_YATFHE_BOOTSTRAPPING_H

#include <vector>
#include "tlwe.h"
#include "trgsw.h"
#include "yatfhe_parameters.h"

using namespace std;

struct BootstrappingKey {
    vector<TrgswDft> bskDft {}; // n
    vector<Trgsw> bsk {}; // n
    int n {};
    int k {};
    int N {};
    int bgBit {};
    int l {};
    int unfolding {};

    explicit BootstrappingKey(const YatfheParameters& parameters) :
        n(parameters.n),
        unfolding(parameters.unfolding),
        bsk(parameters.n, Trgsw(parameters)),
        bskDft(parameters.n,TrgswDft(parameters)) {};
};

void trgswFunctionalBootstrapping(TrgswDft& out, const Tlwe& in, const BootstrappingKey& bsk, Torus msg,
                                  const YatfheParameters& param);

void blindRotate(Trlwe& accum, const BootstrappingKey& bsk, const NegaCyclicTlwe& sample, const YatfheParameters& param);

void muxRotate(Trlwe& res, Trlwe& accum, const TrgswDft& bski, int barai, const YatfheParameters& param);

void trgswMulToTrlwe(Trlwe& accum, const TrgswDft& bski, const YatfheParameters& param);

void trgswEncZero(Trgsw& trgsw, TrgswDft& trgswDft, const YatfheParameters& param, const TrgswKey& trgswKey);

void bootstrappingKeyGenWoUnfolding(BootstrappingKey& bsk, const YatfheParameters& param, const TrgswKey& trgswKey,
                                    const TlweKey& tlweKey);

void bootstrappingKeyGen(BootstrappingKey& bsk, const YatfheParameters& param, const TrgswKey& trgswKey, const TlweKey& tlweKey);

#endif //HLS_YATFHE_BOOTSTRAPPING_H
