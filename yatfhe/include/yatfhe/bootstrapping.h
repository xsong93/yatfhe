//
// Created by Xintong Song on 2023/12/25.
//

#ifndef HLS_YATFHE_BOOTSTRAPPING_H
#define HLS_YATFHE_BOOTSTRAPPING_H

#include <vector>
#include "yatfhe/tlwe.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/cmux.h"
#include "yatfhe/blind_rotate.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/keyswitching.h"
#include "yatfhe/key_patterns.h"

using namespace std;

struct BootstrappingKey {
    vector<TrgswDft> bskDft {};
    vector<Trgsw> bsk {};
    int n {};
    int group {};

    explicit BootstrappingKey(const YatfheParameters& p) : group(p.group) {
        if (group == 1) {
            n = p.n;
            bsk = vector<Trgsw>(n, Trgsw(p));
            bskDft = vector<TrgswDft>(n, TrgswDft(p));
        } else {
            n = p.n / group * (1 << group);
            bsk = vector<Trgsw>(n, Trgsw(p));
            bskDft = vector<TrgswDft>(n, TrgswDft(p));
        }
    };
};

struct BootstrappingKeyCRT {
    std::vector<std::vector<TrgswDft24>> bskCRT {}; // n * d
    std::vector<std::vector<Trgsw8>> bsk8 {}; // n * d
    int n {};
    int d {};

    explicit BootstrappingKeyCRT(const YatfheParameters& param) :
            n(param.n),
            bsk8(param.n, std::vector<Trgsw8>(param.d, Trgsw8(param.dh, param.k, param.N))),
            bskCRT(param.n, std::vector<TrgswDft24>(param.d, TrgswDft24(param.dh, param.k, param.N))) {};
};

void trgswFunctionalBootstrapping(Tlwe& out, const Tlwe& input, BootstrappingKey& bsk, const TlweKeySwitchingKey& ksk, const TorusPolynomial& v, const YatfheParameters& param);

void trgswFunctionalBootstrappingNtt(Tlwe& out, const Tlwe& input, BootstrappingKey& bsk, const TlweKeySwitchingKey& ksk, const TorusPolynomial& v, const YatfheParameters& param);

void trgswFunctionalBootstrappingCRT(Tlwe& out, const Tlwe& input, const BootstrappingKeyCRT& bskCRT, const TlweKeySwitchingKey& ksk, const TorusPolynomial& v, const YatfheParameters& param);

void bootstrappingKeyGenNormal(BootstrappingKey& bsk, TrgswKey& trgswKey, const TlweKey& tlweKey, const YatfheParameters& param);

void bootstrappingKeyGenGroup2(BootstrappingKey& bsk, TrgswKey& trgswKey, const TlweKey& tlweKey, const YatfheParameters& param);

void bootstrappingKeyGenApproxCRT(BootstrappingKeyCRT& bskCRT, TrgswKey& trgswKey, const TlweKey& tlweKey, const YatfheParameters& param);

void bootstrappingKeyGen(BootstrappingKey& bsk, TrgswKey& trgswKey, const TlweKey& tlweKey, const YatfheParameters& param);

void bootstrappingKeyMCRTDecomp(BootstrappingKeyCRT& bskCRT, const BootstrappingKey& bsk, const YatfheParameters& param);

#endif //HLS_YATFHE_BOOTSTRAPPING_H
