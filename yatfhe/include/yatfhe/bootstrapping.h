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

struct BootstrappingKeyInternal {
    vector<vector<TrgswMP>> bsk {};
    vector<vector<TrgswMPDft>> bskDft {};
    int n {};

    explicit BootstrappingKeyInternal(const YatfheParameters& p) :
            bsk(p.n/2, vector<TrgswMP>(2, TrgswMP(p))),
            bskDft(p.n/2, vector<TrgswMPDft>(2, TrgswMPDft(p))),
            n {p.n} {}
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

void functionalBootstrapping(Tlwe& out, const Tlwe& input, const BootstrappingKey& bsk, const TlweKeySwitchingKey& ksk, const TorusPolynomial& v, const YatfheParameters& param);

void functionalBootstrappingNtt(Tlwe& out, const Tlwe& input, const BootstrappingKey& bsk, const TlweKeySwitchingKey& ksk, const TorusPolynomial& v, const YatfheParameters& param);

void functionalBootstrappingCrt(Tlwe& out, const Tlwe& input, const BootstrappingKeyCRT& bskCRT, const TlweKeySwitchingKey& ksk, const TorusPolynomial& v, const YatfheParameters& param);

void genBootstrappingKeyApproxCrt(BootstrappingKeyCRT& bskCRT, TrgswKey& trgswKey, const TlweKey& tlweKey, const YatfheParameters& param);

void genBootstrappingKey(BootstrappingKey& bsk, TrgswKey& trgswKey, const TlweKey& tlweKey, const YatfheParameters& param);

void genBootstrappingKeyInternal(BootstrappingKeyInternal& bsk, const TrgswKey& trgswKey, const TlweKey& tlweKey,
    const YatfheParameters& param);

void decompBootstrappingKeyMcrt(BootstrappingKeyCRT& bskCRT, const BootstrappingKey& bsk, const YatfheParameters& param);

#endif //HLS_YATFHE_BOOTSTRAPPING_H
