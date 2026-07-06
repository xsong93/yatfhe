//
// Created by Xintong Song on 2024/3/8.
//

#ifndef HLS_YATFHE_KEYSWITCHING_H
#define HLS_YATFHE_KEYSWITCHING_H

#include <vector>
#include "yatfhe/tlwe.h"
#include "yatfhe/trlev.h"

using namespace std;

struct TlweKeySwitchingKey {
    vector<vector<Tlwe>> decomposedKsk;
    int nCurrKey {};
    int nTargetKey {};
    int level {};

    explicit TlweKeySwitchingKey(const YatfheParameters& param) :
            nCurrKey(param.N * param.k),
            nTargetKey(param.n),
            level(param.ksLevel),
            decomposedKsk(param.N * param.k, vector<Tlwe>(param.ksLevel, Tlwe(param.n))) {};

    TlweKeySwitchingKey(int nCurr, int nTarget, int l) :
            nCurrKey(nCurr),
            nTargetKey(nTarget),
            level(l),
            decomposedKsk(nCurr, vector<Tlwe>(l, Tlwe(nTarget))) {};
};

struct PrivateKeySwitchingKey {
    vector<Trlev> trlevs;
    Trlev one;
    int level {};

    PrivateKeySwitchingKey(const YatfheParameters& param, const int level) :
        trlevs(vector<Trlev>(param.n, Trlev(param, level))),
        one(param, level),
        level(level) {};
};

void genTlweKeySwitchingKey(TlweKeySwitchingKey& ksk, const TrlweKey& currKey, const TlweKey& targetKey, const YatfheParameters& param);

void genPrivateKeySwitchingKey(PrivateKeySwitchingKey& psk, const TrlweKey& trlweKey, const TlweKey& tlweKey, const YatfheParameters& param);

void switchKeyForTlwe(Tlwe& output, const TlweKeySwitchingKey& ksk, const Tlwe& input, const YatfheParameters& param);

void tlweToTrlwePrivateKeySwitching(Trlwe& out, const Tlwe& in, const PrivateKeySwitchingKey& psk,
                                    const YatfheParameters& param);

#endif //HLS_YATFHE_KEYSWITCHING_H
