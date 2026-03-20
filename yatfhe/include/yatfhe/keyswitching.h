//
// Created by Xintong Song on 2024/3/8.
//

#ifndef HLS_YATFHE_KEYSWITCHING_H
#define HLS_YATFHE_KEYSWITCHING_H

#include <vector>
#include "yatfhe/tlwe.h"

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

void genTlweKeySwitchingKey(TlweKeySwitchingKey& ksk, const TrlweKey& currKey, const TlweKey& targetKey, const YatfheParameters& param);

void switchKeyForTlwe(Tlwe& output, const TlweKeySwitchingKey& ksk, const Tlwe& input, const YatfheParameters& param);

#endif //HLS_YATFHE_KEYSWITCHING_H
