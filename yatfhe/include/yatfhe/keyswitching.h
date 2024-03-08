//
// Created by Xintong Song on 2024/3/8.
//

#ifndef HLS_YATFHE_KEYSWITCHING_H
#define HLS_YATFHE_KEYSWITCHING_H

#include <vector>
#include "tlwe.h"

using namespace std;

struct TlweKeySwitchingKey {
    vector<vector<Tlwe>> decomposedKsk;
    int nCurrKey {};
    int nNewKey {};
    int level {};

    TlweKeySwitchingKey(int nCurr, int nNew, int l) :
            nCurrKey(nCurr),
            nNewKey(nNew),
            level(l),
            decomposedKsk(nCurr, vector<Tlwe>(l, Tlwe(nNew))) {};
};

void genTlweKeySwitchingKey(TlweKeySwitchingKey& ksk, const TlweKey& currKey, const TlweKey& newKey, const YatfheParameters& param);

void tlweKeySwitch(Tlwe& output, const TlweKeySwitchingKey& ksk, const Tlwe& input, const YatfheParameters& param);

#endif //HLS_YATFHE_KEYSWITCHING_H
