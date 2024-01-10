//
// Created by Xintong Song on 2023/12/25.
//

#ifndef HLS_YATFHE_BOOTSTRAP_H
#define HLS_YATFHE_BOOTSTRAP_H

#include <vector>
#include "tlwe.h"
#include "trgsw.h"

struct BootstrappingKey {
//    TrgswDft* bskDft{ new TrgswDft};
    std::vector<TrgswDft> bskDft {}; // n
//    Trgsw* su{ new Trgsw};
    std::vector<Trgsw> bsk {}; // n
    int n {};
    int k {};
    int N {};
    int bgBit {};
    int l {};
    int unfolding {};
};

void trgswEncZero(Trgsw& trgsw, double alpha, const TrgswKey& trgswKey);

void newBootstrappingKeyWoUnfolding(BootstrappingKey& bsk, const TrgswKey& trgswKey, const TlweKey& tlweKey);

void newBootstrappingKey(BootstrappingKey& bsk, const TrgswKey& trgswKey, const TlweKey& tlweKey, int unfolding);

#endif //HLS_YATFHE_BOOTSTRAP_H
