//
// Created by Xintong Song on 2023/12/25.
//

#ifndef HLS_YATFHE_BOOTSTRAP_H
#define HLS_YATFHE_BOOTSTRAP_H

#include <vector>
#include "tlwe.h"
#include "trgsw.h"
#include "yatfhe_parameters.h"

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

    explicit BootstrappingKey(const YatfheParameters& parameters) :
        n(parameters.n),
        unfolding(parameters.unfolding),
        bsk(parameters.n, Trgsw(parameters)),
        bskDft(parameters.n,TrgswDft(parameters)) {};
};

void trgswFunctionalBootstrapping(TrgswDft& out, const Tlwe& in, const BootstrappingKey& bsk, const YatfheParameters& parameters);

void calculateB(std::vector<uint64_t>& coeffsA, std::vector<uint64_t>& coeffsS, std::vector<uint64_t>& coeffsB, int N);

void initCoeffsViaUniformDistribution(std::vector<Torus>& coeffs, int N);

void addGaussianNoiseToCoeffs(std::vector<Torus>& coeffs, int N, double sigma);

void trgswEncZero(Trgsw& trgsw, TrgswDft& trgswDft, const YatfheParameters& yatfheParameters, const TrgswKey& trgswKey);

void newBootstrappingKeyWoUnfolding(BootstrappingKey& bsk, const YatfheParameters& yatfheParameters, const TrgswKey& trgswKey, const TlweKey& tlweKey);

void newBootstrappingKey(BootstrappingKey& bsk, const YatfheParameters& yatfheParameters, const TrgswKey& trgswKey, const TlweKey& tlweKey);

#endif //HLS_YATFHE_BOOTSTRAP_H
