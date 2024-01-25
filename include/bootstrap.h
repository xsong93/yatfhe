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
    std::vector<TrgswDft> bskDft {}; // n
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

void trgswFunctionalBootstrapping(TrgswDft& out, const Tlwe& in, const BootstrappingKey& bsk, Torus msg,
                                  const YatfheParameters& param);

void calModularInnerProduct(LagrangePolynomial& b, std::vector<TorusPolynomial>& a, const std::vector<IntPolynomial>& s,
                            int N, int k);

void modularAccumulate(std::vector<uint64_t>& coeffsB, std::vector<uint64_t>& coeffsA, std::vector<uint64_t>& coeffsS,
                       int N);

void initCoeffsViaUniformDistribution(std::vector<Torus>& coeffs, int N);

void initCoeffsWithGaussianNoise(std::vector<Torus>& coeffs, Torus msg, int N, double sigma);

void trgswEncZero(Trgsw& trgsw, TrgswDft& trgswDft, const YatfheParameters& param, const TrgswKey& trgswKey);

void bootstrappingKeyGenWoUnfolding(BootstrappingKey& bsk, const YatfheParameters& param, const TrgswKey& trgswKey,
                                    const TlweKey& tlweKey);

void bootstrappingKeyGen(BootstrappingKey& bsk, const YatfheParameters& param, const TrgswKey& trgswKey, const TlweKey& tlweKey);

#endif //HLS_YATFHE_BOOTSTRAP_H
