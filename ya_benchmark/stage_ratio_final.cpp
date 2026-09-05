//
// Created by xsong93 on 9/4/2026.
//
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <vector>

#include "yatfhe/blind_rotate.h"
#include "yatfhe/bootstrapping.h"
#include "yatfhe/numeric.h"
#include "yatfhe/tlwe.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yautil/initializer.h"
#include "yautil/multi_threading.h"

using namespace std::chrono;

// Measure the shipped pipeline's per-stage wall times:
// stageA (rotate + GD + 2 moved NTT levels, dispatched as two units of two levels),
// stageB (SS' with the moved levels arriving pre-NTT'd), and EP on the calling thread.
int main(int argc, char** argv) {
    int runs = 10;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a.rfind("--runs=", 0) == 0) runs = std::atoi(a.c_str() + 7);
    }
    YatfheParameters param{};
    initYatfhe(param);
    constexpr int moveNtt = 2, unitsA = 2;
    const auto level = param.lApprox, N = param.N, K = param.k;
    auto& pool = ThreadPool::instance();

    TlweKey tlweKey{param.n, param.lweNoiseB};
    TrgswKey trgswKey{param};
    genTlweKey(tlweKey);
    genTrlweKey(trgswKey.trlweKey);

    TorusPolynomial v{param.N};
    generateTestPolynomialFR(v, param.torusBase, 2 * param.N);

    BootstrappingKeyMPLazyPipeAlt key{param, level, true};
    genBootstrappingKeyMPLazyPipeAlt(key, trgswKey, tlweKey, param);

    Tlwe in{param.n};
    symEncTlwe(in, modSwitchToTorusGeneral(1, 2 * param.torusBase, LWE_Q), tlweKey);
    ScaledTlwe s{2 * param.N, param.n};
    rescaleTlweToNewMod(s, in);

    // ping-pong buffers, as in the pipeline
    std::vector holders(level, Trlwe{param});
    std::vector decompA0(level, std::vector(param.l, std::vector(K, DecompPolynomial{N})));
    std::vector decompA1(level, std::vector(param.l, std::vector(K, DecompPolynomial{N})));
    std::vector decompADft0(level, std::vector(moveNtt, std::vector(K, NttPolynomial{N})));
    std::vector decompADft1(level, std::vector(moveNtt, std::vector(K, NttPolynomial{N})));
    std::vector b0(level, TorusPolynomial{N});
    std::vector b1(level, TorusPolynomial{N});
    std::vector expanded0(level, TrgswMPDft{param, level});
    std::vector expanded1(level, TrgswMPDft{param, level});
    int keyIndex = 0, rotateBy = s.a[1];

    auto* currDecompA = &decompA0; auto* nextDecompA = &decompA1;
    auto* currDecompADft = &decompADft0; auto* nextDecompADft = &decompADft1;
    auto* currB = &b0; auto* nextB = &b1;
    auto* currExpanded = &expanded0; auto* nextExpanded = &expanded1;

    const auto stageA = [&](const int u) {
        for (int l = u; l < level; l += unitsA) {
            rotateTrlweMinusOneBPlusOne(holders[l], (*nextB)[l], key.bskPrime[keyIndex].cPrime[l], rotateBy,
                                        static_cast<Torus>(1) << (param.torusBits - (l + 1) * param.radixBits));
            gadgetDecomposeTrlweA((*nextDecompA)[l], holders[l].a, param);
            for (int gl = 0; gl < moveNtt; gl++) {
                for (int k1 = 0; k1 < param.k; k1++) {
                    NttHexl::applyNtt((*nextDecompADft)[l][gl][k1], (*nextDecompA)[l][gl][k1]);
                }
            }
        }
    };

    const auto stageB = [&](const int l) {
        clearTrlwe((*nextExpanded)[l].cPrime[l]);
        for (auto &item : (*nextExpanded)[l].c[l]) clearTrlwe(item);
        const auto K2 = param.k;
        const auto L = param.l;
        thread_local NttPolynomial aDft;
        if (aDft.N != param.N) aDft = NttPolynomial{param.N};
        auto& cDft = (*nextExpanded)[l].c[l];
        auto& cPrimeDft = (*nextExpanded)[l].cPrime[l];
        auto& decompA = (*currDecompA)[l];
        auto& decompADft = (*currDecompADft)[l];
        for (auto gl = 0; gl < L; gl++) {
            auto& s2l = key.s2Dft.trlweDfts[gl];
            for (auto k1 = 0; k1 < K2; k1++) {
                const NttPolynomial* aSrc;
                if (gl < moveNtt) {
                    aSrc = &decompADft[gl][k1];
                } else {
                    NttHexl::applyNtt(aDft, decompA[gl][k1]);
                    aSrc = &aDft;
                }
                NttHexl::calModularInnerProductNtt(cPrimeDft.a[k1], *aSrc, NttHexl::getNttGadgetRecomper(gl));
                for (auto k2 = 0; k2 < K2; k2++) {
                    NttHexl::calModularInnerProductNtt(cDft[k1].a[k2], *aSrc, s2l.a[k2]);
                }
                NttHexl::calModularInnerProductNtt(cDft[k1].b, *aSrc, s2l.b);
            }
        }
        NttHexl::applyNtt(cPrimeDft.b, (*currB)[l]);
        for (auto k1 = 0; k1 < K2; k1++) {
            for (auto k2 = 0; k2 < K2; k2++) {
                NttHexl::addNttPolynomial(cDft[k1].a[k2], cDft[k1].a[k2], cPrimeDft.b);
            }
        }
    };

    const auto ep = [&](Trlwe& acc, TrgswMPDft& gsw) {
        externalProductTrgswMPNttInPlace(acc, gsw, level, param);
    };

    const int reps = 400;
    std::vector<double> ra, rb, re;
    Trlwe acc{param.k, param.N};
    for (int run = 0; run < runs; run++) {
        // warm up
        for (int r = 0; r < 20; r++) {
            TaskGroup g;
            pool.run(g, unitsA, stageA);
            g.wait();
            pool.run(g, level, stageB);
            g.wait();
        }
        auto t0 = steady_clock::now();
        for (int r = 0; r < reps; r++) {
            keyIndex = r % (param.n - 1);
            TaskGroup g;
            pool.run(g, unitsA, stageA);
            g.wait();
        }
        auto t1 = steady_clock::now();
        ra.push_back(duration<double, std::micro>(t1 - t0).count() / reps);

        t0 = steady_clock::now();
        for (int r = 0; r < reps; r++) {
            TaskGroup g;
            pool.run(g, level, stageB);
            g.wait();
        }
        t1 = steady_clock::now();
        rb.push_back(duration<double, std::micro>(t1 - t0).count() / reps);

        t0 = steady_clock::now();
        for (int r = 0; r < reps; r++) {
            ep(acc, (*nextExpanded)[r % level]);
        }
        t1 = steady_clock::now();
        re.push_back(duration<double, std::micro>(t1 - t0).count() / reps);

        // flip buffers each run to keep the code paths symmetric
        std::swap(currDecompA, nextDecompA);
        std::swap(currDecompADft, nextDecompADft);
        std::swap(currB, nextB);
        std::swap(currExpanded, nextExpanded);
    }

    auto med = [](std::vector<double> x) {
        std::sort(x.begin(), x.end());
        return x[x.size() / 2];
    };
    std::vector<double> ratios;
    for (size_t i = 0; i < ra.size(); i++) ratios.push_back(rb[i] / re[i]);
    std::sort(ratios.begin(), ratios.end());
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "shipped pipeline stages, " << runs << " runs:\n";
    std::cout << "  stageA (rotate+GD+2 NTT, 2 units): " << med(ra) << " us\n";
    std::cout << "  stageB (SS', 6 own NTT levels)   : " << med(rb) << " us\n";
    std::cout << "  EP (external product)             : " << med(re) << " us\n";
    std::cout << "  ratio SS'/EP = " << med(ratios) << "  range "
              << ratios.front() << "--" << ratios.back() << "\n";
    return 0;
}
