#include <chrono>
#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>

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

int main() {
    YatfheParameters param{};
    initYatfhe(param);
    const auto level = param.lApprox, N = param.N, K = param.k;

    std::cout << "n=" << param.n << " N=" << N << " l=" << param.l
              << " lApprox=" << level
              << "   stageA/stageB run " << level
              << " tasks in parallel\n\n";

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

    auto& pool = ThreadPool::instance();
    std::vector holders(level, Trlwe{param});
    std::vector decompA(level, std::vector(param.l, std::vector(K, DecompPolynomial{N})));
    std::vector bb(level, TorusPolynomial{N});
    TrgswMPDft expanded{param, level};
    auto& s2 = key.s2Dft;
    int keyIndex = 0, rotateBy = s.a[1];

    const auto stageA = [&](const int l) {
        rotateTrlweMinusOneBPlusOne(holders[l], bb[l], key.bskPrime[keyIndex].cPrime[l], rotateBy,
            static_cast<Torus>(1) << (param.torusBits - (l + 1) * param.radixBits));
        gadgetDecomposeTrlweA(decompA[l], holders[l].a, param);
    };

    const auto stageB = [&](const int l) {
        clearTrlwe(expanded.cPrime[l]);
        for (auto& item : expanded.c[l]) {
            clearTrlwe(item);
        }
        switchTrlweToSecretEmbeddingNttMix(expanded.c[l], expanded.cPrime[l], decompA[l], bb[l], s2, param);
    };

    const int reps = 400;

    // Warm-up and measure stageA
    for (int r = 0; r < 20; ++r) {
        TaskGroup g;
        pool.run(g, level, stageA);
        g.wait();
    }
    auto t0 = steady_clock::now();
    for (int r = 0; r < reps; ++r) {
        keyIndex = r % (param.n - 1);
        TaskGroup g;
        pool.run(g, level, stageA);
        g.wait();
    }
    auto t1 = steady_clock::now();
    double A = duration<double, std::micro>(t1 - t0).count() / reps;

    // Measure stageB
    for (int r = 0; r < 20; ++r) {
        TaskGroup g;
        pool.run(g, level, stageB);
        g.wait();
    }
    t0 = steady_clock::now();
    for (int r = 0; r < reps; ++r) {
        TaskGroup g;
        pool.run(g, level, stageB);
        g.wait();
    }
    t1 = steady_clock::now();
    double B = duration<double, std::micro>(t1 - t0).count() / reps;

    // Measure clear operations
    const auto clearsOnly = [&](const int l) {
        clearTrlwe(expanded.cPrime[l]);
        for (auto& item : expanded.c[l]) {
            clearTrlwe(item);
        }
    };
    for (int r = 0; r < 20; ++r) {
        TaskGroup g;
        pool.run(g, level, clearsOnly);
        g.wait();
    }
    t0 = steady_clock::now();
    for (int r = 0; r < reps; ++r) {
        TaskGroup g;
        pool.run(g, level, clearsOnly);
        g.wait();
    }
    t1 = steady_clock::now();
    double Cl = duration<double, std::micro>(t1 - t0).count() / reps;

    // Measure thread‑pool dispatch overhead
    const auto noop = [&](const int) {};
    for (int r = 0; r < 20; ++r) {
        TaskGroup g;
        pool.run(g, level, noop);
        g.wait();
    }
    t0 = steady_clock::now();
    for (int r = 0; r < reps; ++r) {
        TaskGroup g;
        pool.run(g, level, noop);
        g.wait();
    }
    t1 = steady_clock::now();
    double Ov = duration<double, std::micro>(t1 - t0).count() / reps;

    // move L NTTs from stageB to stageA
    std::vector dftBuf(level, std::vector(param.l, NttPolynomial{N}));

    const auto stageAp = [&](const int l) {
        rotateTrlweMinusOneBPlusOne(holders[l], bb[l], key.bskPrime[keyIndex].cPrime[l], rotateBy,
            1 << (param.torusBits - (l + 1) * param.radixBits));
        gadgetDecomposeTrlweA(decompA[l], holders[l].a, param);
        for (int j = 0; j < param.l; ++j) {
            NttHexl::applyNtt(dftBuf[l][j], decompA[l][j][0]);
        }
    };

    // Warm‑up
    for (int r = 0; r < 20; ++r) {
        TaskGroup g;
        pool.run(g, level, stageAp);
        g.wait();
    }
    t0 = steady_clock::now();
    for (int r = 0; r < reps; ++r) {
        keyIndex = r % (param.n - 1);
        TaskGroup g;
        pool.run(g, level, stageAp);
        g.wait();
    }
    t1 = steady_clock::now();
    double Ap = duration<double, std::micro>(t1 - t0).count() / reps;

    // stageBp: pointwise products only
    const auto stageBp = [&](const int l) {
        clearTrlwe(expanded.cPrime[l]);
        for (auto& item : expanded.c[l]) clearTrlwe(item);
        for (int j = 0; j < param.l; ++j) {
            const auto& aD = dftBuf[l][j];
            NttHexl::calModularInnerProductNtt(expanded.cPrime[l].a[0], aD,NttHexl::getNttGadgetRecomper(j));
            NttHexl::calModularInnerProductNtt(expanded.c[l][0].a[0], aD, s2.trlweDfts[j].a[0]);
            NttHexl::calModularInnerProductNtt(expanded.c[l][0].b, aD, s2.trlweDfts[j].b);
        }
        NttHexl::applyNtt(expanded.cPrime[l].b, bb[l]);
    };

    for (int r = 0; r < 20; ++r) {
        TaskGroup g;
        pool.run(g, level, stageBp);
        g.wait();
    }
    t0 = steady_clock::now();
    for (int r = 0; r < reps; ++r) {
        TaskGroup g;
        pool.run(g, level, stageBp);
        g.wait();
    }
    t1 = steady_clock::now();
    double Bp = duration<double, std::micro>(t1 - t0).count() / reps;

    // EP
    Trlwe acc{param};
    genNoiselessTrlweSample(acc, v, s);
    for (int r = 0; r < 20; ++r) {
        externalProductTrgswMPNttInPlace(acc, expanded, level, param);
    }
    t0 = steady_clock::now();
    for (int r = 0; r < reps; ++r) {
        externalProductTrgswMPNttInPlace(acc, expanded, level, param);
    }
    t1 = steady_clock::now();
    double E = duration<double, std::micro>(t1 - t0).count() / reps;

    // output
    std::cout << "stageA  NS' : rotate + gadgetDecomposeTrlweA " << A << " us\n";
    std::cout << "stageB  SS' : switchTrlweToSecretEmbeddingNttMix " << B << " us\n";
    std::cout << "EP          : externalProductTrgswMPNttInPlace  " << E << " us\n";
    std::cout << "\n  ratio SS'/EP = " << B / E << "  -> " << (B < E ? "CLAIM HOLDS" : "CLAIM FAILS") << "\n";
    std::cout << "pipeline floor = max(A,B,EP) = " << std::max({A, B, E}) << " us/component  (stage: "
              << ((A >= B && A >= E) ? "NS'" : ((B >= E) ? "SS'" : "EP")) << ")\n";

    std::cout << "\n  move the L NTTs from stageB into stageA\n";
    std::cout << "stageA' rotate + GD + " << param.l << " NTTs" << Ap << " us   (was " << A << ")\n";
    std::cout << "stageB' pointwise products only " << Bp << " us   (was " << B << ")\n";
    std::cout << "EP " << E << " us\n";
    std::cout << "pipeline floor " << std::max({A, B, E})
              << " -> " << std::max({Ap, Bp, E}) << " us/component ("
              << 100 * (1 - std::max({Ap, Bp, E}) / std::max({A, B, E})) << "% faster)\n";
    std::cout << "SS'/EP ratio " << B / E << " -> " << Bp / E
              << "claim: " << (Bp < E ? "WOULD HOLD" : "still fails") << "\n";

    std::cout << "\nSS' BREAKDOWN\n";
    std::cout << "clearTrlwe calls alone " << Cl << " us\n";
    std::cout << "thread-pool dispatch+join (noop) " << Ov << " us\n";
    std::cout << "arithmetic (9 NTT + 24 pointwise) " << B - Cl << " us  [derived]\n";
    std::cout << "=> removing the clears would give SS' ~ "
              << B - (Cl - Ov) << " us, ratio " << (B - (Cl - Ov)) / E << "\n";

    return 0;
}