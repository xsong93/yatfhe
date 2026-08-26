//
// Created by xsong93 on 08/26/2026.
//

#include <chrono>
#include <cstdio>
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

    std::cout << "n=" << param.n << " N=" << param.N << " l=" << param.l
              << " lApprox=" << param.lApprox
              << "  (SS decomposes to l, EP to lApprox)\n";

    TlweKey tlweKey{param.n, param.lweNoiseB};
    TrgswKey trgswKey{param};
    genTlweKey(tlweKey);
    genTrlweKey(trgswKey.trlweKey);

    BootstrappingKeyWWL24 key{param, param.lApprox};
    genBootstrappingKeyWWL24(key, trgswKey, tlweKey, param);

    const auto level = key.bskDft[0].l;
    auto& pool = ThreadPool::instance();
    std::vector cRows(level, std::vector(param.k, TrlweDft(param.k, param.N)));
    Trlwe acc{param}, tmp{param};
    TorusPolynomial v{param.N};
    generateTestPolynomialFR(v, param.torusBase, 2 * param.N);
    Tlwe in{param.n};
    symEncTlwe(in, modSwitchToTorusGeneral(1, 2 * param.torusBase, LWE_Q), tlweKey);
    ScaledTlwe s{2 * param.N, param.n};
    rescaleTlweToNewMod(s, in);
    genNoiselessTrlweSample(acc, v, s);

    const int reps = 400;

    // SS stage
    int keyIndex = 0;
    const auto schemeSwitch = [&](const int l) {
        switchTrlweToSecretEmbeddingNtt(cRows[l], key.bskDft[keyIndex].cPrime[l], key.s2Dft, param);
    };

    for (int r = 0; r < 20; r++) {
        TaskGroup g;
        pool.run(g, level, schemeSwitch);
        g.wait();
    }
    auto t0 = steady_clock::now();
    for (int r = 0; r < reps; r++) {
        keyIndex = r % (param.n - 1);
        TaskGroup g;
        pool.run(g, level, schemeSwitch);
        g.wait();
    }
    auto t1 = steady_clock::now();
    double ss = duration<double, std::micro>(t1 - t0).count() / reps;

    // EP stage
    for (int r = 0; r < 20; r++) {
        externalProductSplitNttInPlace(tmp, cRows, key.bskDft[0].cPrime, level,
                                       param);
    }
    t0 = steady_clock::now();
    for (int r = 0; r < reps; r++) {
        externalProductSplitNttInPlace(tmp, cRows,
                                       key.bskDft[r % (param.n - 1)].cPrime, level,
                                       param);
    }
    t1 = steady_clock::now();
    double ep = duration<double, std::micro>(t1 - t0).count() / reps;

    std::cout << std::fixed << std::setprecision(2);
    std::cout << "T_SS  (scheme switch, " << level
              << " parallel tasks, inner L=" << param.l << ") : "
              << std::setw(7) << ss << " us\n";
    std::cout << "T_EP  (external product, level=" << level << ") : "
              << std::setw(7) << ep << " us\n";
    std::cout << "ratio T_SS / T_EP = " << ss / ep << "\n";
    std::cout << "pipeline throughput floor = max(SS,EP) = "
              << std::max(ss, ep) << " us/component\n";
    std::cout << "ideal (fully hidden) = EP = " << ep << " us/component\n";

    return 0;
}