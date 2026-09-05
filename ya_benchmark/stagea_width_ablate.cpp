//
// Created by xsong93 on 9/4/2026.
//
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>

#include "yatfhe/blind_rotate.h"
#include "yatfhe/bootstrapping.h"
#include "yatfhe/numeric.h"
#include "yatfhe/tlwe.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yautil/initializer.h"

#include "include/bench_out.h"
#include <nlohmann/json.hpp>

using std::chrono::duration_cast;
using std::chrono::microseconds;
using std::chrono::steady_clock;
using json = nlohmann::json;

static json summarize(std::vector<double> v) {
    json o;
    o["count"] = v.size();
    if (v.empty()) return o;
    std::sort(v.begin(), v.end());
    const double mean = std::accumulate(v.begin(), v.end(), 0.0) / v.size();
    double sq = 0.0;
    for (const double x : v) sq += (x - mean) * (x - mean);
    o["mean_ms"] = mean;
    o["sd_ms"] = std::sqrt(sq / v.size());
    o["median_ms"] = v[v.size() / 2];
    return o;
}

static json paired(const std::vector<double>& u2, const std::vector<double>& u4) {
    const auto m = std::min(u2.size(), u4.size());
    std::vector<double> d(m);
    for (size_t i = 0; i < m; i++) d[i] = u2[i] - u4[i];  // positive = 2 units faster
    return summarize(d);
}

// Paired A/B of the stageA dispatch width on the SHIPPED split (moveNtt=2):
// blindRotateLazyPipeAltNttSweep(..., 2, 2) vs (..., 2, 4).
int main(int argc, char** argv) {
    int reps = 2000;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a.rfind("--reps=", 0) == 0) reps = std::atoi(a.c_str() + 7);
        if (a.rfind("--out=", 0) == 0) yabench::setBenchOutDir(a.substr(6));
    }
    YatfheParameters param{};
    initYatfhe(param);
    printf("n=%d k=%d N=%d l=%d lApprox=%d reps=%d\n",
           param.n, param.k, param.N, param.l, param.lApprox, reps);

    TlweKey tlweKey{param.n, param.lweNoiseB};
    TrgswKey trgswKey{param};
    genTlweKey(tlweKey);
    genTrlweKey(trgswKey.trlweKey);
    TrlweKey& trlweKey = trgswKey.trlweKey;

    BootstrappingKeyMPLazyPipeAlt bsk{param, param.lApprox, true};
    genBootstrappingKeyMPLazyPipeAlt(bsk, trgswKey, tlweKey, param);

    TorusPolynomial v{param.N};
    generateTestPolynomialFR(v, param.torusBase, 2 * param.N);
    vector<NttPolynomial> gdVntt;
    prepareGdV(gdVntt, v, param);

    ScaledTlwe sTlwe{2 * param.N, param.n};
    for (int i = 0; i < sTlwe.n; i++) sTlwe.a[i] = genIntUniformDist(-2 * param.N, 2 * param.N);
    sTlwe.b = genIntUniformDist(-2 * param.N, 2 * param.N);

    const auto run = [&](const int unitsA, Trlwe& acc) {
        blindRotateLazyPipeAltNttSweep(acc, bsk, sTlwe, v, gdVntt, 2, unitsA, param);
    };

    // one-shot correctness
    {
        Trlwe accA{param.k, param.N}, accB{param.k, param.N};
        run(2, accA);
        run(4, accB);
        IntPolynomial decA{param.N}, decB{param.N};
        symDecTrlweToInt(decA, accA, trlweKey, param.torusBase);
        symDecTrlweToInt(decB, accB, trlweKey, param.torusBase);
        const bool ok = decA.coeffs == decB.coeffs;
        printf("correctness: %s\n", ok ? "PASS" : "FAIL");
        if (!ok) return 1;
    }

    {
        Trlwe acc{param.k, param.N};
        for (int k = 0; k < 5; k++) {
            run(2, acc);
            run(4, acc);
        }
    }

    std::vector<double> u2, u4;
    u2.reserve(reps);
    u4.reserve(reps);
    Trlwe acc{param.k, param.N};
    for (int rep = 0; rep < reps; rep++) {
        const bool u2First = (rep % 2 == 0);
        for (int step = 0; step < 2; step++) {
            const bool runU2 = (step == 0) == u2First;
            auto t0 = steady_clock::now();
            run(runU2 ? 2 : 4, acc);
            auto t1 = steady_clock::now();
            (runU2 ? u2 : u4).push_back(duration_cast<microseconds>(t1 - t0).count() / 1000.0);
        }
    }

    json out;
    out["reps"] = reps;
    out["unitsA2_ms"] = summarize(u2);
    out["unitsA4_ms"] = summarize(u4);
    out["paired_delta_ms"] = paired(u2, u4);

    printf("unitsA=2: med %.3f mean %.3f sd %.3f\n", out["unitsA2_ms"]["median_ms"].get<double>(),
           out["unitsA2_ms"]["mean_ms"].get<double>(), out["unitsA2_ms"]["sd_ms"].get<double>());
    printf("unitsA=4: med %.3f mean %.3f sd %.3f\n", out["unitsA4_ms"]["median_ms"].get<double>(),
           out["unitsA4_ms"]["mean_ms"].get<double>(), out["unitsA4_ms"]["sd_ms"].get<double>());
    printf("paired delta (u2-u4): med %.4f mean %.4f\n",
           out["paired_delta_ms"]["median_ms"].get<double>(), out["paired_delta_ms"]["mean_ms"].get<double>());

    std::ofstream ofs(yabench::benchOutPath("stagea_width_results.json"));
    ofs << out.dump(2) << "\n";
    return 0;
}
