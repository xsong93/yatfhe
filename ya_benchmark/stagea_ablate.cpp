//
// Created by xsong93 on 09/3/2026.
//
#include "include/bench_out.h"
#include <nlohmann/json.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <numeric>

#include "yatfhe/tlwe.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/bootstrapping.h"
#include "yatfhe/blind_rotate.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/numeric.h"
#include "yautil/initializer.h"
#include "yautil/tool.h"

using std::chrono::steady_clock;
using std::chrono::duration_cast;
using std::chrono::microseconds;
using json = nlohmann::json;

static double median(std::vector<double> v) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
}

static json summarize(const std::vector<double>& v) {
    json o;
    o["count"] = v.size();
    if (v.empty()) return o;
    double mean = std::accumulate(v.begin(), v.end(), 0.0) / v.size();
    double sq = 0.0;
    for (const double x : v) sq += (x - mean) * (x - mean);
    o["mean_ms"] = mean;
    o["sd_ms"] = std::sqrt(sq / v.size());
    o["median_ms"] = median(v);
    return o;
}

static json paired(const std::vector<double>& a, const std::vector<double>& b) {
    const auto m = std::min(a.size(), b.size());
    std::vector<double> d(m);
    for (size_t i = 0; i < m; i++) d[i] = a[i] - b[i];
    json o = summarize(d);
    double mean = o["mean_ms"].get<double>();
    double se = o["sd_ms"].get<double>() / std::sqrt(static_cast<double>(m));
    o["se_ms"] = se;
    return o;
}

// Paired A/B for the NS' dispatch width
int main(int argc, char** argv) {
    int reps = 300;
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

    // one-shot correctness: both arms must decrypt to the same plaintext
    {
        Trlwe accA{param.k, param.N}, accB{param.k, param.N};
        blindRotateLazyPipeAltNtt(accA, bsk, sTlwe, v, gdVntt, param);
        blindRotateLazyPipeAltNttStageA2(accB, bsk, sTlwe, v, gdVntt, param);
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
            blindRotateLazyPipeAltNtt(acc, bsk, sTlwe, v, gdVntt, param);
            blindRotateLazyPipeAltNttStageA2(acc, bsk, sTlwe, v, gdVntt, param);
        }
    }

    std::vector<double> a4, a2;
    a4.reserve(reps);
    a2.reserve(reps);
    Trlwe acc{param.k, param.N};
    for (int rep = 0; rep < reps; rep++) {
        const bool a4First = (rep % 2 == 0);
        for (int step = 0; step < 2; step++) {
            const bool runA4 = (step == 0) == a4First;
            auto t0 = steady_clock::now();
            if (runA4) {
                blindRotateLazyPipeAltNtt(acc, bsk, sTlwe, v, gdVntt, param);
            } else {
                blindRotateLazyPipeAltNttStageA2(acc, bsk, sTlwe, v, gdVntt, param);
            }
            auto t1 = steady_clock::now();
            (runA4 ? a4 : a2).push_back(duration_cast<microseconds>(t1 - t0).count() / 1000.0);
        }
    }

    json out;
    out["reps"] = reps;
    // Arms: a4[] holds blindRotateLazyPipeAltNtt = Sweep(moveNtt=2, unitsA=2)(the shipped warm path)
    // a2[] holds blindRotateLazyPipeAltNttStageA2 = rotate+GD only, moveNtt=0, unitsA=2.
    // The paired delta isolates the NTT split (m=0 vs m=2), not the dispatch width.
    out["shipped_m2_u2_ms"] = summarize(a4);
    out["m0_u2_ms"] = summarize(a2);
    out["paired_delta_ms"] = paired(a2, a4);  // positive = shipped (m2) faster
    const auto& d = out["paired_delta_ms"];
    printf("shipped(m2,u2) %.3f  m0,u2 %.3f  delta(m0-shipped) %.4f +/- %.4f ms (medians: %.3f / %.3f)\n",
           out["shipped_m2_u2_ms"]["median_ms"].get<double>(),
           out["m0_u2_ms"]["median_ms"].get<double>(),
           d["mean_ms"].get<double>(), d["se_ms"].get<double>(),
           out["shipped_m2_u2_ms"]["median_ms"].get<double>(),
           out["m0_u2_ms"]["median_ms"].get<double>());

    const std::string path = yabench::benchOutPath("stagea_results.json");
    std::ofstream ofs(path);
    ofs << out.dump(2) << std::endl;
    printf("wrote %s\n", path.c_str());
    return 0;
}
