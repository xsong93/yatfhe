//
// Created by xsong93 on 9/4/2026.
//
#include "include/bench_out.h"
#include <nlohmann/json.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <numeric>

using std::chrono::steady_clock;
using std::chrono::duration_cast;
using std::chrono::microseconds;

#include "yatfhe/tlwe.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/bootstrapping.h"
#include "yatfhe/blind_rotate.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/numeric.h"
#include "yautil/initializer.h"
#include "yautil/tool.h"

using json = nlohmann::json;

static double median(std::vector<double> v) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
}

//
// NTT-split sweep: scan moveNtt = 0..8
// Usage: ntt_split_ablate [--reps=N] [--out=DIR]
//
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

    const int M = 9;   // moveNtt = 0..8
    int unitsA = 2;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a.rfind("--unitsA=", 0) == 0) unitsA = std::atoi(a.c_str() + 9);
    }
    printf("unitsA=%d\n", unitsA);

    // correctness: every split must decrypt to the same plaintext as m=0
    {
        IntPolynomial ref{param.N};
        {
            Trlwe acc{param.k, param.N};
            blindRotateLazyPipeAltNttSweep(acc, bsk, sTlwe, v, gdVntt, 0, 2, param);
            symDecTrlweToInt(ref, acc, trlweKey, param.torusBase);
        }
        for (int m = 1; m < M; m++) {
            Trlwe acc{param.k, param.N};
            blindRotateLazyPipeAltNttSweep(acc, bsk, sTlwe, v, gdVntt, m, unitsA, param);
            IntPolynomial dec{param.N};
            symDecTrlweToInt(dec, acc, trlweKey, param.torusBase);
            const bool ok = dec.coeffs == ref.coeffs;
            printf("correctness moveNtt=%d: %s\n", m, ok ? "PASS" : "FAIL");
            if (!ok) return 1;
        }
    }

    std::vector<std::vector<double>> ts(M);
    for (auto& t : ts) t.reserve(reps);
    Trlwe acc{param.k, param.N};
    for (int rep = 0; rep < reps; rep++) {
        for (int m = 0; m < M; m++) {
            auto t0 = steady_clock::now();
            blindRotateLazyPipeAltNttSweep(acc, bsk, sTlwe, v, gdVntt, m, unitsA, param);
            auto t1 = steady_clock::now();
            ts[m].push_back(duration_cast<microseconds>(t1 - t0).count() / 1000.0);
        }
    }

    json out;
    out["reps"] = reps;
    double best = 1e9;
    int bestM = 0;
    for (int m = 0; m < M; m++) {
        out[std::to_string(m)] = {{"median_ms", median(ts[m])}};
        const double med = median(ts[m]);
        printf("moveNtt=%d  warm median %.3f ms  (vs m=0: %+.4f ms)\n",
               m, med, med - median(ts[0]));
        if (med < best) { best = med; bestM = m; }
    }
    printf("best: moveNtt=%d at %.3f ms\n", bestM, best);

    const std::string path = yabench::benchOutPath("ntt_split_results.json");
    std::ofstream ofs(path);
    ofs << out.dump(2) << std::endl;
    printf("wrote %s\n", path.c_str());
    return 0;
}
