//
// Created by xsong93 on 09/3/2026.
//

#include <nlohmann/json.hpp>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <numeric>
#include "include/bench_out.h"
#include "yatfhe/tlwe.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/bootstrapping.h"
#include "yatfhe/blind_rotate.h"
#include "yautil/ya_serializer.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/numeric.h"
#include "yautil/initializer.h"
#include "yautil/tool.h"

using std::chrono::steady_clock;
using std::chrono::duration_cast;
using std::chrono::microseconds;
using json = nlohmann::json;

static double median(std::vector<double> v) {
    if (v.empty()) {
        return 0.0;
    }
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
}

static json summarize(const std::vector<double>& v) {
    json o;
    o["count"] = v.size();
    if (v.empty()) {
        return o;
    }
    double mean = std::accumulate(v.begin(), v.end(), 0.0) / v.size();
    double sq = 0.0;
    for (const double x : v) {
        sq += (x - mean) * (x - mean);
    }
    o["mean_ms"] = mean;
    o["sd_ms"] = std::sqrt(sq / v.size());
    o["median_ms"] = median(v);
    return o;
}

// paired restructured-minus-plain deltas with SE
static json paired(const std::vector<double>& plain, const std::vector<double>& restruct) {
    const auto m = std::min(plain.size(), restruct.size());
    std::vector<double> d(m);
    for (size_t i = 0; i < m; i++) {
        d[i] = restruct[i] - plain[i];
    }
    json o = summarize(d);
    double mean = o["mean_ms"].get<double>();
    double se = o["sd_ms"].get<double>() / std::sqrt(static_cast<double>(m));
    o["se_ms"] = se;
    o["mean_plus_se_ms"] = mean + se;
    o["mean_minus_se_ms"] = mean - se;
    return o;
}

// paired micro-benchmark for the first-accumulation-step restructure
int main(int argc, char** argv) {
    int reps = 300;
    std::string keyFile;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a.rfind("--reps=", 0) == 0) {
            reps = std::atoi(a.c_str() + 7);
        }
        if (a.rfind("--out=", 0) == 0) {
            yabench::setBenchOutDir(a.substr(6));
        }
        if (a.rfind("--keyfile=", 0) == 0) {
            keyFile = a.substr(10);
        }
    }

    YatfheParameters param{};
    initYatfhe(param);
    printf("n=%d k=%d N=%d l=%d lApprox=%d reps=%d\n", param.n, param.k, param.N, param.l, param.lApprox, reps);

    // client side
    TlweKey tlweKey{param.n, param.lweNoiseB};
    TrgswKey trgswKey{param};
    genTlweKey(tlweKey);
    genTrlweKey(trgswKey.trlweKey);
    TrlweKey& trlweKey = trgswKey.trlweKey;

    BootstrappingKeyMPLazyPipeAlt bsk{param, param.lApprox, true};
    bool keyLoaded = false;
    {
        std::ifstream probe(keyFile);
        if (!keyFile.empty() && probe.good()) {
            probe.close();
            deserializeBskLazyPipeAlt(bsk, keyFile, param.n);
            keyLoaded = true;
            printf("key loaded from %s\n", keyFile.c_str());
        } else {
            genBootstrappingKeyMPLazyPipeAlt(bsk, trgswKey, tlweKey, param);
            if (!keyFile.empty()) {
                serializeBskLazyPipeAlt(bsk, keyFile);
                printf("key generated and written to %s\n", keyFile.c_str());
            }
        }
    }

    TorusPolynomial v{param.N};
    generateTestPolynomialFR(v, param.torusBase, 2 * param.N);

    vector<NttPolynomial> gdVntt;
    prepareGdV(gdVntt, v, param);

    ScaledTlwe sTlwe{2 * param.N, param.n};
    for (int i = 0; i < sTlwe.n; i++) {
        sTlwe.a[i] = genIntUniformDist(-2 * param.N, 2 * param.N);
    }
    sTlwe.b = genIntUniformDist(-2 * param.N, 2 * param.N);

    // both arms must decrypt to the same plaintext
    // skipped when the key was loaded from a file
    const bool loadedFromFile = !keyFile.empty() && keyLoaded;
    if (!loadedFromFile) {
        Trlwe accA{param.k, param.N}, accB{param.k, param.N};
        blindRotateLazyPipeAltNtt(accA, bsk, sTlwe, v, gdVntt, param);
        blindRotateLazyPipeAltNoFirstNtt(accB, bsk, sTlwe, v, param);
        IntPolynomial decA{param.N}, decB{param.N};
        symDecTrlweToInt(decA, accA, trlweKey, param.torusBase);
        symDecTrlweToInt(decB, accB, trlweKey, param.torusBase);
        const bool ok = decA.coeffs == decB.coeffs;
        printf("correctness: %s\n", ok ? "PASS" : "FAIL");
        if (!ok) {
            return 1;
        }
    }

    // warm up both arms
    {
        Trlwe acc{param.k, param.N};
        for (int k = 0; k < 5; k++) {
            blindRotateLazyPipeAltNtt(acc, bsk, sTlwe, v, gdVntt, param);
            blindRotateLazyPipeAltNoFirstNtt(acc, bsk, sTlwe, v, param);
        }
    }

    std::vector<double> plainSteady, restructSteady;
    plainSteady.reserve(reps);
    restructSteady.reserve(reps);

    Trlwe acc{param.k, param.N};
    for (int rep = 0; rep < reps; rep++) {
        const bool restructGoesFirst = (rep % 2 == 0);   // alternate order to cancel pair-order bias

        for (int step = 0; step < 2; step++) {
            const bool runRestruct = (step == 0) == restructGoesFirst;
            auto t0 = steady_clock::now();
            if (runRestruct) {
                blindRotateLazyPipeAltNtt(acc, bsk, sTlwe, v, gdVntt, param);
            } else {
                blindRotateLazyPipeAltNoFirstNtt(acc, bsk, sTlwe, v, param);
            }
            auto t1 = steady_clock::now();
            (runRestruct ? restructSteady : plainSteady).push_back(duration_cast<microseconds>(t1 - t0).count() / 1000.0);
        }
    }

    json out;
    out["reps"] = reps;
    out["steady"] = {
        {"plain_ms", summarize(plainSteady)},
        {"restructured_ms", summarize(restructSteady)},
        {"paired_delta_ms", paired(plainSteady, restructSteady)},
    };

    auto printCell = [](const char* name, const json& cell) {
        const auto& d = cell["paired_delta_ms"];
        printf("%-28s plain %.3f  restruct %.3f  delta %.4f +/- %.4f ms (medians: %.3f / %.3f)\n",
               name, cell["plain_ms"]["median_ms"].get<double>(),
               cell["restructured_ms"]["median_ms"].get<double>(),
               d["mean_ms"].get<double>(), d["se_ms"].get<double>(),
               cell["plain_ms"]["median_ms"].get<double>(),
               cell["restructured_ms"]["median_ms"].get<double>());
    };
    printCell("steady", out["steady"]);

    const std::string path = yabench::benchOutPath("firststep_fixedkey_results.json");
    std::ofstream ofs(path);
    ofs << out.dump(2) << std::endl;
    printf("wrote %s\n", path.c_str());
    return 0;
}
