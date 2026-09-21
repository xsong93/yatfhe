//
// Created by xsong93 on 08/25/2026.
//

#include "include/bench_out.h"
#include <nlohmann/json.hpp>

#include <algorithm>
#include <numeric>
#include <chrono>

using std::chrono::steady_clock;
using std::chrono::duration_cast;
using std::chrono::microseconds;
#include <sys/stat.h>

#include "yatfhe/tlwe.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/bootstrapping.h"
#include "yatfhe/blind_rotate.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/numeric.h"
#include "yautil/initializer.h"
#include "yautil/tool.h"
#include "yautil/ya_serializer.h"

using json = nlohmann::json;

static long fileBytes(const std::string& f) {
    struct stat st{};
    return stat(f.c_str(), &st) == 0 ? st.st_size : -1;
}

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
    o["min_ms"] = v.front();
    o["median_ms"] = v[v.size() / 2];
    o["max_ms"] = v.back();
    return o;
}

template <typename F>
static std::vector<double> measureTime(const int reps, F&& body) {
    std::vector<double> out;
    out.reserve(reps);
    for (int i = 0; i < reps; ++i) {
        const auto t0 = steady_clock::now();
        body();
        const auto t1 = steady_clock::now();
        out.push_back(duration_cast<microseconds>(t1 - t0).count() / 1000.0);
    }
    return out;
}

int main(int argc, char** argv) {
    int warmReps = 30;
    int coldReps = 15;
    std::string keyPrefix;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a.rfind("--warmreps=", 0) == 0)
            warmReps = std::atoi(a.c_str() + 11);
        if (a.rfind("--coldreps=", 0) == 0)
            coldReps = std::atoi(a.c_str() + 11);
        if (a.rfind("--out=", 0) == 0)
            yabench::setBenchOutDir(a.substr(6));
        if (a.rfind("--keyprefix=", 0) == 0)
            keyPrefix = a.substr(12);
    }
    // With --keyprefix, key files are read from the given directory and never regenerated or rewritten.
    const auto keyPath = [&](const std::string& name) {
        return keyPrefix.empty() ? name : (keyPrefix + "/" + name);
    };

    YatfheParameters param{};
    initYatfhe(param);
    printf("n=%d k=%d N=%d T=%d b=%d l=%d lApprox=%d  warmReps=%d coldReps=%d\n",
           param.n, param.k, param.N, param.torusBits, param.radixBits, param.l, param.lApprox,
           warmReps, coldReps);

    // ---- client side ----
    TlweKey tlweKey{param.n, param.lweNoiseB};
    TrgswKey trgswKey{param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    TlweKeySwitchingKey ksKey{param};
    genTlweKey(tlweKey);
    genTrlweKey(trlweKey);
    TlweKey tlweKsKey = tlweKey;
    genTlweKeySwitchingKey(ksKey, trlweKey, tlweKsKey, param);
    TorusPolynomial v{param.N};
    generateTestPolynomialFR(v, param.torusBase, 2 * param.N);

    const int p = param.torusBase;
    Integer pt = -1;
    const int slot = (pt % p + p) % p;
    Torus mu = modSwitchToTorusGeneral(slot, 2 * p, LWE_Q);
    Tlwe input{param.n};
    symEncTlwe(input, mu, tlweKey);
    ScaledTlwe sTlwe{2 * param.N, param.n};
    rescaleTlweToNewMod(sTlwe, input);

    Trlwe acc{param};
    json results;
    results["case"] = json::array();

    auto record = [&](const std::string& id, const std::string& label, const std::string& notes,
                      const std::string& file, const std::vector<double>& warm,
                      const std::vector<double>& cold, const std::vector<double>& ilv) {
        json r;
        r["id"] = id;
        r["label"] = label;
        r["enables"] = notes;
        r["key_bytes"] = fileBytes(file);
        r["key_mb"] = fileBytes(file) / 1048576.0;
        r["warm"] = summarize(warm);
        if (!cold.empty()) r["cold"] = summarize(cold);
        if (!ilv.empty()) r["interleaved"] = summarize(ilv);
        results["case"].push_back(r);
    };

    // TFHE / GINX baseline
    {
        const std::string f = keyPath("COMP_A_GINX.bin");
        if (keyPrefix.empty()) {
            BootstrappingKeyMP key{param, param.lApprox};
            genBootstrappingKeyMP(key, trgswKey, tlweKey, param);
            serializeBskMP(key, f);
        }
        BootstrappingKeyMP srv;
        deserializeBskMP(srv, f, param.n);
        auto warm = measureTime(warmReps, [&] {
            genNoiselessTrlweSample(acc, v, sTlwe);
            blindRotateJP22Ntt(acc, srv, sTlwe, param);
        });
        auto cold = measureTime(coldReps, [&] {
            clearFileCache(f);
            BootstrappingKeyMP s2;
            deserializeBskMP(s2, f, param.n);
            genNoiselessTrlweSample(acc, v, sTlwe);
            blindRotateJP22Ntt(acc, s2, sTlwe, param);
        });
        record("A", "TFHE (GINX) baseline", "none", f, warm, cold, {});
    }

    // WWL+24
    {
        const std::string f = keyPath("COMP_C_WWL24.bin");
        if (keyPrefix.empty()) {
            BootstrappingKeyWWL24 key{param, param.lApprox};
            genBootstrappingKeyWWL24(key, trgswKey, tlweKey, param);
            serializeBskWWL24(key, f);
        }
        BootstrappingKeyWWL24 srv;
        deserializeBskWWL24(srv, f, param.n);
        auto warm = measureTime(warmReps, [&] {
            genNoiselessTrlweSample(acc, v, sTlwe);
            blindRotateWWL24Ntt(acc, srv, sTlwe, param);
        });
        auto cold = measureTime(coldReps, [&] {
            clearFileCache(f);
            BootstrappingKeyWWL24 s2;
            deserializeBskWWL24(s2, f, param.n);
            genNoiselessTrlweSample(acc, v, sTlwe);
            blindRotateWWL24Ntt(acc, s2, sTlwe, param);
        });
        record("C", "WWL+24 succinct key", "compact key + on-demand scheme switching", f, warm, cold, {});
    }

    // on-the-fly NTT
    {
        const std::string f = keyPath("COMP_D_WWL24_ALT.bin");
        if (keyPrefix.empty()) {
            BootstrappingKeyWWL24Alt key{param, param.lApprox};
            genBootstrappingKeyWWL24Alt(key, trgswKey, tlweKey, param);
            serializeBskWWL24Alt(key, f);
        }
        BootstrappingKeyWWL24Alt srv;
        deserializeBskWWL24Alt(srv, f, param.n);
        auto warm = measureTime(warmReps, [&] {
            genNoiselessTrlweSample(acc, v, sTlwe);
            blindRotateWWL24AltNtt(acc, srv, sTlwe, param);
        });
        auto cold = measureTime(coldReps, [&] {
            clearFileCache(f);
            BootstrappingKeyWWL24Alt s2;
            deserializeBskWWL24Alt(s2, f, param.n);
            genNoiselessTrlweSample(acc, v, sTlwe);
            blindRotateWWL24AltNtt(acc, s2, sTlwe, param);
        });
        record("D", "+ on-the-fly NTT (plaintext key)", "avoids NTT key expansion", f, warm, cold, {});
    }

    // restructuring + naive parallel expansion
    {
        const std::string f = keyPath("COMP_E_PAR_LAZY.bin");
        if (keyPrefix.empty()) {
            BootstrappingKeyMPLazy key{param, param.lApprox, true, true};
            genBootstrappingKeyMPLazy(key, trgswKey, tlweKey, v, param);
            serializeBskMPLazy(key, f);
        }
        BootstrappingKeyMPLazyPipeAlt s2Src{param, param.lApprox, true};
        genBootstrappingKeyMPLazyPipeAlt(s2Src, trgswKey, tlweKey, param);
        Trlwe out{param};

        BootstrappingKeyMPLazy srvWarm;
        deserializeBskMPLazy(srvWarm, f, param.n);
        auto warm = measureTime(warmReps, [&] {
            blindRotateLazyMTNtt(out, srvWarm.bskFirst, srvWarm.bskTrim, srvWarm.bskDecompA, sTlwe, v, s2Src.s2Dft, param);
        });
        auto cold = measureTime(coldReps, [&] {
            clearFileCache(f);
            BootstrappingKeyMPLazy srv;
            deserializeBskMPLazy(srv, f, param.n);
            blindRotateLazyMTNtt(out, srv.bskFirst, srv.bskTrim, srv.bskDecompA, sTlwe, v, s2Src.s2Dft, param);
        });
        record("E", "parallel expansion (Alg. 3), LUT-dependent lazy key", "expansion off critical path via threads", f, warm, cold, {});
    }

    // OURS = pipeline + on-the-fly NTT
    {
        const std::string f = keyPath("COMP_G_PIPE_ALT.bin");
        if (keyPrefix.empty()) {
            BootstrappingKeyMPLazyPipeAlt key{param, param.lApprox, true};
            genBootstrappingKeyMPLazyPipeAlt(key, trgswKey, tlweKey, param);
            serializeBskLazyPipeAlt(key, f);
        }
        BootstrappingKeyMPLazyPipeAlt srv;
        deserializeBskLazyPipeAlt(srv, f, param.n);
        Trlwe out{param};
        vector<NttPolynomial> gdVntt;
        prepareGdV(gdVntt, v, param);
        auto warm = measureTime(warmReps, [&] {
            blindRotateLazyPipeAltNtt(out, srv, sTlwe, v, gdVntt, param);
        });
        auto cold = measureTime(coldReps, [&] {
            clearFileCache(f);
            BootstrappingKeyMPLazyPipeAlt s2;
            deserializeBskLazyPipeAlt(s2, f, param.n);
            blindRotateLazyPipeAltNtt(out, s2, sTlwe, v, gdVntt, param);
        });
        auto ilv = measureTime(coldReps, [&] {
            clearFileCache(f);
            BootstrappingKeyMPLazyPipeAlt s2;
            blindRotateLazyPipeAltInitNtt(out, s2, sTlwe, v, f, gdVntt, param);
        });
        record("G", "OURS: pipeline + on-the-fly NTT", "full design", f, warm, cold, ilv);
    }

    const std::string compPath = yabench::benchOutPath("comp_results.json");
    std::ofstream out(compPath);
    out << results.dump(2) << std::endl;
    return 0;
}
