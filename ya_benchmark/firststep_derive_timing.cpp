//
// Created by xsong93 on 09/3/2026.
//
#include <chrono>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <cmath>
#include "yatfhe/blind_rotate.h"
#include "yautil/initializer.h"
#include "yautil/tool.h"
#include "include/bench_out.h"
#include <fstream>
#include <nlohmann/json.hpp>
#include "yatfhe/trgsw.h"
#include "yatfhe/tlwe.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/bootstrapping.h"
#include "yatfhe/polynomial.h"
#include "yautil/tool.h"

using namespace std;
using namespace std::chrono;

static void report(const char* name, vector<double> v) {
    sort(v.begin(), v.end());
    double mean = 0; for (double x : v) mean += x; mean /= v.size();
    double sd = 0; for (double x : v) sd += (x - mean) * (x - mean);
    sd = sqrt(sd / v.size());
    printf("%-34s n=%zu  median %7.3f us  mean %7.3f us  sd %6.3f us  se %5.3f us\n",
           name, v.size(), v[v.size() / 2] * 1000, mean * 1000, sd * 1000, sd / sqrt((double)v.size()) * 1000);
}

int main(int argc, char** argv) {
    int reps = 20000;
    for (int i = 1; i < argc; ++i) {
        const string a = argv[i];
        if (a.rfind("--reps=", 0) == 0) {
            reps = atoi(a.c_str() + 7);
        }
        if (a.rfind("--out=", 0) == 0) {
            yabench::setBenchOutDir(a.substr(6));
        }
    }
    YatfheParameters param{};
    initYatfhe(param);

    TlweKey tlweKey{param.n, param.lweNoiseB};
    TrgswKey trgswKey{param};
    genTlweKey(tlweKey);
    genTrlweKey(trgswKey.trlweKey);

    BootstrappingKeyMPLazyPipeAlt bsk{param, param.lApprox, true};
    genBootstrappingKeyMPLazyPipeAlt(bsk, trgswKey, tlweKey, param);

    TorusPolynomial v{param.N};
    generateTestPolynomialFR(v, param.torusBase, 2 * param.N);
    vector<NttPolynomial> gdVntt;
    prepareGdV(gdVntt, v, param);

    printf("l=%d lApprox=%d k=%d N=%d | bskPrime %zu, cPrime[0] a %zu b %d, gdVntt %zu (N=%d)\n",
           param.l, param.lApprox, param.k, param.N,
           bsk.bskPrime.size(), bsk.bskPrime[0].cPrime.size(),
           bsk.bskPrime[0].cPrime.empty() ? -1 : bsk.bskPrime[0].cPrime[0].b.N,
           gdVntt.size(), gdVntt.empty() ? -1 : gdVntt[0].N);

    // warm the allocator and the NTT tables
    Trlwe out{param.k, param.N};
    for (int i = 0; i < 200; i++) {
        deriveFirstComponentNtt(out, bsk.bskPrime[0].cPrime, gdVntt, param);
    }

    Torus sink = 0;
    vector<double> t;
    t.reserve(reps);
    for (int i = 0; i < reps; i++) {
        const auto& lev = bsk.bskPrime[i % bsk.bskPrime.size()].cPrime;
        const auto t0 = steady_clock::now();
        deriveFirstComponentNtt(out, lev, gdVntt, param);
        const auto t1 = steady_clock::now();
        sink += out.b.coeffs[i % param.N];
        __asm__ __volatile__("" ::: "memory");
        t.push_back(duration_cast<nanoseconds>(t1 - t0).count() / 1e9);
    }
    report("deriveFirstComponentNtt", t);
    double perCallUs = 0.0;
    {
        const auto a0 = steady_clock::now();
        for (int i = 0; i < reps; i++) {
            const auto& lev = bsk.bskPrime[i % bsk.bskPrime.size()].cPrime;
            deriveFirstComponentNtt(out, lev, gdVntt, param);
            sink += out.b.coeffs[i % param.N] & 1;
        }
        const auto a1 = steady_clock::now();
        const double total_us = duration_cast<nanoseconds>(a1 - a0).count() / 1000.0;
        perCallUs = total_us / reps;
        printf("%-34s total %9.1f ms for %d calls -> %7.3f us per call\n",
               "deriveFirstComponentNtt (loop wall clock)", total_us / 1000.0, reps, perCallUs);
    }
    {
        nlohmann::json j;
        j["reps"] = reps;
        j["derive_us_per_call"] = perCallUs;
        j["sink"] = (unsigned long long)sink;
        const std::string path = yabench::benchOutPath("derive_timing.json");
        std::ofstream ofs(path);
        ofs << j.dump(2) << std::endl;
        printf("wrote %s\n", path.c_str());
    }

    // calibration: one NTT pair is what the derivation does four times over
    TrlweDft dft{param.k, param.N};
    vector<double> u;
    u.reserve(reps / 10);
    for (int i = 0; i < reps / 10; i++) {
        const auto t0 = steady_clock::now();
        NttHexl::applyNttForAB(dft, out);
        const auto t1 = steady_clock::now();
        sink += dft.b.coeffs[i % param.N];
        __asm__ __volatile__("" ::: "memory");
        u.push_back(duration_cast<nanoseconds>(t1 - t0).count() / 1e9);
    }
    report("applyNttForAB (2 NTTs, calibration)", u);
    fprintf(stderr, "(checksum %llu)\n", (unsigned long long)sink);
    printf("(stage_ratio_final reports the SS' stage at 25.19 us on this machine)\n");
    return 0;
}
