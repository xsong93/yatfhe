#include <benchmark/benchmark.h>
#include <filesystem>
#include "yatfhe/blind_rotate.h"
#include "yatfhe/bootstrapping.h"
#include "yatfhe/numeric.h"
#include "yatfhe/tlwe.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yautil/initializer.h"
#include "yautil/ya_serializer.h"

class BlindRotateBenchmark : public benchmark::Fixture {
public:
    BlindRotateBenchmark() = default;

    void SetUp(const benchmark::State& state) override {
        param = YatfheParameters{};
        param.N = 1024;
        param.batchSize = 3;
        param.tasksPerThread = 1;
        initYatfhe(param);

        // key gen
        TlweKey tlweKey{param.n, param.lweStdDev};
        TrgswKey trgswKey{param};
        TrlweKey& trlweKey = trgswKey.trlweKey;
        TlweKeySwitchingKey ksKey{param};
        genTlweKey(tlweKey);
        genTrlweKey(trlweKey);
        TlweKey tlweKsKey = tlweKey;
        tlweKsKey.sigma = param.rlweStdDev;
        genTlweKeySwitchingKey(ksKey, trlweKey, tlweKsKey, param);
        s2Dft = TrlevDft{param, param.l};
        symEncTrlevWithKeyNtt(s2Dft, trlweKey, trlweKey.s, true, param);
        v = TorusPolynomial{param.N};
        generateTestPolynomial(v, param.torusBase, 2 * param.N);

        bskMPOpt = BootstrappingKeyMPOpt{param, param.lApprox, false};
        genBootstrappingKeyMPOpt(bskMPOpt, trgswKey, tlweKey, v, param);
        bskMPLazyPipe = BootstrappingKeyMPLazyPipe{param, param.lApprox, true, true};
        genBootstrappingKeyMPLazyPipe(bskMPLazyPipe, trgswKey, tlweKey, v, param);
        bskMPLazyOpt = BootstrappingKeyMPLazy{param, param.lApprox, true, true};
        genBootstrappingKeyMPLazy(bskMPLazyOpt, trgswKey, tlweKey, v, param);


        // data gen
        Integer pt = 3;
        Torus mu = modSwitchToTorusGeneral(pt, param.torusBase, LWE_Q);
        Tlwe input{param.n};
        symEncTlwe(input, mu, tlweKey);
        sTlwe = ScaledTlwe{param.N * 2, param.n};
        rescaleTlweToNewMod(sTlwe, input);
        acc = Trlwe{param.k, param.N};
        genNoiselessTrlweSample(acc, v, sTlwe);

        one = TrgswMPDft{param};
        encryptTrgswMPNtt(one, 1, trgswKey, 0, param);
        out = Trlwe{param};
    }

protected:
    YatfheParameters param;
    BootstrappingKeyMPOpt bskMPOpt;
    BootstrappingKeyMPLazyPipe bskMPLazyPipe;
    BootstrappingKeyMPLazy bskMPLazyOpt;
    TorusPolynomial v;
    Trlwe acc;
    Trlwe out;
    ScaledTlwe sTlwe;
    TrlevDft s2Dft;
    TrgswMPDft one;
};

BENCHMARK_DEFINE_F(BlindRotateBenchmark, GINX_OPT)(benchmark::State& state) {
    for (auto _ : state) {
        deserializeBskMPOpt(bskMPOpt, "BSK_GINX_OPT.bin", param.n);
        blindRotateOptNtt(out, bskMPOpt.bskFirst, bskMPOpt.bskDft, sTlwe, v, param);
    }
    std::filesystem::remove("BSK_GINX_OPT.bin");
}

BENCHMARK_DEFINE_F(BlindRotateBenchmark, OURS_NAIVE)(benchmark::State& state) {
    for (auto _ : state) {
        deserializeBskMPLazy(bskMPLazyOpt, "BSK_PAL_LAZY.bin", param.n);
        blindRotateLazyMTNtt(out, bskMPLazyOpt.bskFirst, bskMPLazyOpt.bskTrim, bskMPLazyOpt.bskDecompA,
                             sTlwe, v, s2Dft, param);
        bskMPLazyOpt.initialized = true;
        bskMPLazyOpt.bskDecompA.clear();
    }
    std::filesystem::remove("BSK_PAL_LAZY.bin");
}

BENCHMARK_DEFINE_F(BlindRotateBenchmark, OURS_PIPELINE)(benchmark::State& state) {
    for (auto _ : state) {
        blindRotateLazyPipeSerializationNtt(out, bskMPLazyPipe.bskFirst, bskMPLazyPipe.bskDft,
                                            bskMPLazyPipe.bskDecompA, sTlwe, v, s2Dft, one,
                                            "BSK_PIPE.bin", param);
        bskMPLazyPipe.initialized = true;
        bskMPLazyPipe.bskDecompA.clear();
    }
    std::filesystem::remove("BSK_PIPE.bin");
}

BENCHMARK_REGISTER_F(BlindRotateBenchmark, GINX_OPT)
    ->Unit(benchmark::kMicrosecond)
    ->Iterations(500)
    ->UseRealTime()
    ->MeasureProcessCPUTime();
BENCHMARK_REGISTER_F(BlindRotateBenchmark, OURS_NAIVE)
    ->Unit(benchmark::kMicrosecond)
    ->Iterations(500)
    ->UseRealTime()
    ->MeasureProcessCPUTime();
BENCHMARK_REGISTER_F(BlindRotateBenchmark, OURS_PIPELINE)
    ->Unit(benchmark::kMicrosecond)
    ->Iterations(500)
    ->UseRealTime()
    ->MeasureProcessCPUTime();

BENCHMARK_MAIN();