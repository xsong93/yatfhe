#include <benchmark/benchmark.h>
#include "yatfhe/bootstrapping.h"
#include "yatfhe/tlwe.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yautil/initializer.h"
#include "yautil/ya_serializer.h"
#include "yatfhe/blind_rotate.h"
#include "yautil/tool.h"

class ReadKeyBenchmark : public benchmark::Fixture {
public:
    ReadKeyBenchmark() = default;

    void SetUp(const benchmark::State& state) override {
        initYatfhe(param);

        // key gen
        TrlweKey& trlweKey = trgswKey.trlweKey;
        genTlweKey(tlweKey);
        genTrlweKey(trlweKey);
        symEncTrlevWithKeyNtt(s2Dft, trlweKey, trlweKey.s, true, param);
        generateTestPolynomial(v, param.torusBase, 2 * param.N);

        // data gen
        Integer pt = 3;
        Torus mu = modSwitchToTorusGeneral(pt, param.torusBase, LWE_Q);
        Tlwe input{param.n};
        symEncTlwe(input, mu, tlweKey);
        rescaleTlweToNewMod(sTlwe, input);
        genNoiselessTrlweSample(acc, v, sTlwe);
        encryptTrgswMPNtt(one, 1, trgswKey, 0, param);
    }

protected:
    YatfheParameters param{};
    TorusPolynomial v{param.N};
    TlweKey tlweKey{param.n, param.lweStdDev};
    TrgswKey trgswKey{param};
    Trlwe acc{param};
    Trlwe out{param};
    ScaledTlwe sTlwe{param.N * 2, param.n};
    TrlevDft s2Dft{param, param.l};
    TrgswMPDft one{param};
};

BENCHMARK_DEFINE_F(ReadKeyBenchmark, GINX)(benchmark::State& state) {
    BootstrappingKeyMP bskMP{param, param.lApprox};
    genBootstrappingKeyMP(bskMP, trgswKey, tlweKey, param);
    serializeBskMP(bskMP, "BSK_GINX.bin");
    for (auto _ : state) {
        state.PauseTiming();
        clearFileCache();
        state.ResumeTiming();
        BootstrappingKeyMP bskMPServer;
        deserializeBskMP(bskMPServer, "BSK_GINX.bin", param.n);
        blindRotateJP22Ntt(acc, bskMPServer.bskDft, sTlwe, param);
    }
}

BENCHMARK_DEFINE_F(ReadKeyBenchmark, GINX_OPT)(benchmark::State& state) {
    BootstrappingKeyMPOpt bskMPOpt{param, param.lApprox, false};
    genBootstrappingKeyMPOpt(bskMPOpt, trgswKey, tlweKey, v, param);
    serializeBskMPOpt(bskMPOpt, "BSK_GINX_OPT.bin");
    for (auto _ : state) {
        state.PauseTiming();
        clearFileCache();
        state.ResumeTiming();
        BootstrappingKeyMPOpt bskMPOptServer;
        deserializeBskMPOpt(bskMPOptServer, "BSK_GINX_OPT.bin", param.n);
        blindRotateOptNtt(out, bskMPOptServer.bskFirst, bskMPOptServer.bskDft,sTlwe, v, param);
    }
}

BENCHMARK_DEFINE_F(ReadKeyBenchmark, LAZY_PIPE)(benchmark::State& state) {
    BootstrappingKeyMPLazyPipe bskMPLazyPipe{param, param.lApprox, true, true};
    genBootstrappingKeyMPLazyPipe(bskMPLazyPipe, trgswKey, tlweKey, v, param);
    serializeBskLazyPipe(bskMPLazyPipe, "BSK_PIPE.bin");
    for (auto _ : state) {
        state.PauseTiming();
        clearFileCache();
        state.ResumeTiming();
        BootstrappingKeyMPLazyPipe bskMPLazyServer;
        deserializeBskLazyPipe(bskMPLazyServer, "BSK_PIPE.bin", param.n);
        blindRotateLazyPipeNtt(out, bskMPLazyServer.bskFirst, bskMPLazyServer.bskDft,
                               bskMPLazyServer.bskDecompA, sTlwe, v, s2Dft, one, param);
    }
}

BENCHMARK_REGISTER_F(ReadKeyBenchmark, GINX)
    ->Unit(benchmark::kMicrosecond)
    ->Iterations(100)
    ->UseRealTime();
BENCHMARK_REGISTER_F(ReadKeyBenchmark, GINX_OPT)
    ->Unit(benchmark::kMicrosecond)
    ->Iterations(100)
    ->UseRealTime();
BENCHMARK_REGISTER_F(ReadKeyBenchmark, LAZY_PIPE)
    ->Unit(benchmark::kMicrosecond)
    ->Iterations(100)
    ->UseRealTime();

BENCHMARK_MAIN();