#include <benchmark/benchmark.h>
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

        // client side
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

        v = TorusPolynomial{param.N};
        generateTestPolynomial(v, param.torusBase, 2 * param.N);

        bskMP = BootstrappingKeyMP{param, param.lApprox};
        genBootstrappingKeyMP(bskMP, trgswKey, tlweKey, param);
        bskMPOpt = BootstrappingKeyMPOpt{param, param.lApprox, false};
        genBootstrappingKeyMPOpt(bskMPOpt, trgswKey, tlweKey, v, param);

        // data gen
        Integer pt = 3;
        Torus mu = modSwitchToTorusGeneral(pt, param.torusBase, LWE_Q);
        Tlwe input{param.n};
        symEncTlwe(input, mu, tlweKey);
        sTlwe = ScaledTlwe{param.N * 2, param.n};
        rescaleTlweToNewMod(sTlwe, input);
        acc = Trlwe{param.k, param.N};
        genNoiselessTrlweSample(acc, v, sTlwe);

        out = Trlwe{param};
    }

protected:
    YatfheParameters param;
    BootstrappingKeyMP bskMP;
    BootstrappingKeyMPOpt bskMPOpt;
    TorusPolynomial v;
    Trlwe acc;
    Trlwe out;
    ScaledTlwe sTlwe;
    TrlevDft s2Dft;
    TrgswMPDft one;
};

BENCHMARK_DEFINE_F(BlindRotateBenchmark, GINX)(benchmark::State& state) {
    for (auto _ : state) {
        blindRotateJP22Ntt(acc, bskMP.bskDft, sTlwe, param);
    }
}

BENCHMARK_DEFINE_F(BlindRotateBenchmark, GINX_OPT)(benchmark::State& state) {
    for (auto _ : state) {
        blindRotateOptNtt(out, bskMPOpt.bskFirst, bskMPOpt.bskDft, sTlwe, v, param);
    }
}

BENCHMARK_REGISTER_F(BlindRotateBenchmark, GINX)
    ->Unit(benchmark::kMicrosecond)
    ->Iterations(500);
BENCHMARK_REGISTER_F(BlindRotateBenchmark, GINX_OPT)
    ->Unit(benchmark::kMicrosecond)
    ->Iterations(500);

BENCHMARK_MAIN();