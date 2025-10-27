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
        s2Dft = TrlevDft{param, param.l};
        symEncTrlevWithKeyNtt(s2Dft, trlweKey, trlweKey.s, true, param);
        v = TorusPolynomial{param.N};
        generateTestPolynomial(v, param.torusBase, 2 * param.N);

        bskMP = BootstrappingKeyMP{param, param.lApprox};
        genBootstrappingKeyMP(bskMP, trgswKey, tlweKey, param);
        bskMPOpt = BootstrappingKeyMPOpt{param, param.lApprox, false};
        genBootstrappingKeyMPOpt(bskMPOpt, trgswKey, tlweKey, v, param);
        bskMPLazy = BootstrappingKeyMPOpt{param, param.lApprox, true};
        genBootstrappingKeyMPOpt(bskMPLazy, trgswKey, tlweKey, v, param);
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
    BootstrappingKeyMP bskMP;
    BootstrappingKeyMPPreRot bskPre;
    BootstrappingKeyMPOpt bskMPOpt;
    BootstrappingKeyMPOpt bskMPLazy;
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
        auto bskMPOptServer = bskMPOpt;
        {
                vector<char> buffer;
                std::ostringstream oss1(std::ios::binary);
                serialize(bskMPOptServer.bskFirst[0], oss1);
                const std::string &data1 = oss1.str();
                buffer.insert(buffer.end(), data1.begin(), data1.end());
                for (auto i = 0; i < param.n - 1; i++) {
                    std::ostringstream oss(std::ios::binary);
                    serialize(bskMPOptServer.bskDft[i][0], oss);
                    const std::string &data = oss.str();
                    buffer.insert(buffer.end(), data.begin(), data.end());
                }
                std::ofstream outFile2("bsk_serialized_GINX_opt_bsk.bin", std::ios::binary | std::ios::app);
                outFile2.write(buffer.data(), buffer.size());
                buffer.clear();
        }
        blindRotateOptNtt(out, bskMPOptServer.bskFirst, bskMPOptServer.bskDft, sTlwe, v, param);
    }
    std::filesystem::remove("bsk_serialized_GINX_opt_bsk.bin");
}

BENCHMARK_DEFINE_F(BlindRotateBenchmark, OURS_NAIVE)(benchmark::State& state) {
    BootstrappingKeyMPOpt bskMPLazyServer{param, param.lApprox, true, true};
    for (auto _ : state) {
        {
            {
                bskMPLazyServer.bskFirst = bskMPLazyOpt.bskFirst;
                for (auto i = 0; i < param.n - 1; i++) {
                    bskMPLazyServer.bskDft[i] = bskMPLazyOpt.bskTrim[i];
                }
            }
            blindRotateLazyMTNtt(out, bskMPLazyServer.bskFirst, bskMPLazyServer.bskDft, bskMPLazyOpt.bskDecompA, sTlwe, v, s2Dft, param);
            {
                vector<char> buffer;
                std::ostringstream oss1(std::ios::binary);
                serialize(bskMPLazyServer.bskFirst[0], oss1);
                const std::string &data1 = oss1.str();
                buffer.insert(buffer.end(), data1.begin(), data1.end());
                for (auto i = 0; i < param.n - 1; i++) {
                    std::ostringstream oss(std::ios::binary);
                    serialize(bskMPLazyServer.bskDft[i][0], oss);
                    const std::string &data = oss.str();
                    buffer.insert(buffer.end(), data.begin(), data.end());
                }
                std::ofstream outFile2("bsk_serialized_LAZY_MT_bsk.bin", std::ios::binary | std::ios::app);
                outFile2.write(buffer.data(), buffer.size());
                buffer.clear();
            }
            bskMPLazyServer.initialized = true;
        }
    }
    std::filesystem::remove("bsk_serialized_LAZY_MT_bsk.bin");
}

BENCHMARK_DEFINE_F(BlindRotateBenchmark, OURS_PIPELINE)(benchmark::State& state) {
    BootstrappingKeyMPOpt bskMPLazyServer{param, param.lApprox, true, true};
    for (auto _ : state) {
        {
            {
                bskMPLazyServer.bskFirst = bskMPLazyPipe.bskFirst;
                bskMPLazyServer.bskDft[0] = bskMPLazyPipe.bskFull[0];
                bskMPLazyServer.bskDft[1] = bskMPLazyPipe.bskFull[1];
                for (auto i = 2; i < param.n - 1; i++) {
                    bskMPLazyServer.bskDft[i] = bskMPLazyPipe.bskTrim[i - 2];
                }
            }
            blindRotateLazyPipeSerializationNtt(out, bskMPLazyServer.bskFirst, bskMPLazyServer.bskDft, bskMPLazyPipe.bskDecompA, sTlwe, v, s2Dft, one, param);
            bskMPLazyServer.initialized = true;
        }
    }
    std::filesystem::remove("bsk_serialized_PIPE.bin");
}

BENCHMARK_REGISTER_F(BlindRotateBenchmark, GINX_OPT)
    ->Unit(benchmark::kMicrosecond)
    ->Iterations(500);
BENCHMARK_REGISTER_F(BlindRotateBenchmark, OURS_NAIVE)
    ->Unit(benchmark::kMicrosecond)
    ->Iterations(500);
BENCHMARK_REGISTER_F(BlindRotateBenchmark, OURS_PIPELINE)
    ->Unit(benchmark::kMicrosecond)
    ->Iterations(500);

BENCHMARK_MAIN();