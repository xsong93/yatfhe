#include "yatfhe/tlwe.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/bootstrapping.h"
#include "yatfhe/blind_rotate.h"
#include "yautil/time_counter.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/numeric.h"
#include "yautil/cache_manager.h"
#include "yautil/cache_workload_generator.h"
#include "yautil/initializer.h"
#include "yautil/tool.h"
#include "yautil/ya_serializer.h"

int main(int argc, char **argv) {
    YatfheParameters param{};
    initYatfhe(param);
    printf("n:%d, k:%d, N:%d, b:%d, l:%d\n", param.n, param.k, param.N, param.radixBits, param.l);

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
    TorusPolynomial v {param.N};
    generateTestPolynomial(v, param.torusBase, 2 * param.N);

    // data gen
    Integer pt = 3;
    cout << "decPre: " << pt << endl;
    Torus mu = modSwitchToTorusGeneral(pt, param.torusBase, LWE_Q);
    Tlwe input{param.n};
    symEncTlwe(input, mu, tlweKey);
    ScaledTlwe sTlwe {param.N * 2, param.n};
    rescaleTlweToNewMod(sTlwe, input);
    Trlwe acc{param};
    genNoiselessTrlweSample(acc, v, sTlwe);

    SimpleCacheManager cache(512, 10);
    CacheWorkloadGenerator workload(1.0);
    auto accessPattern = workload.generateAccessPattern(1000);
    printArray(accessPattern, "access pattern");

    std::vector<double> iterationTimesUs;
    for (auto i = 1; i <= 50; i++) {
        std::string file = DiskReader::generateLazyKeyFilename(i);
        auto start = std::chrono::high_resolution_clock::now();

        rescaleTlweToNewMod(sTlwe, input);
        auto* bskServer = cache.getLazyKeySimple(accessPattern[i]);
        Trlwe out{param};
        if (bskServer != nullptr) {
            blindRotateLazyPipeAltNtt(out, bskServer->bskFirst, bskServer->bskPrime,bskServer->s2Dft, sTlwe, v, param);
        } else {
            BootstrappingKeyMPLazyPipeAlt bsk;
            blindRotateLazyPipeAltInitNtt(out, bsk.bskFirst, bsk.bskPrime,bsk.s2Dft, sTlwe,
                v, file, param);
            cache.putLazyKey(accessPattern[i], bsk);
        }

        Tlwe tmp{ksKey.nCurrKey}, output{param.n};
        extractTlweFromTrlwe(tmp, out, param.driftPhase);
        switchKeyForTlwe(output, ksKey, tmp, param);

        auto end = std::chrono::high_resolution_clock::now();
        auto elapsedUs = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

        iterationTimesUs.push_back(static_cast<double>(elapsedUs));
        clearFileCache();
    }
    printArray(iterationTimesUs, "iterationTimesUs");
    auto stats = cache.getStats();
    std::cout << "Total requests: " << stats.totalRequests << std::endl;
    std::cout << "Hit rate: " << stats.lazyHitRate() * 100 << "%" << std::endl;

    return 0;
}
