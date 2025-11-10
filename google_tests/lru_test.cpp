//
// Created by Xintong Song on 2024/3/19.
//
#include "gtest/gtest.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/numeric.h"
#include "yautil/tool.h"
#include "yautil/cache_manager.h"
#include "yautil/initializer.h"
#include "yautil/zipf.h"
#include "yautil/cache_workload_generator.h"

TEST(LRU_TEST, LRU) {
    YatfheParameters param{};
    initYatfhe(param);

    int cap = 5;
    int maxKeys = 50;

    bool isWrite = false;
    if (isWrite) {
        for (auto i = 1; i <= maxKeys; ++i) {
            // key gen
            TrgswKey trgswKey{param};
            TlweKey tlweKey{param.n, param.lweStdDev};
            TorusPolynomial v{param.N};
            TrlweKey &trlweKey = trgswKey.trlweKey;
            TlweKeySwitchingKey ksKey{param};
            genTlweKey(tlweKey);
            genTrlweKey(trlweKey);
            TlweKey tlweKsKey = tlweKey;
            tlweKsKey.sigma = param.rlweStdDev;
            genTlweKeySwitchingKey(ksKey, trlweKey, tlweKsKey, param);
            generateTestPolynomial(v, param.torusBase, 2 * param.N);
            BootstrappingKeyMP bskMP{param, param.lApprox};
            genBootstrappingKeyMP(bskMP, trgswKey, tlweKey, param);
            auto file = DiskReader::generateGinxKeyFilename(i);
            serializeBskMP(bskMP, file);
            cout << "write GINX " << to_string(i) << endl;
        }

        for (auto i = 1; i <= maxKeys; ++i) {
            // key gen
            TrgswKey trgswKey{param};
            TlweKey tlweKey{param.n, param.lweStdDev};
            TorusPolynomial v{param.N};
            TrlweKey &trlweKey = trgswKey.trlweKey;
            TlweKeySwitchingKey ksKey{param};
            genTlweKey(tlweKey);
            genTrlweKey(trlweKey);
            TlweKey tlweKsKey = tlweKey;
            tlweKsKey.sigma = param.rlweStdDev;
            genTlweKeySwitchingKey(ksKey, trlweKey, tlweKsKey, param);
            generateTestPolynomial(v, param.torusBase, 2 * param.N);
            BootstrappingKeyMPLazyPipeAlt bskMPLazyPipeAlt{param, param.lApprox, true};
            symEncTrlevWithKeyNtt(bskMPLazyPipeAlt.s2Dft, trgswKey.trlweKey, trgswKey.trlweKey.s, true, param);
            genBootstrappingKeyMPLazyPipeAlt(bskMPLazyPipeAlt, trgswKey, tlweKey, v, param);
            auto file = DiskReader::generateGinxKeyFilename(i);
            serializeBskLazyPipeAlt(bskMPLazyPipeAlt, file);
            cout << "write LAZY " << to_string(i) << endl;
        }
    }

    std::cout << "=== Testing LRU Cache ===" << std::endl;

    SimpleCacheManager cache(param.n, cap);


    std::cout << "\n1. Testing basic functionality..." << std::endl;
    auto key1 = cache.getGinxKey(8);
    assert(key1.n == param.n);

    auto key2 = cache.getLazyKey(10);
    assert(key2.n == param.n);


    std::cout << "\n2. Testing cache hits..." << std::endl;
    auto key1_again = cache.getGinxKey(8);
    auto key2_again = cache.getLazyKey(10);

    auto stats = cache.getStats();
    std::cout << "Total requests: " << stats.totalRequests << std::endl;
    std::cout << "GINX hits: " << stats.ginxHits << std::endl;
    std::cout << "Lazy hits: " << stats.lazyHits << std::endl;
    std::cout << "Overall hit rate: " << stats.overallHitRate() * 100 << "%" << std::endl;


    std::cout << "\n3. Testing LRU eviction..." << std::endl;
    for (int i = 1; i <= 6; ++i) {
        cache.getGinxKey(i);
    }

    auto mp_keys = cache.getGinxCachedKeys();
    std::cout << "Cached MP keys: ";
    for (int key : mp_keys) {
        std::cout << key << " ";
    }
    std::cout << std::endl;


    std::cout << "\n4. Testing preloading..." << std::endl;
    cache.preload_mp_keys({30, 40, 50});
    cache.preload_lazy_keys({35, 45});

    std::cout << "Test completed successfully!" << std::endl;
}

TEST(LRU_TEST, LRU_ZIPF) {
    SimpleCacheManager cache(512, 10);
    CacheWorkloadGenerator workload(1.0);

    auto accessPattern = workload.generateAccessPattern(1000);
    printArray(accessPattern, "access pattern");

    std::cout << "Testing cache with Zipf distribution (s=0.8)" << std::endl;

    for (int id : accessPattern) {
        try {
            cache.getLazyKey(id);
        } catch (const std::exception& e) {
            std::cerr << "Error accessing ID " << id << ": " << e.what() << std::endl;
        }
    }

    auto stats = cache.getStats();
    std::cout << "Total requests: " << stats.totalRequests << std::endl;
    std::cout << "Hit rate: " << stats.lazyHitRate() * 100 << "%" << std::endl;
}