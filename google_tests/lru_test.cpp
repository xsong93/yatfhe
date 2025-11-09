//
// Created by Xintong Song on 2024/3/19.
//
#include "gtest/gtest.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/numeric.h"
#include "yautil/tool.h"
#include "yautil/cache_manager.h"
#include "yautil/initializer.h"

TEST(LRU_TEST, LRU) {
    YatfheParameters param{};
    initYatfhe(param);

    int cap = 5;
    int maxKeys = 50;

    for (auto i = 1; i <= maxKeys; ++i) {
        // key gen
        TrgswKey trgswKey{param};
        TlweKey tlweKey{param.n, param.lweStdDev};
        TorusPolynomial v{param.N};
        TrlweKey& trlweKey = trgswKey.trlweKey;
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
        cout <<"write GINX " << to_string(i) << endl;
    }

    for (auto i = 1; i <= maxKeys; ++i) {
        // key gen
        TrgswKey trgswKey{param};
        TlweKey tlweKey{param.n, param.lweStdDev};
        TorusPolynomial v{param.N};
        TrlweKey& trlweKey = trgswKey.trlweKey;
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
        cout <<"write LAZY " << to_string(i) << endl;
    }

    std::cout << "=== Testing LRU Cache ===" << std::endl;

    // 创建缓存管理器
    SimpleCacheManager cache(param.n, cap);

    // 测试1: 基本功能
    std::cout << "\n1. Testing basic functionality..." << std::endl;
    auto key1 = cache.getGinxKey(8);
    assert(key1.n == param.n);

    auto key2 = cache.getLazyKey(10);
    assert(key2.n == param.n);

    // 测试2: 缓存命中
    std::cout << "\n2. Testing cache hits..." << std::endl;
    auto key1_again = cache.getGinxKey(8);
    auto key2_again = cache.getLazyKey(10);

    auto stats = cache.get_stats();
    std::cout << "Total requests: " << stats.total_requests << std::endl;
    std::cout << "GINX hits: " << stats.mp_hits << std::endl;
    std::cout << "Lazy hits: " << stats.lazy_hits << std::endl;
    std::cout << "Overall hit rate: " << stats.overall_hit_rate() * 100 << "%" << std::endl;

    // 测试3: LRU驱逐
    std::cout << "\n3. Testing LRU eviction..." << std::endl;
    // 加载6个不同的键，应该触发驱逐
    for (int i = 1; i <= 6; ++i) {
        cache.getGinxKey(i);
    }

    auto mp_keys = cache.get_mp_cached_keys();
    std::cout << "Cached MP keys: ";
    for (int key : mp_keys) {
        std::cout << key << " ";
    }
    std::cout << std::endl;

    // 测试4: 预加载
    std::cout << "\n4. Testing preloading..." << std::endl;
    cache.preload_mp_keys({30, 40, 50});
    cache.preload_lazy_keys({35, 45});

    std::cout << "Test completed successfully!" << std::endl;
}