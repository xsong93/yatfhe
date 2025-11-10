//
// Created by xsong93 on 2025/11/9.
//

#ifndef BASE_CACHE_MANAGER_H
#define BASE_CACHE_MANAGER_H

#include "yatfhe/bootstrapping.h"
#include "yautil/lru_cache.h"
#include "yautil/ya_serializer.h"
#include "yautil/tool.h"

class DiskReader {
private:
    int fixedN;

public:
    explicit DiskReader(const int n) : fixedN(n) {}

    BootstrappingKeyMP readGinx(const int id) {
        const std::string filename = generateGinxKeyFilename(id);
        printMsg(filename, "Loading GINX key from disk");

        BootstrappingKeyMP bsk;
        deserializeBskMP(bsk, filename, fixedN);
        return bsk;
    }

    BootstrappingKeyMPLazyPipeAlt readLazy(const int id) {
        std::string filename = generateLazyKeyFilename(id);
        printMsg(filename, "Loading Lazy key from disk");

        BootstrappingKeyMPLazyPipeAlt bskLazy;
        deserializeBskLazyPipeAlt(bskLazy, filename, fixedN);
        return bskLazy;
    }

    static std::string generateGinxKeyFilename(const int id) {
        return "BSK_GINX_" + std::to_string(id) + ".bin";
    }

    static std::string generateLazyKeyFilename(const int id) {
        return "BSK_PIPE_" + std::to_string(id) + ".bin";
    }
};

class SimpleCacheManager {
private:
    static constexpr size_t DEFAULT_CAPACITY = 25;

    SimpleLRUCache<int, BootstrappingKeyMP> ginxCache;
    SimpleLRUCache<int, BootstrappingKeyMPLazyPipeAlt> lazyCache;
    DiskReader diskReader;

    size_t totalRequests = 0;
    size_t ginxHits = 0;
    size_t lazyHits = 0;

public:
    explicit SimpleCacheManager(const int n, const size_t capacity = DEFAULT_CAPACITY)
        : ginxCache(capacity, [this](const int id) { return diskReader.readGinx(id); })
        , lazyCache(capacity, [this](const int id) { return diskReader.readLazy(id); })
        , diskReader(n) {}

    BootstrappingKeyMP getGinxKey(int id) {
        ++totalRequests;

        if (id < 0 || id > 50) {
            throw std::out_of_range("ID must be between 0 and 50");
        }

        auto result = ginxCache.get(id);
        if (result.second) {
            ++ginxHits;
        }
        return result.first;
    }

    BootstrappingKeyMPLazyPipeAlt getLazyKey(int id) {
        ++totalRequests;

        if (id < 0 || id > 50) {
            throw std::out_of_range("ID must be between 0 and 50");
        }

        auto result = lazyCache.get(id);
        if (result.second) {
            ++lazyHits;
        }
        return result.first;
    }

    void preload_mp_keys(const std::vector<int>& ids) {
        for (int id : ids) {
            if (id >= 0 && id <= 50) {
                ginxCache.get(id);
            }
        }
    }

    void preload_lazy_keys(const std::vector<int>& ids) {
        for (int id : ids) {
            if (id >= 0 && id <= 50) {
                lazyCache.get(id);
            }
        }
    }

    void put_mp_key(int id, const BootstrappingKeyMP& key) {
        ginxCache.put(id, key);
    }

    void put_lazy_key(int id, const BootstrappingKeyMPLazyPipeAlt& key) {
        lazyCache.put(id, key);
    }

    struct Stats {
        size_t totalRequests;
        size_t ginxHits;
        size_t lazyHits;
        size_t ginxCacheSize;
        size_t lazyCacheSize;
        size_t mpCacheCapacity;
        size_t lazyCacheCapacity;

        double overallHitRate() const {
            return totalRequests > 0 ? static_cast<double>(ginxHits + lazyHits) / totalRequests : 0.0;
        }

        double ginxHitRate() const {
            return ginxHits > 0 ? static_cast<double>(ginxHits) / (ginxHits + (totalRequests - ginxHits - lazyHits)) : 0.0;
        }

        double lazyHitRate() const {
            return lazyHits > 0 ? static_cast<double>(lazyHits) / (lazyHits + (totalRequests - ginxHits - lazyHits)) : 0.0;
        }
    };

    Stats getStats() const {
        return Stats{
                totalRequests,
                ginxHits,
                lazyHits,
                ginxCache.size(),
                lazyCache.size(),
                ginxCache.capacity(),
                lazyCache.capacity()
        };
    }

    std::vector<int> getGinxCachedKeys() const {
        return ginxCache.get_keys();
    }

    std::vector<int> getLazyCachedKeys() const {
        return lazyCache.get_keys();
    }

    void clearAll() {
        ginxCache.clear();
        lazyCache.clear();
        totalRequests = 0;
        ginxHits = 0;
        lazyHits = 0;
    }

    void resize(size_t new_capacity) {
        ginxCache = SimpleLRUCache<int, BootstrappingKeyMP>(
            new_capacity, [this](int id) { return diskReader.readGinx(id); });
        lazyCache = SimpleLRUCache<int, BootstrappingKeyMPLazyPipeAlt>(
            new_capacity, [this](int id) { return diskReader.readLazy(id); });
    }
};


#endif //BASE_CACHE_MANAGER_H