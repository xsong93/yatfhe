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

     void readGinx(BootstrappingKeyMP& bsk, const int id) {
        const std::string filename = generateGinxKeyFilename(id);
        printMsg(filename, "Loading GINX key from disk");

        deserializeBskMP(bsk, filename, fixedN);
    }

    void readLazy(BootstrappingKeyMPLazyPipeAlt&bsk, const int id) {
        std::string filename = generateLazyKeyFilename(id);
        printMsg(filename, "Loading Lazy key from disk");

        deserializeBskLazyPipeAlt(bsk, filename, fixedN);
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
    static constexpr size_t DEFAULT_CAPACITY{25};

    SimpleLRUCache<int, BootstrappingKeyMP> ginxCache;
    SimpleLRUCache<int, BootstrappingKeyMPLazyPipeAlt> lazyCache;
    DiskReader diskReader;

    size_t ginxRequest{};
    size_t lazyRequest{};
    size_t ginxHits{};
    size_t lazyHits{};

public:
    explicit SimpleCacheManager(const int n, const size_t capacity = DEFAULT_CAPACITY)
        : ginxCache(capacity, [this](BootstrappingKeyMP& result, const int id) {
            diskReader.readGinx(result, id);
        })
        , lazyCache(capacity, [this](BootstrappingKeyMPLazyPipeAlt& result, const int id) {
            diskReader.readLazy(result, id);
        })
        , diskReader(n) {}

    SimpleCacheManager(const int n, const size_t ginxCapacity, const size_t lazyCapacity)
            : ginxCache(ginxCapacity, [this](BootstrappingKeyMP& result, const int id) {
        diskReader.readGinx(result, id);
    })
            , lazyCache(lazyCapacity, [this](BootstrappingKeyMPLazyPipeAlt& result, const int id) {
                diskReader.readLazy(result, id);
            })
            , diskReader(n) {}

    BootstrappingKeyMP& getGinxKey(const int id) {
        ++ginxRequest;

        auto [fst, snd] = ginxCache.get(id);
        if (snd) {
            ++ginxHits;
        }
        return fst;
    }

    BootstrappingKeyMPLazyPipeAlt& getLazyKey(const int id) {
        ++lazyRequest;

        auto [fst, snd] = lazyCache.get(id);
        if (snd) {
            ++lazyHits;
        }
        return fst;
    }

    BootstrappingKeyMPLazyPipeAlt* getLazyKeySimple(const int id) {
        ++lazyRequest;

        auto* result = lazyCache.getSimple(id);
        if (result != nullptr) {
            ++lazyHits;
        }
        return result;
    }

    void preloadGinxKeys(const std::vector<int>& ids) {
        for (int id : ids) {
            if (id >= 0 && id <= 50) {
                ginxCache.get(id);
            }
        }
    }

    void preloadLazyKeys(const std::vector<int>& ids) {
        for (int id : ids) {
            if (id >= 0 && id <= 50) {
                lazyCache.get(id);
            }
        }
    }

    void putMpKey(const int id, const BootstrappingKeyMP& key) {
        ginxCache.put(id, key);
    }

    void putLazyKey(const int id, const BootstrappingKeyMPLazyPipeAlt& key) {
        lazyCache.put(id, key);
    }

    struct Stats {
        size_t ginxHits;
        size_t ginxRequest;
        size_t lazyHits;
        size_t lazyRequest;
        size_t ginxCacheSize;
        size_t lazyCacheSize;
        size_t mpCacheCapacity;
        size_t lazyCacheCapacity;

        double ginxHitRate() const {
            return ginxHits > 0 ? static_cast<double>(ginxHits) / ginxRequest : 0.0;
        }

        double lazyHitRate() const {
            return lazyHits > 0 ? static_cast<double>(lazyHits) / lazyRequest : 0.0;
        }
    };

    Stats getStats() const {
        return Stats{
                ginxHits,
                ginxRequest,
                lazyHits,
                lazyRequest,
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
        ginxHits = 0;
        ginxRequest = 0;
        lazyHits = 0;
        lazyRequest = 0;
    }

    void resize(const size_t new_capacity) {
        ginxCache = SimpleLRUCache<int, BootstrappingKeyMP>(
            new_capacity, [this](BootstrappingKeyMP& bsk, const int id) {
                diskReader.readGinx(bsk, id);
            });
        lazyCache = SimpleLRUCache<int, BootstrappingKeyMPLazyPipeAlt>(
            new_capacity, [this](BootstrappingKeyMPLazyPipeAlt& bsk, const int id) {
                diskReader.readLazy(bsk, id);
            });
    }
};


#endif //BASE_CACHE_MANAGER_H