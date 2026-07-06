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

    void readWWL24(BootstrappingKeyWWL24&bsk, const int id) {
        std::string filename = generateWWL24KeyFilename(id);
        printMsg(filename, "Loading WWL+24 key from disk");

        deserializeBskWWL24(bsk, filename, fixedN);
    }

    static std::string generateGinxKeyFilename(const int id) {
        return "BSK_GINX_" + std::to_string(id) + ".bin";
    }

    static std::string generateLazyKeyFilename(const int id) {
        return "BSK_LAZY_" + std::to_string(id) + ".bin";
    }
    static std::string generateWWL24KeyFilename(const int id) {
        return "BSK_WWL+24_" + std::to_string(id) + ".bin";
    }
};

class SimpleCacheManager {
private:
    static constexpr size_t DEFAULT_CAPACITY{25};

    SimpleLRUCache<int, BootstrappingKeyMP> ginxCache;
    SimpleLRUCache<int, BootstrappingKeyMPLazyPipeAlt> lazyCache;
    SimpleLRUCache<int, BootstrappingKeyWWL24> wwl24Cache;
    DiskReader diskReader;

    size_t ginxRequest{};
    size_t lazyRequest{};
    size_t wwl24Request{};
    size_t ginxHits{};
    size_t lazyHits{};
    size_t wwl24Hits{};

public:
    explicit SimpleCacheManager(const int n, const size_t capacity = DEFAULT_CAPACITY)
        : ginxCache(capacity, [this](BootstrappingKeyMP& result, const int id) {
            diskReader.readGinx(result, id);
        })
        , lazyCache(capacity, [this](BootstrappingKeyMPLazyPipeAlt& result, const int id) {
            diskReader.readLazy(result, id);
        })
        , wwl24Cache(capacity, [this](BootstrappingKeyWWL24& result, const int id) {
            diskReader.readWWL24(result, id);
        })
        , diskReader(n) {}

    SimpleCacheManager(const int n, const size_t ginxCapacity, const size_t lazyCapacity, const size_t wwl24Capacity)
            : ginxCache(ginxCapacity, [this](BootstrappingKeyMP& result, const int id) {
                diskReader.readGinx(result, id);
            })
            , lazyCache(lazyCapacity, [this](BootstrappingKeyMPLazyPipeAlt& result, const int id) {
                diskReader.readLazy(result, id);
            })
            , wwl24Cache(wwl24Capacity, [this](BootstrappingKeyWWL24& result, const int id) {
                diskReader.readWWL24(result, id);
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

    const BootstrappingKeyMP* getGinxKeySimple(const int id) {
        ++ginxRequest;

        const auto* result = ginxCache.getSimple(id);
        if (result != nullptr) {
            ++ginxHits;
        }
        return result;
    }

    BootstrappingKeyWWL24& getWWL24Key(const int id) {
        ++wwl24Request;

        auto [fst, snd] = wwl24Cache.get(id);
        if (snd) {
            ++wwl24Hits;
        }
        return fst;
    }

    BootstrappingKeyWWL24* getWWL24KeySimple(const int id) {
        ++wwl24Request;

        auto* result = wwl24Cache.getSimple(id);
        if (result != nullptr) {
            ++wwl24Hits;
        }
        return result;
    }

    BootstrappingKeyMPLazyPipeAlt& getLazyKey(const int id) {
        ++lazyRequest;

        auto [fst, snd] = lazyCache.get(id);
        if (snd) {
            ++lazyHits;
        }
        return fst;
    }

    const BootstrappingKeyMPLazyPipeAlt* getLazyKeySimple(const int id) {
        ++lazyRequest;

        const auto* result = lazyCache.getSimple(id);
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

    void putMpKey(const int id, BootstrappingKeyMP&& key) {
        ginxCache.put(id, std::move(key));
    }

    void putWWL24Key(const int id, BootstrappingKeyWWL24&& key) {
        wwl24Cache.put(id, std::move(key));
    }

    void putLazyKey(const int id, BootstrappingKeyMPLazyPipeAlt&& key) {
        lazyCache.put(id, std::move(key));
    }

    struct Stats {
        size_t ginxHits;
        size_t ginxRequest;
        size_t wwl24Hits;
        size_t wwl24Request;
        size_t lazyHits;
        size_t lazyRequest;
        size_t ginxCacheSize;
        size_t wwl24CacheSize;
        size_t lazyCacheSize;
        size_t mpCacheCapacity;
        size_t wwl24CacheCapacity;
        size_t lazyCacheCapacity;

        double ginxHitRate() const {
            return ginxHits > 0 ? static_cast<double>(ginxHits) / ginxRequest : 0.0;
        }

        double wwl24HitRate() const {
            return wwl24Hits > 0 ? static_cast<double>(wwl24Hits) / wwl24Request : 0.0;
        }

        double lazyHitRate() const {
            return lazyHits > 0 ? static_cast<double>(lazyHits) / lazyRequest : 0.0;
        }
    };

    Stats getStats() const {
        return Stats{
                ginxHits,
                ginxRequest,
                wwl24Hits,
                wwl24Request,
                lazyHits,
                lazyRequest,
                ginxCache.size(),
                wwl24Cache.size(),
                lazyCache.size(),
                ginxCache.capacity(),
                wwl24Cache.capacity(),
                lazyCache.capacity()
        };
    }

    std::vector<int> getGinxCachedKeys() const {
        return ginxCache.get_keys();
    }

    std::vector<int> getWWL24CachedKeys() const {
        return wwl24Cache.get_keys();
    }

    std::vector<int> getLazyCachedKeys() const {
        return lazyCache.get_keys();
    }


    void resetStats() {
        ginxHits = 0;
        ginxRequest = 0;
        wwl24Hits = 0;
        wwl24Request = 0;
        lazyHits = 0;
        lazyRequest = 0;
    }

    void clearAll() {
        ginxCache.clear();
        wwl24Cache.clear();
        lazyCache.clear();
        ginxHits = 0;
        ginxRequest = 0;
        wwl24Hits = 0;
        wwl24Request = 0;
        lazyHits = 0;
        lazyRequest = 0;
    }

    void resize(const size_t new_capacity) {
        ginxCache = SimpleLRUCache<int, BootstrappingKeyMP>(
            new_capacity, [this](BootstrappingKeyMP& bsk, const int id) {
                diskReader.readGinx(bsk, id);
        });
        wwl24Cache = SimpleLRUCache<int, BootstrappingKeyWWL24>(
            new_capacity, [this](BootstrappingKeyWWL24& bsk, const int id) {
                diskReader.readWWL24(bsk, id);
        });
        lazyCache = SimpleLRUCache<int, BootstrappingKeyMPLazyPipeAlt>(
            new_capacity, [this](BootstrappingKeyMPLazyPipeAlt& bsk, const int id) {
                diskReader.readLazy(bsk, id);
        });
    }
};


#endif //BASE_CACHE_MANAGER_H