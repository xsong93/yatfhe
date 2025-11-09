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
    int fixed_n_;

public:
    explicit DiskReader(const int n) : fixed_n_(n) {}

    BootstrappingKeyMP readGinx(const int id) {
        const std::string filename = generateGinxKeyFilename(id);
        printMsg("Loading BootstrappingKeyMP from disk: ", filename);

        BootstrappingKeyMP bsk;
        deserializeBskMP(bsk, filename, fixed_n_);
        return bsk;
    }

    BootstrappingKeyMPLazyPipeAlt readLazy(const int id) {
        std::string filename = generateLazyKeyFilename(id);
        printMsg("Loading BootstrappingKeyMPLazyPipeAlt from disk: ", filename);

        BootstrappingKeyMPLazyPipeAlt bskLazy;
        deserializeBskLazyPipeAlt(bskLazy, filename, fixed_n_);
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

    SimpleLRUCache<int, BootstrappingKeyMP> mp_cache_;
    SimpleLRUCache<int, BootstrappingKeyMPLazyPipeAlt> lazy_cache_;
    DiskReader disk_reader_;

    size_t total_requests_ = 0;
    size_t mp_hits_ = 0;
    size_t lazy_hits_ = 0;

public:
    explicit SimpleCacheManager(const int n, const size_t capacity = DEFAULT_CAPACITY)
        : mp_cache_(capacity, [this](const int id) { return disk_reader_.readGinx(id); })
        , lazy_cache_(capacity, [this](const int id) { return disk_reader_.readLazy(id); })
        , disk_reader_(n) {}

    BootstrappingKeyMP getGinxKey(int id) {
        ++total_requests_;

        if (id < 0 || id > 50) {
            throw std::out_of_range("ID must be between 0 and 50");
        }

        auto result = mp_cache_.get(id);
        ++mp_hits_;
        return result;
    }

    BootstrappingKeyMPLazyPipeAlt getLazyKey(int id) {
        ++total_requests_;

        if (id < 0 || id > 50) {
            throw std::out_of_range("ID must be between 0 and 50");
        }

        auto result = lazy_cache_.get(id);
        ++lazy_hits_;
        return result;
    }

    void preload_mp_keys(const std::vector<int>& ids) {
        for (int id : ids) {
            if (id >= 0 && id <= 50) {
                mp_cache_.get(id);
            }
        }
    }

    void preload_lazy_keys(const std::vector<int>& ids) {
        for (int id : ids) {
            if (id >= 0 && id <= 50) {
                lazy_cache_.get(id);
            }
        }
    }

    void put_mp_key(int id, const BootstrappingKeyMP& key) {
        mp_cache_.put(id, key);
    }

    void put_lazy_key(int id, const BootstrappingKeyMPLazyPipeAlt& key) {
        lazy_cache_.put(id, key);
    }

    struct Stats {
        size_t total_requests;
        size_t mp_hits;
        size_t lazy_hits;
        size_t mp_cache_size;
        size_t lazy_cache_size;
        size_t mp_cache_capacity;
        size_t lazy_cache_capacity;

        double overall_hit_rate() const {
            return total_requests > 0 ? static_cast<double>(mp_hits + lazy_hits) / total_requests : 0.0;
        }

        double mp_hit_rate() const {
            return mp_hits > 0 ? static_cast<double>(mp_hits) / (mp_hits + (total_requests - mp_hits - lazy_hits)) : 0.0;
        }

        double lazy_hit_rate() const {
            return lazy_hits > 0 ? static_cast<double>(lazy_hits) / (lazy_hits + (total_requests - mp_hits - lazy_hits)) : 0.0;
        }
    };

    Stats get_stats() const {
        return Stats{
            total_requests_,
            mp_hits_,
            lazy_hits_,
            mp_cache_.size(),
            lazy_cache_.size(),
            mp_cache_.capacity(),
            lazy_cache_.capacity()
        };
    }

    std::vector<int> get_mp_cached_keys() const {
        return mp_cache_.get_keys();
    }

    std::vector<int> get_lazy_cached_keys() const {
        return lazy_cache_.get_keys();
    }

    void clear_all() {
        mp_cache_.clear();
        lazy_cache_.clear();
        total_requests_ = 0;
        mp_hits_ = 0;
        lazy_hits_ = 0;
    }

    void resize(size_t new_capacity) {
        mp_cache_ = SimpleLRUCache<int, BootstrappingKeyMP>(
            new_capacity, [this](int id) { return disk_reader_.readGinx(id); });
        lazy_cache_ = SimpleLRUCache<int, BootstrappingKeyMPLazyPipeAlt>(
            new_capacity, [this](int id) { return disk_reader_.readLazy(id); });
    }
};


#endif //BASE_CACHE_MANAGER_H