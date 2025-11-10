//
// Created by xsong93 on 2025/11/9.
//

#ifndef BASE_LRU_CACHE_H
#define BASE_LRU_CACHE_H

#include <unordered_map>
#include <list>
#include <functional>
#include <string>
#include <fstream>
#include <iostream>

template<typename KeyType, typename ValueType>
class SimpleLRUCache {
private:
    struct CacheNode {
        KeyType key;
        ValueType value;
        CacheNode(KeyType k, ValueType v) : key(std::move(k)), value(std::move(v)) {}
    };

    using NodeList = std::list<CacheNode>;
    using NodeMap = std::unordered_map<KeyType, typename NodeList::iterator>;

    size_t capacity_;
    NodeList node_list_;
    NodeMap node_map_;
    std::function<ValueType(const KeyType&)> loader_;

public:

    SimpleLRUCache(size_t capacity, std::function<ValueType(const KeyType&)> loader)
        : capacity_(capacity), loader_(std::move(loader)) {}

    std::pair<ValueType, bool> get(const KeyType& key) {
        auto it = node_map_.find(key);
        if (it != node_map_.end()) {
            node_list_.splice(node_list_.begin(), node_list_, it->second);
            return {it->second->value, true};
        }

        ValueType value = loader_(key);

        put(key, value);
        return {value, false};
    }

    void put(const KeyType& key, const ValueType& value) {
        auto it = node_map_.find(key);
        if (it != node_map_.end()) {
            it->second->value = value;
            node_list_.splice(node_list_.begin(), node_list_, it->second);
            return;
        }

        if (node_list_.size() >= capacity_) {
            evict();
        }

        node_list_.emplace_front(key, value);
        node_map_[key] = node_list_.begin();
    }

    bool contains(const KeyType& key) const {
        return node_map_.find(key) != node_map_.end();
    }

    size_t size() const {
        return node_list_.size();
    }

    size_t capacity() const {
        return capacity_;
    }

    void clear() {
        node_list_.clear();
        node_map_.clear();
    }

    std::vector<KeyType> get_keys() const {
        std::vector<KeyType> keys;
        for (const auto& node : node_list_) {
            keys.push_back(node.key);
        }
        return keys;
    }

private:

    void evict() {
        if (node_list_.empty()) return;

        auto last_node = node_list_.end();
        --last_node;

        printMsg(last_node->key, "Evicting key");
        node_map_.erase(last_node->key);
        node_list_.pop_back();
    }
};

#endif //BASE_LRU_CACHE_H