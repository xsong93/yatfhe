//
// Created by xsong93 on 08/27/2026.
//
#pragma once

#include <cstdlib>
#include <filesystem>
#include <string>

namespace yabench {

inline std::string& benchOutDirRef() {
    static std::string dir = [] {
        const char* env = std::getenv("YA_OUT_DIR");
        return std::string(env ? env : ".");
    }();
    return dir;
}

inline void setBenchOutDir(const std::string& dir) {
    benchOutDirRef() = dir.empty() ? "." : dir;
}

inline std::string benchOutPath(const std::string& filename) {
    const std::string& dir = benchOutDirRef();
    if (dir.empty() || dir == ".") {
        return filename;
    }
    std::error_code ec;
    std::filesystem::create_directories(dir, ec);
    if (ec) {
        return filename;
    }
    return (std::filesystem::path(dir) / filename).string();
}

}  // namespace yabench
