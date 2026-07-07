//
// ChaCha20-based CSPRNG.
//
#ifndef HLS_YATFHE_CSPRNG_H
#define HLS_YATFHE_CSPRNG_H

#include <array>
#include <cstdint>
#include <cstddef>

// Seeded once from OS entropy rather than reseeded per call, since ChaCha20's
// keystream is itself cryptographically secure -- unlike mt19937, whose 624-word state can be
// reconstructed from a modest number of its outputs regardless of how it was seeded.
class ChaCha20Rng {
public:
    using result_type = uint32_t;

    ChaCha20Rng();

    result_type operator()();

    static constexpr result_type min() { return 0; }
    static constexpr result_type max() { return UINT32_MAX; }

private:
    void refill();

    std::array<uint32_t, 16> state_{};
    std::array<uint32_t, 16> block_{};
    size_t blockPos_ = 16; // force a refill on first use
};

#endif //HLS_YATFHE_CSPRNG_H
