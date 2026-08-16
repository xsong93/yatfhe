//
// ChaCha20-based CSPRNG.
//
#ifndef YATFHE_CSPRNG_H
#define YATFHE_CSPRNG_H

#include <array>
#include <cstdint>
#include <cstddef>

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

#endif //YATFHE_CSPRNG_H
