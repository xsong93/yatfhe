#include "yatfhe/csprng.h"
#include <random>

#if defined(__linux__)
#include <sys/random.h>
#endif

namespace {

constexpr uint32_t kChaChaConstants[4] = {0x61707865, 0x3320646e, 0x79622d32, 0x6b206574};

inline uint32_t rotl(uint32_t x, int n) {
    return (x << n) | (x >> (32 - n));
}

void quarterRound(uint32_t &a, uint32_t &b, uint32_t &c, uint32_t &d) {
    a += b; d ^= a; d = rotl(d, 16);
    c += d; b ^= c; b = rotl(b, 12);
    a += b; d ^= a; d = rotl(d, 8);
    c += d; b ^= c; b = rotl(b, 7);
}

void chacha20Block(const std::array<uint32_t, 16> &in, std::array<uint32_t, 16> &out) {
    auto x = in;
    for (int i = 0; i < 10; i++) {
        quarterRound(x[0], x[4], x[8], x[12]);
        quarterRound(x[1], x[5], x[9], x[13]);
        quarterRound(x[2], x[6], x[10], x[14]);
        quarterRound(x[3], x[7], x[11], x[15]);
        quarterRound(x[0], x[5], x[10], x[15]);
        quarterRound(x[1], x[6], x[11], x[12]);
        quarterRound(x[2], x[7], x[8], x[13]);
        quarterRound(x[3], x[4], x[9], x[14]);
    }
    for (int i = 0; i < 16; i++) {
        out[i] = x[i] + in[i];
    }
}

// Fills `count` 32-bit words with OS-provided entropy: getrandom(2) on Linux (the kernel CSPRNG,
// no /dev/urandom fd exhaustion issues), falling back to std::random_device elsewhere.
void fillFromOsEntropy(uint32_t *words, size_t count) {
#if defined(__linux__)
    size_t filled = 0;
    while (filled < count) {
        ssize_t n = getrandom(words + filled, (count - filled) * sizeof(uint32_t), 0);
        if (n > 0) {
            filled += static_cast<size_t>(n) / sizeof(uint32_t);
        }
    }
#else
    std::random_device rd;
    for (size_t i = 0; i < count; i++) {
        words[i] = rd();
    }
#endif
}

} // namespace

ChaCha20Rng::ChaCha20Rng() {
    state_[0] = kChaChaConstants[0];
    state_[1] = kChaChaConstants[1];
    state_[2] = kChaChaConstants[2];
    state_[3] = kChaChaConstants[3];
    fillFromOsEntropy(&state_[4], 8);  // 256-bit key
    state_[12] = 0;                    // 64-bit block counter
    state_[13] = 0;
    fillFromOsEntropy(&state_[14], 2); // 64-bit nonce
}

void ChaCha20Rng::refill() {
    chacha20Block(state_, block_);
    if (++state_[12] == 0) {
        ++state_[13];
    }
    blockPos_ = 0;
}

ChaCha20Rng::result_type ChaCha20Rng::operator()() {
    if (blockPos_ == 16) {
        refill();
    }
    return block_[blockPos_++];
}
