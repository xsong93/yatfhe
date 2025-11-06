//
// Created by Xintong Song on 2023/12/25.
//

#ifndef HLS_YATFHE_BOOTSTRAPPING_H
#define HLS_YATFHE_BOOTSTRAPPING_H

#include <vector>
#include "yatfhe/tlwe.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/keyswitching.h"

using namespace std;

struct BootstrappingKey {
    vector<TrgswDft> bskDft {};
    vector<Trgsw> bsk {};
    int n {};
    int group {};

    explicit BootstrappingKey(const YatfheParameters& p) : group(p.group) {
        if (group == 1) {
            n = p.n;
            bsk = vector<Trgsw>(n, Trgsw(p));
            bskDft = vector<TrgswDft>(n, TrgswDft(p));
        } else {
            n = p.n / group * (1 << group);
            bsk = vector<Trgsw>(n, Trgsw(p));
            bskDft = vector<TrgswDft>(n, TrgswDft(p));
        }
    }

    BootstrappingKey(const YatfheParameters& p, const int level) : group(p.group) {
        if (group == 1) {
            n = p.n;
            bsk = vector<Trgsw>(n, Trgsw(p, level));
            bskDft = vector<TrgswDft>(n, TrgswDft(p, level));
        } else {
            n = p.n / group * (1 << group);
            bsk = vector<Trgsw>(n, Trgsw(p, level));
            bskDft = vector<TrgswDft>(n, TrgswDft(p, level));
        }
    }
};

struct BootstrappingKeyMP {
    vector<vector<TrgswMPDft>> bskDft {};
    vector<vector<TrgswMP>> bsk {};
    int n {};
    int group {};
    bool isHalf{false};

    explicit BootstrappingKeyMP() = default;

    explicit BootstrappingKeyMP(const YatfheParameters& p) : group(p.group) {
#ifdef TERNARY
        if (group == 1) {
            n = p.n;
        } else {
            n = p.n / group * (1 << group);
        }
        bsk = vector(n, vector(2, TrgswMP(p, p.l)));
        bskDft = vector(n, vector(2, TrgswMPDft(p, p.l)));
#else
        if (group == 1) {
            n = p.n;
        } else {
            n = p.n / group * (1 << group);
        }
        bsk = vector(n, vector(1, TrgswMP(p, p.l)));
        bskDft = vector(n, vector(1, TrgswMPDft(p, p.l)));
#endif
    }

    BootstrappingKeyMP(const YatfheParameters& p, const int level) : group(p.group) {
#ifdef TERNARY
        if (group == 1) {
            n = p.n;
        } else {
            n = p.n / group * (1 << group);
        }
        bsk = vector(n, vector(2, TrgswMP(p, level)));
        bskDft = vector(n, vector(2, TrgswMPDft(p, level)));
#else
        if (group == 1) {
            n = p.n;
        } else {
            n = p.n / group * (1 << group);
        }
        bsk = vector(n, vector(1, TrgswMP(p, level)));
        bskDft = vector(n, vector(1, TrgswMPDft(p, level)));
#endif
    }
    BootstrappingKeyMP(const YatfheParameters& p, const int level, bool half) : group(p.group), isHalf(half) {
#ifdef TERNARY
        if (group == 1) {
            n = p.n;
        } else {
            n = p.n / group * (1 << group);
        }
        bsk = vector(n, vector(2, TrgswMP(p, level)));
        bskDft = vector(n, vector(2, TrgswMPDft(p, level)));
#else
        if (group == 1) {
            n = p.n;
        } else {
            n = p.n / group * (1 << group);
        }
        bsk = vector(n, vector(1, TrgswMP(p, level, isHalf)));
        bskDft = vector(n, vector(1, TrgswMPDft(p, level, isHalf)));
#endif
    }
};

struct BootstrappingKeyMPOpt {
    vector<Trlwe> bskFirst{};
    vector<vector<TrgswMPDft>> bskDft {};
    int n {};
    int group {};

    explicit BootstrappingKeyMPOpt() = default;

    BootstrappingKeyMPOpt(const YatfheParameters& p, const int level, bool half) : group(p.group) {
#ifdef TERNARY
        if (group == 1) {
            n = p.n;
        } else {
            n = p.n / group * (1 << group);
        }
        bskFirst = vector(2, Trlwe{p.k, p.N});
        bskDft = vector(n - 1, vector(2, TrgswMPDft(p, level, isHalf)));
#else
        if (group == 1) {
            n = p.n;
        } else {
            n = p.n / group * (1 << group);
        }
        bskFirst = vector(1, Trlwe{p.k, p.N});
        bskDft = vector(n - 1, vector(1, TrgswMPDft(p, level, half)));
#endif
    }

    BootstrappingKeyMPOpt(const YatfheParameters& p, const int level, bool half, const bool isOnlyB) : group(p.group) {
#ifdef TERNARY
        if (group == 1) {
            n = p.n;
        } else {
            n = p.n / group * (1 << group);
        }
        bskFirst = vector(2, Trlwe{p.k, p.N});
        bskDft = vector(n - 1, vector(2, TrgswMPDft(p, level, isHalf)));
#else
        if (group == 1) {
            n = p.n;
        } else {
            n = p.n / group * (1 << group);
        }
        bskFirst = vector(1, Trlwe{p.k, p.N});
        bskDft = vector(n - 1, vector(1, TrgswMPDft(p, level, half, isOnlyB)));
#endif
    }
};

struct BootstrappingKeyMPLazy {
    vector<Trlwe> bskFirst{};
    vector<vector<TrgswMPDft>> bskTrim{};
    vector<vector<vector<vector<vector<DecompPolynomial>>>>> bskDecompA{};
    int n{};
    int level{};
    int group{};
    bool initialized{false};

    explicit BootstrappingKeyMPLazy() = default;

    BootstrappingKeyMPLazy(const YatfheParameters& p, const int level, const bool isHalf, const bool isOnlyB) : level(level), group(p.group) {
#ifdef TERNARY
        if (group == 1) {
            n = p.n;
        } else {
            n = p.n / group * (1 << group);
        }
        bskFirst = vector(2, Trlwe{p.k, p.N});
        bskTrim = vector(n - 1, vector(2, TrgswMPDft{p, level, isHalf, isOnlyB}));
        bskDecompA = vector(n - 1, vector(2, vector(level, vector(p.l, vector(p.k, DecompPolynomial{p.N})))));
#else
        if (group == 1) {
            n = p.n;
        } else {
            n = p.n / group * (1 << group);
        }
        bskFirst = vector(1, Trlwe{p.k, p.N});
        bskTrim = vector(n - 1, vector(1, TrgswMPDft{p, level, isHalf, isOnlyB}));
        bskDecompA = vector(n - 1, vector(1, vector(level, vector(p.l, vector(p.k, DecompPolynomial{p.N})))));
#endif
    }
};

struct BootstrappingKeyMPLazyPipe {
    vector<Trlwe> bskFirst{};
    vector<vector<TrgswMPDft>> bskDft{};
    vector<vector<vector<vector<vector<DecompPolynomial>>>>> bskDecompA{};
    int n{};
    int level{};
    int group{};
    bool initialized{false};

    explicit BootstrappingKeyMPLazyPipe() = default;

    BootstrappingKeyMPLazyPipe(const YatfheParameters& p, const int level, const bool isHalf, const bool isOnlyB) : level(level), group(p.group) {
#ifdef TERNARY
        if (group == 1) {
            n = p.n;
        } else {
            n = p.n / group * (1 << group);
        }
        bskFirst = vector(2, Trlwe{p.k, p.N});
        bskFull = vector(2, vector(2, TrgswMPDft{p, level}));
        bskTrim = vector(n - 3, vector(2, TrgswMPDft{p, level, isHalf, isOnlyB}));
        bskDecompA = vector(n - 3, vector(2, vector(level, vector(p.l, vector(p.k, DecompPolynomial{p.N})))));
#else
        if (group == 1) {
            n = p.n;
        } else {
            n = p.n / group * (1 << group);
        }
        bskFirst = vector(1, Trlwe{p.k, p.N});
        auto bskFull = vector(2, vector(1, TrgswMPDft{p, level}));
        auto bskTrim = vector(n - 3, vector(1, TrgswMPDft{p, level, isHalf, isOnlyB}));
        bskDft.reserve(bskFull.size() + bskTrim.size());
        bskDft.insert(bskDft.end(), bskFull.begin(), bskFull.end());
        bskDft.insert(bskDft.end(), bskTrim.begin(), bskTrim.end());
        bskDecompA = vector(n - 3, vector(1, vector(level, vector(p.l, vector(p.k, DecompPolynomial{p.N})))));
#endif
    }
};

struct BootstrappingKeyMPLazyPipeAlt {
    vector<Trlwe> bskFirst{};
    vector<TrgswMPDft> bskSecond;
    vector<vector<TrgswMP>> bskPrime;
    int n{};
    int level{};
    int group{};
    bool initialized{false};

    explicit BootstrappingKeyMPLazyPipeAlt() = default;

    BootstrappingKeyMPLazyPipeAlt(const YatfheParameters& p, const int level, const bool isHalf) : level(level), group(p.group) {
#ifdef TERNARY
        if (group == 1) {
            n = p.n;
        } else {
            n = p.n / group * (1 << group);
        }
        bskFirst = vector(2, Trlwe{p.k, p.N});
        bskSecond = vector(2, TrgswMPDft{p, level});
        bskPrime = vector(p.n - 2, vector(2, TrgswMP{p, level, isHalf}));
#else
        if (group == 1) {
            n = p.n;
        } else {
            n = p.n / group * (1 << group);
        }
        bskFirst = vector(1, Trlwe{p.k, p.N});
        bskSecond = vector(1, TrgswMPDft{p, level});
        bskPrime = vector(p.n - 2, vector(1, TrgswMP{p, level, isHalf}));
#endif
    }
};

struct BootstrappingKeyMPPreRot {
    vector<Trlwe> bskFirst{};
    vector<vector<TrgswMPDft>> bskDft{};
    int l{};
    int n{};

    explicit BootstrappingKeyMPPreRot() = default;

    explicit BootstrappingKeyMPPreRot(const YatfheParameters& p) :
#ifdef TERNARY
            bskFirst(3, Trlwe{p.k, p.N}),
            bskDft(p.n - 1, vector(3, TrgswMPDft(p, p.l))),
#else
            bskFirst(2, Trlwe{p.k, p.N}),
            bskDft(p.n - 1, vector(2, TrgswMPDft(p, p.l))),
#endif
            l{p.l},
            n {p.n} {}

    BootstrappingKeyMPPreRot(const YatfheParameters& p, const int l) :
#ifdef TERNARY
            bskFirst(3, Trlwe{p.k, p.N}),
            bskDft(p.n - 1, vector(3, TrgswMPDft(p, l))),
#else
            bskFirst(2, Trlwe{p.k, p.N}),
            bskDft(p.n - 1, vector(2, TrgswMPDft(p, l))),
#endif
            l{l},
            n {p.n} {}
};

struct BootstrappingKeyInternal {
    vector<vector<TrgswMP>> bsk {};
    vector<vector<TrgswMPDft>> bskDft {};
    int n {};

    explicit BootstrappingKeyInternal(const YatfheParameters& p) :
            bsk(p.n/2, vector<TrgswMP>(2, TrgswMP(p))),
            bskDft(p.n/2, vector<TrgswMPDft>(2, TrgswMPDft(p))),
            n {p.n} {}
};

struct BootstrappingKeyCRT {
    std::vector<std::vector<TrgswDft24>> bskCRT {}; // n * d
    std::vector<std::vector<Trgsw8>> bsk8 {}; // n * d
    int n {};
    int d {};

    explicit BootstrappingKeyCRT(const YatfheParameters& param) :
            n(param.n),
            bsk8(param.n, std::vector<Trgsw8>(param.d, Trgsw8(param.dh, param.k, param.N))),
            bskCRT(param.n, std::vector<TrgswDft24>(param.d, TrgswDft24(param.dh, param.k, param.N))) {};
};

void functionalBootstrapping(Tlwe& out, const Tlwe& input, const BootstrappingKey& bsk,
                             const TlweKeySwitchingKey& ksk, const TorusPolynomial& v, const YatfheParameters& param);

void functionalBootstrappingNtt(Tlwe& out, const Tlwe& input, const BootstrappingKey& bsk,
                                const TlweKeySwitchingKey& ksk, const TorusPolynomial& v, const YatfheParameters& param);

void functionalBootstrappingCrt(Tlwe& out, const Tlwe& input, const BootstrappingKeyCRT& bskCRT,
                                const TlweKeySwitchingKey& ksk, const TorusPolynomial& v, const YatfheParameters& param);

void genBootstrappingKeyApproxCrt(BootstrappingKeyCRT& bskCRT, TrgswKey& trgswKey, const TlweKey& tlweKey,
                                  const YatfheParameters& param);

void genBootstrappingKey(BootstrappingKey& bsk, TrgswKey& trgswKey, const TlweKey& tlweKey, const YatfheParameters& param);

void genBootstrappingKeyMP(BootstrappingKeyMP& bsk, TrgswKey& trgswKey, const TlweKey& tlweKey,
                           const YatfheParameters& param);

void genBootstrappingKeyMPLazy(BootstrappingKeyMPLazy& bsk, TrgswKey& trgswKey, const TlweKey& tlweKey,
                               const TorusPolynomial& v, const YatfheParameters& param);

void genBootstrappingKeyMPOpt(BootstrappingKeyMPOpt& bsk, TrgswKey& trgswKey, const TlweKey& tlweKey,
                              const TorusPolynomial& v, const YatfheParameters& param);

void genBootstrappingKeyMPLazyPipe(BootstrappingKeyMPLazyPipe& bsk, TrgswKey& trgswKey, const TlweKey& tlweKey,
                                   const TorusPolynomial& v, const YatfheParameters& param);

void genBootstrappingKeyMPLazyPipeAlt(BootstrappingKeyMPLazyPipeAlt& bsk, const TrgswKey& trgswKey, const TlweKey& tlweKey,
                                      const TorusPolynomial& v, const YatfheParameters& param);

void genBootstrappingKeyMPFixedNoise(BootstrappingKeyMP& bsk, TrgswKey& trgswKey, const TlweKey& tlweKey, Torus noise,
                                     const YatfheParameters& param);

void genBootstrappingKeyMPPreRot(BootstrappingKeyMPPreRot& bsk, const TrgswKey& trgswKey, const TlweKey& tlweKey,
                                 const TorusPolynomial& v, int batchSize, const YatfheParameters& param);

void genBootstrappingKeyMPPreRotTernaryFixedNoise(BootstrappingKeyMPPreRot& bsk, const TrgswKey& trgswKey, const TlweKey& tlweKey,
                                                  const TorusPolynomial& v, int batchSize, Torus noise, const YatfheParameters& param);

void decompBootstrappingKeyMcrt(BootstrappingKeyCRT& bskCRT, const BootstrappingKey& bsk, const YatfheParameters& param);

#endif //HLS_YATFHE_BOOTSTRAPPING_H
