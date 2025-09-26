#include "yatfhe/tlwe.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/bootstrapping.h"
#include "yatfhe/blind_rotate.h"
#include "yautil/time_counter.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/numeric_functions.h"
#include "yautil/initializer.h"
#include "yautil/multi_threading.h"

int main(int argc, char **argv) {
    YatfheParameters param{};
    initYatfhe(param);

    // key gen
    TlweKey tlweKey{param.n, param.lweStdDev};
    TrgswKey trgswKey{param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    BootstrappingKey bsKey{param};
    TlweKeySwitchingKey ksKey{param};
    genTlweKey(tlweKey);
    genTrlweKey(trlweKey);
    TlweKey tlweKsKey = tlweKey;
    tlweKsKey.sigma = param.rlweStdDev;
    genTlweKeySwitchingKey(ksKey, trlweKey, tlweKsKey, param);

    TorusPolynomial v {param.N};
    generateTestPolynomial(v, param.torusBase, 2 * param.N);

    BootstrappingKeyInternalAsymOpt bskAsymOpt{param};
    genBootstrappingKeyInternalAsymOpt(bskAsymOpt, trgswKey, tlweKey, v, param);

    TrlevDft s2Dft(param);
    symEncTrlevWithKeyNtt(s2Dft, trlweKey, trlweKey.s, true, param);

    // data gen
    Integer pt = 3;
    cout << "decPre: " << pt << endl;
    Torus mu = modSwitchToTorus32(pt, param.torusBase);
    Tlwe input{param.n};
    symEncTlwe(input, mu, tlweKey);
    ScaledTlwe sTlwe {param.N * 2, param.n};
    rescaleTlweFromTorus32(sTlwe, input);

    Trlwe out{param.k, param.N};
    Tlwe tmp {ksKey.nCurrKey};
    Tlwe output {param.n};

    // rot
    for (auto i = 20; i <= 25; i++) {
        BENCH_CUSTOM("Batch size " + to_string(i), {
            auto bsk = bskAsymOpt;
            blindRotateInternalPairWiseAsymOptNtt(out, bsk.bsk, bsk.bskLast, bsk.bskDft, sTlwe, s2Dft, i, param);
        }, 100)
    }
    extractTlweFromTrlwe(tmp, out, param.driftPhase);
    switchKeyForTlwe(output, ksKey, tmp, param);
    auto decAft = symDecTlweToInt(output, tlweKey, param.torusBase);
    cout << "decAft: "<< decAft << endl;
    cout << "err:" << calTlweError(output, tlweKey, mu) << endl;
    return 0;
}
