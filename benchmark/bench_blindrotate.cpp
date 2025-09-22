#include "yatfhe/tlwe.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/bootstrapping.h"
#include "yatfhe/blind_rotate.h"
#include "yautil/time_counter.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/numeric_functions.h"
#include "yautil/initializer.h"
#include "yatfhe/ntt_hexl.h"

int main(int argc, char **argv) {
    YatfheParameters param{};
    initYatfhe(param);

    // key gen
    TlweKey tlweKey{param.n, param.lweStdDev};
    genTlweKey(tlweKey);
    TrgswKey trgswKey{param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    genTrlweKey(trlweKey);
    BootstrappingKeyInternal bsk {param};
    genBootstrappingKeyInternal(bsk, trgswKey, tlweKey, param);
    BootstrappingKey bskNor {param};
    genBootstrappingKey(bskNor, trgswKey, tlweKey, param);
    BootstrappingKeyInternalAsym bskAsym {param};
    genBootstrappingKeyInternalAsym(bskAsym, trgswKey, tlweKey, param);

    TorusPolynomial v {param.N};
    generateTestPolynomial(v, param.torusBase, 2 * param.N);

    BootstrappingKeyInternalAsymOpt bskAsymOpt {param};
    genBootstrappingKeyInternalAsymOpt(bskAsymOpt, trgswKey, tlweKey, v, param);

    Trlev s2(param); // (a-s, as+e)
    TrlevDft s2Dft(param);
    encTrlevSingleSample(s2, trlweKey, 0, param);
    for (auto i = 0; i < param.l; i++) {
        for (auto j = 0; j < param.k; j++) {
            subTorusPolynomial(s2.trlwes[i].a[j], s2.trlwes[i].a[j], trlweKey.s[j]);
            NttHexl::applyNtt(s2Dft.trlweDfts[i].a[j], s2.trlwes[i].a[j]);
        }
    }

    // data gen
    Integer pt = 3;
    Torus mu = modSwitchToTorus32(pt, param.torusBase);
    Tlwe input{param.n};
    symEncTlwe(input, mu, tlweKey);
    ScaledTlwe sTlwe {param.N * 2, param.n};
    rescaleTlweFromTorus32(sTlwe, input);
    Trlwe acc{param.k, param.N};
    genNoiselessTrlweSample(acc, v, sTlwe);
    Trlwe out{param.k, param.N};


    // rot
    BENCH100("blindRotateNtt", blindRotateNtt(acc, bskNor.bskDft, sTlwe, param);)
//    COUNT_TIME("blindRotateInternalPireWiseNtt", blindRotateInternalPairWiseNtt(acc, bsk.bsk, bsk.bskDft, sTlwe, param);)
    BENCH100("blindRotateInternalPairWiseAsymNtt", blindRotateInternalPairWiseAsymNtt(acc, bskAsym.bsk, bskAsym.bskDft, sTlwe, s2Dft, param);)
    BENCH100("blindRotateInternalPairWiseAsymOptNtt", blindRotateInternalPairWiseAsymOptNtt(out, bskAsymOpt.bsk, bskAsymOpt.bskLast, bskAsymOpt.bskDft, sTlwe, s2Dft, param);)

    return 0;
}
