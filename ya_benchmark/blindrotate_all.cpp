#include "yatfhe/tlwe.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/bootstrapping.h"
#include "yatfhe/blind_rotate.h"
#include "yautil/time_counter.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/numeric.h"
#include "yautil/initializer.h"

int main(int argc, char **argv) {
    YatfheParameters param{};
    param.N = 1024;
    param.batchSize = 16;
    param.tasksPerThread = 15;
    initYatfhe(param);
    printf("n:%d, k:%d, N:%d, b:%d, l:%d\n", param.n, param.k, param.N, param.radixBits, param.l);

    // key gen
    TlweKey tlweKey{param.n, param.lweStdDev};
    TrgswKey trgswKey{param};
    TrlweKey& trlweKey = trgswKey.trlweKey;
    TlweKeySwitchingKey ksKey{param};
    genTlweKey(tlweKey);
    genTrlweKey(trlweKey);
    TlweKey tlweKsKey = tlweKey;
    tlweKsKey.sigma = param.rlweStdDev;
    genTlweKeySwitchingKey(ksKey, trlweKey, tlweKsKey, param);
    TrlevDft s2Dft{param, param.l};
    symEncTrlevWithKeyNtt(s2Dft, trlweKey, trlweKey.s, true, param);

    BootstrappingKeyMP bskMP{param, param.lApprox};
    genBootstrappingKeyMP(bskMP, trgswKey, tlweKey, param);

    TorusPolynomial v {param.N};
    generateTestPolynomial(v, param.torusBase, 2 * param.N);

    // BootstrappingKeyMPPreRot bskPre{param, param.lApprox};
    // genBootstrappingKeyMPPreRot(bskPre, trgswKey, tlweKey, v, param.batchSize, param);

    BootstrappingKeyMPOpt bskMPOpt{param, param.lApprox, false};
    genBootstrappingKeyMPOpt(bskMPOpt, trgswKey, tlweKey, v, param);

    BootstrappingKeyMPOpt bskMPLazy{param, param.lApprox, true};
    genBootstrappingKeyMPOpt(bskMPLazy, trgswKey, tlweKey, v, param);

    // data gen
    Integer pt = 3;
    cout << "decPre: " << pt << endl;
    Torus mu = modSwitchToTorusGeneral(pt, param.torusBase, LWE_Q);
    Tlwe input{param.n};
    symEncTlwe(input, mu, tlweKey);
    ScaledTlwe sTlwe {param.N * 2, param.n};
    rescaleTlweToNewMod(sTlwe, input);
    Trlwe acc{param.k, param.N};
    genNoiselessTrlweSample(acc, v, sTlwe);

    Trlwe out{param.k, param.N};
    Tlwe tmp{ksKey.nCurrKey};
    Tlwe output {param.n};
    vector<TrgswMPDft> helper(param.n-1, TrgswMPDft{param});
    TrgswMPDft gv{param};
    encryptTrgswMPMultiNtt(gv, v.coeffs, trgswKey, param);

    // rot
    COUNT_TIME("blindRotateJP22Ntt single thread", blindRotateJP22Ntt(acc, bskMP.bskDft, sTlwe, param);)
//    COUNT_TIME("blindRotateWithPreRotNtt single thread", blindRotateWithPreRotNtt(out, bskPre.bskFirst, bskPre.bskDft, sTlwe, param);)
//    COUNT_TIME("blindRotateJP22NttMT multiple threads", blindRotateJP22NttMT(acc, bskMP.bskDft, sTlwe, param);)
//    COUNT_TIME("blindRotateWithPreRotNtt multiple threads", blindRotateWithPreRotNttMT(out, bskPre.bskFirst, bskPre.bskDft, sTlwe, param);)
    COUNT_TIME("blindRotateOptNtt", blindRotateOptNtt(out, bskMPOpt.bskFirst, bskMPOpt.bskDft, sTlwe, v, param);)
    COUNT_TIME("blindRotateLazyNtt", blindRotateLazyNtt(out, bskMPLazy.bskFirst, bskMPLazy.bskDft, bskMPLazy.initialized, sTlwe, v, s2Dft, param);)

    extractTlweFromTrlwe(tmp, out, param.driftPhase);
    switchKeyForTlwe(output, ksKey, tmp, param);
    auto decAft = symDecTlweToInt(output, tlweKey, param.torusBase);
    cout << "decAft: "<< decAft << endl;
    cout << "err:" << calTlweError(output, tlweKey, mu) << endl;
    return 0;
}
