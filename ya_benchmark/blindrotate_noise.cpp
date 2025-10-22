#include "yatfhe/tlwe.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/bootstrapping.h"
#include "yatfhe/blind_rotate.h"
#include "yautil/time_counter.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/numeric_functions.h"
#include "yautil/initializer.h"
#include "yautil/tool.h"

int main(int argc, char **argv) {
    YatfheParameters param{};
    param.n = 2;
    // param.N = 512;
    // param.N = 8;
    param.lApprox = 4;
    param.batchSize = 40;
    initYatfhe(param);
    printf("n:%d, k:%d, N:%d, b:%d, l:%d\n", param.n, param.k, param.N, param.radixBits, param.l);

    int64_t noise0 = 1;
    int64_t noise1 = 1;
    int loop = 0;
    while (loop++ < 10000) {
        cout << "loop: " << loop << ", diff: " << noise0-noise1<< endl;
        // key gen
        TlweKey tlweKey{param};
        TrgswKey trgswKey{param};
        TrlweKey& trlweKey = trgswKey.trlweKey;
        TlweKeySwitchingKey ksKey{param};
        genTlweKey(tlweKey);
        tlweKey.s = vector(param.n, 1);
        genTrlweKey(trlweKey);
        TlweKey tlweKsKey = tlweKey;
        tlweKsKey.sigma = param.rlweStdDev;
        genTlweKeySwitchingKey(ksKey, trlweKey, tlweKsKey, param);

        BootstrappingKeyMP bskMP{param};
        genBootstrappingKeyMP(bskMP, trgswKey, tlweKey, param);

        TorusPolynomial v{param.N};
        generateTestPolynomial(v, param.torusBase, 2 * param.N);

        BootstrappingKeyMPPreRot bskPre{param};
        genBootstrappingKeyMPPreRot(bskPre, trgswKey, tlweKey, v, param.batchSize, param);


        // data gen
        Integer pt = 0;
        Torus mu = modSwitchToTorus32(pt, param.torusBase);
        Tlwe input{param.n};
        symEncTlwe(input, mu, tlweKey);
        ScaledTlwe sTlwe {param.N * 2, param.n};
        rescaleTlweToNewMod(sTlwe, input);
        Trlwe acc{param.k, param.N};
        genNoiselessTrlweSample(acc, v, sTlwe);
        Trlwe out{param.k, param.N};
        Tlwe tmp {ksKey.nCurrKey};
        Tlwe output {param.n};
        TorusPolynomial outRlwe{param.N};
        TorusPolynomial outRlwe2{param.N};

        // rot
        blindRotateJP22NttMT(acc, bskMP.bskDft, sTlwe, param);
        blindRotateWithPreRotNttMT(out, bskPre.bskFirst, bskPre.bskDft, sTlwe, param);

        symDecTrlweWoRounding(outRlwe, acc, trlweKey);
        symDecTrlweWoRounding(outRlwe2, out, trlweKey);

        // extractTlweFromTrlwe(tmp, out, param.driftPhase);
        // switchKeyForTlwe(output, ksKey, tmp, param);
        noise0 += std::abs(outRlwe.coeffs[0]);
        //
        // extractTlweFromTrlwe(tmp, acc, param.driftPhase);
        // switchKeyForTlwe(output, ksKey, tmp, param);
        noise1 += std::abs(outRlwe2.coeffs[0]);
    }
    return 0;
}
