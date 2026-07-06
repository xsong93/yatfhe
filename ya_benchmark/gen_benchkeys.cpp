
#include "yatfhe/bootstrapping.h"
#include "yatfhe/tlwe.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yautil/initializer.h"
#include "yautil/ya_serializer.h"

int main(){
    YatfheParameters param{};
    initYatfhe(param);

    TrgswKey trgswKey{param};
    TlweKey tlweKey{param.n, param.lweStdDev};
    TorusPolynomial v{param.N};

    // key gen
    TrlweKey& trlweKey = trgswKey.trlweKey;
    TlweKeySwitchingKey ksKey{param};
    genTlweKey(tlweKey);
    genTrlweKey(trlweKey);
    TlweKey tlweKsKey = tlweKey;
    tlweKsKey.sigma = param.rlweStdDev;
    genTlweKeySwitchingKey(ksKey, trlweKey, tlweKsKey, param);
    generateTestPolynomial(v, param.torusBase, 2 * param.N);

    for (int i = 1; i <= 50; i++) {
        BootstrappingKeyWWL24 bskWWL24{param, param.lApprox};
        genBootstrappingKeyWWL24(bskWWL24, trgswKey, tlweKey, param);
        string file;
        file.append("BSK_WWL+24_")
                .append(to_string(i))
                .append(".bin");
        serializeBskWWL24(bskWWL24, file);
    }
    for (int i = 1; i <= 50; i++) {
        BootstrappingKeyMPLazyPipeAlt bskMPLazyPipeAlt{param, param.lApprox, true};
        symEncTrlevWithKeyNtt(bskMPLazyPipeAlt.s2Dft, trgswKey.trlweKey, trgswKey.trlweKey.s, true, param);
        genBootstrappingKeyMPLazyPipeAlt(bskMPLazyPipeAlt, trgswKey, tlweKey, v, param);
        string file;
        file.append("BSK_LAZY_")
                .append(to_string(i))
                .append(".bin");
        serializeBskLazyPipeAlt(bskMPLazyPipeAlt, file);
    }
    for (int i = 1; i <= 50; i++) {
        BootstrappingKeyMP bskMP{param, param.lApprox};
        genBootstrappingKeyMP(bskMP, trgswKey, tlweKey, param);
        string file;
        file.append("BSK_GINX_")
                .append(to_string(i))
                .append(".bin");
        serializeBskMP(bskMP, file);
    }
}