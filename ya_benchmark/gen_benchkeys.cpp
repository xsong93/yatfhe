
#include <iostream>

#include "yatfhe/bootstrapping.h"
#include "yatfhe/tlwe.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yautil/initializer.h"
#include "yautil/ya_serializer.h"

int main(int argc, char** argv){
    int users = 50;
    int from = 1;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a.rfind("--users=", 0) == 0) users = std::atoi(a.c_str() + 8);
        if (a.rfind("--from=", 0) == 0)  from  = std::atoi(a.c_str() + 7);
    }
    if (users < 1) users = 1;
    if (from < 1) from = 1;
    std::cout << "generating ids " << from << ".." << users << std::endl;

    YatfheParameters param{};
    initYatfhe(param);

    TrgswKey trgswKey{param};
    TlweKey tlweKey{param.n, param.lweNoiseB};
    TorusPolynomial v{param.N};

    // key gen
    TrlweKey& trlweKey = trgswKey.trlweKey;
    TlweKeySwitchingKey ksKey{param};
    genTlweKey(tlweKey);
    genTrlweKey(trlweKey);
    TlweKey tlweKsKey = tlweKey;
    genTlweKeySwitchingKey(ksKey, trlweKey, tlweKsKey, param);
    generateTestPolynomialFR(v, param.torusBase, 2 * param.N);

    for (int i = from; i <= users; i++) {
        BootstrappingKeyWWL24 bskWWL24{param, param.lApprox};
        genBootstrappingKeyWWL24(bskWWL24, trgswKey, tlweKey, param);
        string file;
        file.append("BSK_WWL+24_")
                .append(to_string(i))
                .append(".bin");
        serializeBskWWL24(bskWWL24, file);
        cout << file << " generated." << endl;
    }
    for (int i = from; i <= users; i++) {
        BootstrappingKeyMPLazyPipeAlt bskMPLazyPipeAlt{param, param.lApprox, true};
        symEncTrlevWithKeyNtt(bskMPLazyPipeAlt.s2Dft, trgswKey.trlweKey, trgswKey.trlweKey.s, true, param);
        genBootstrappingKeyMPLazyPipeAlt(bskMPLazyPipeAlt, trgswKey, tlweKey, param);
        string file;
        file.append("BSK_LAZY_")
                .append(to_string(i))
                .append(".bin");
        serializeBskLazyPipeAlt(bskMPLazyPipeAlt, file);
        cout << file << " generated." << endl;
    }
    for (int i = from; i <= users; i++) {
        BootstrappingKeyMP bskMP{param, param.lApprox};
        genBootstrappingKeyMP(bskMP, trgswKey, tlweKey, param);
        string file;
        file.append("BSK_GINX_")
                .append(to_string(i))
                .append(".bin");
        serializeBskMP(bskMP, file);
        cout << file << " generated." << endl;
    }
}