#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/tlwe.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/trlgsw.h"
#include "yatfhe/numeric.h"
#include "yatfhe/ntt.h"
#include "yatfhe/ntt24.h"
#include "yautil/tool.h"
#include "yautil/initializer.h"
#include "yautil/time_counter.h"
#include "yatfhe/ntt_hexl.h"

int main() {
    YatfheParameters param{};
    initYatfhe(param);

// key gen
    TrgswKey trgswKey{param};
    TrlweKey &trlweKey = trgswKey.trlweKey;
    genTrlweKey(trlweKey);

// enc s^2
    Trlev s2p(param, param.l);
    TrlevDft s2pDft(param, param.l);
    symEncTrlevWithKey(s2p, trlweKey, trlweKey.s, true, param);
    symEncTrlevWithKeyNtt(s2pDft, trlweKey, trlweKey.s, true, param);

// trgsw enc
    TrgswMPDft in1Dft{param, param.lApprox};
    TrgswMP in1{param, param.lApprox, true};
    Integer mu1 = 1;
    encryptTrgswMP(in1, mu1, trgswKey, 0, param);
    {
        // encryptTrgswMPNtt samples cPrime's "a" natively/uniformly over the wide qNtt domain,
        // which switchTrlweToSecretEmbeddingNtt cannot correctly INTT+decompose (it needs a
        // genuinely bounded Torus-domain "a" forward-transformed into NTT domain instead).
        TorusPolynomial muPolyFix{param.N};
        Trlwe scratch{param.k, param.N};
        for (auto lvl = 0; lvl < in1Dft.l; lvl++) {
            muPolyFix.coeffs[0] = static_cast<Torus>(mu1) << (param.torusBits - (lvl + 1) * param.radixBits);
            symEncTrlweMultiSampleNtt(scratch, in1Dft.cPrime[lvl], trgswKey.trlweKey, muPolyFix.coeffs);
        }
    }

    IntPolynomial dec{param.N};
    decryptTrgswMPNtt(dec, in1Dft, param, trgswKey, false);
    printArray(dec.coeffs, "mu1");

// init data
    TrgswMPDft tmpMp{param, param.lApprox};
    tmpMp.cPrime = in1Dft.cPrime;
    TrgswMPDft tmpMp2{param, param.lApprox};
    vector decompA(param.lApprox, vector(param.l, vector(param.k, DecompPolynomial{param.N})));
    for (auto l = 0; l < param.lApprox; l++) {
        NttHexl::applyNtt(tmpMp2.cPrime[l].b, in1.cPrime[l].b);
    }
    for (auto l0 = 0; l0 < param.lApprox; l0++) {
        auto &a = in1.cPrime[l0].a;
        for (auto k = 0; k < param.k; k++) {
            for (auto j = 0; j < param.N; j++) {
                DecomposedData d{param.l};
                gadgetDecompose(d, a[k].coeffs[j], param);
                for (auto l = 0; l < param.l; l++) {
                    decompA[l0][l][k].coeffs[j] = d.value[l] * d.sign;
                }
            }
        }
    }
    TrgswMPDft tmpMp3{param, param.lApprox};
    TrgswMP t1{param, param.lApprox};
    t1.cPrime = in1.cPrime;

// scheme switching
    for (auto l = 0; l < param.lApprox; l++) {
        COUNT_TIME("switchTrlweToSecretEmbedding", switchTrlweToSecretEmbedding(t1.c[l], in1.cPrime[l], s2p, param);)
        COUNT_TIME("switchTrlweToSecretEmbeddingNtt",
                   switchTrlweToSecretEmbeddingNtt(tmpMp.c[l], tmpMp.cPrime[l], s2pDft, param);)
        COUNT_TIME("switchTrlweToSecretEmbeddingNttOpt",
                   switchTrlweToSecretEmbeddingNttOpt(tmpMp2.c[l], tmpMp2.cPrime[l], decompA[l], s2pDft, param);)
        COUNT_TIME("switchTrlweToSecretEmbeddingNttMix",
                   switchTrlweToSecretEmbeddingNttMix(tmpMp3.c[l], tmpMp3.cPrime[l], decompA[l], in1.cPrime[l].b,
                                                      s2pDft, param);)
    }

// result validation
    Trlwe trlwe{param.k, param.N};
    Trlwe res{param.k, param.N};
    TorusPolynomial mu{param.N};
    IntPolynomial plainMult{param.N};
    symEncTrlweSingleSample(trlwe, trlweKey, modSwitchToTorus32(3, param.torusBase), 0);

// switchTrlweToSecretEmbedding
    externalProductTrgswMP(res, t1, trlwe, param.lApprox, param);
    symDecTrlweToInt(dec, res, trlweKey, param.torusBase);
    printArray(dec.coeffs, "dec");
    cout << "switchTrlweToSecretEmbedding: " << (dec.coeffs[0] == 3) << endl;

// switchTrlweToSecretEmbeddingNtt
    COUNT_TIME("externalProductTrgswMPNtt", externalProductTrgswMPNtt(res, tmpMp, trlwe, param.lApprox, param);)
    symDecTrlweToInt(dec, res, trlweKey, param.torusBase);
    printArray(dec.coeffs, "dec");
    cout << "switchTrlweToSecretEmbeddingNtt: "  << (dec.coeffs[0] == 3) << endl;

// switchTrlweToSecretEmbeddingNttOpt
    externalProductTrgswMPNtt(res, tmpMp2, trlwe, param.lApprox, param);
    symDecTrlweToInt(dec, res, trlweKey, param.torusBase);
    printArray(dec.coeffs, "dec");
    cout << "switchTrlweToSecretEmbeddingNttOpt: "  << (dec.coeffs[0] == 3) << endl;

// switchTrlweToSecretEmbeddingNttMix
    externalProductTrgswMPNtt(res, tmpMp3, trlwe, param.lApprox, param);
    symDecTrlweToInt(dec, res, trlweKey, param.torusBase);
    printArray(dec.coeffs, "dec");
    cout << "switchTrlweToSecretEmbeddingNttMix: "  << (dec.coeffs[0] == 3) << endl;
}