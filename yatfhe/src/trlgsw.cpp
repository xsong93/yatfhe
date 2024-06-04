//
// Created by Xintong Song on 2024/5/9.
//
#include "yatfhe/trgsw.h"
#include "yatfhe/trlgsw.h"
#include "yatfhe/ntt14.h"

void trgswEncryptNtt14(Trgsw& trgsw, TrlgswDft14& trlgswDft14, const YatfheParameters& param, const TrgswKey& trgswKey, const Integer mu) {
    trgswEncZero(trgsw, param, trgswKey);
    trgswAddInteger(trgsw, mu, param);
    Trlgsw trlgsw {param};
    for (auto lvl = 0; lvl < param.l; lvl++) {
        for (auto row = 0; row < param.k + 1; row++) {
            DecomposedTrlwe decomposedTrlwe {param, param.l2};
            gadgetDecomposeTrlwe(decomposedTrlwe, trgsw.trlweSamples[lvl][row], param);
            for (auto lvl2 = 0; lvl2 < param.l2; lvl2++) {
                for (auto j = 0; j < param.N; j++) {
                    for (auto col = 0; col < param.k; col++) {
                        trlgsw.trgsws[lvl2].trlweSamples[lvl][row].a[col].coeffs[j] = decomposedTrlwe.rlwes[lvl2].a[col].coeffs[j];
                    }
                    trlgsw.trgsws[lvl2].trlweSamples[lvl][row].b.coeffs[j] = decomposedTrlwe.rlwes[lvl2].b.coeffs[j];
                }
                applyNttForAB14(trlgswDft14.trgswDfts[lvl2].trlweDftSamples[lvl][row], trlgsw.trgsws[lvl2].trlweSamples[lvl][row]);
//                printTrgsw(trlgsw.trgsws[lvl2], "lvl2_" + to_string(lvl2));
            }
        }
    }
}
