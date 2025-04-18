//
// Created by Xintong Song on 2024/5/9.
//
#include "yatfhe/trgsw.h"
#include "yatfhe/trlgsw.h"
#include "yatfhe/ntt14.h"

using namespace NttNative14;

void trlgswEncryptNtt14(Trlgsw& trlgsw, TrlgswDft14& trlgswDft14, const YatfheParameters& param, const TrgswKey& trgswKey, const Integer mu) {
    Trgsw trgsw {param};
    trgswEncZero(trgsw, param, trgswKey);
    trgswAddInteger(trgsw, mu, 0, param);
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
                applyNttForAB(trlgswDft14.trgswDfts[lvl2].trlweDftSamples[lvl][row],
                              trlgsw.trgsws[lvl2].trlweSamples[lvl][row]);
//                printTrgsw(trlgsw.trgsws[lvl2], "lvl2_" + to_string(lvl2));
            }
        }
    }
}


//todo: debug
void trlgswExternalProduct(Trlwe& output, const Trlgsw& trlgsw, Trlwe& trlweInput, const YatfheParameters& param) {
    const auto k = trlweInput.k;
    const auto level1 = param.l;
    const auto level2 = param.l2;

    DecomposedTrlwe decomposedTrlwe {param};
    DecomposedTrlwe decomposedTrlweRes {param, level2};

    gadgetDecomposeTrlwe(decomposedTrlwe, trlweInput, param);

//#pragma omp parallel for collapse(2) private(out)
//    for (auto lvl2 = 0; lvl2 < level2; lvl2++) {
//        for (auto lvl = 0; lvl < level1; lvl++) {
//            for (auto col = 0; col < k + 1; col++) {
//                auto &curr = (col < k) ? decomposedTrlwe.rlwes[lvl].a[col] : decomposedTrlwe.rlwes[lvl].b;
//                for (auto col2 = 0; col2 < k + 1; col2++) {
//                    auto &out = (col2 < k) ? decomposedTrlweRes.rlwes[lvl2].a[col2] : decomposedTrlweRes.rlwes[lvl2].b;
//                    auto &curr2 = (col2 < k) ? trlgsw.trgsws[lvl2].trlweSamples[lvl][col].a[col2]
//                                             : trlgsw.trgsws[lvl2].trlweSamples[lvl][col].b;
//                    polynomialMulAccNaiveI32(out, curr, curr2);
//                }
//            }
//        }
//    }
    for (auto lvl = 0; lvl < level1; lvl++) {
        for (auto col = 0; col < k + 1; col++) {
            auto &curr = (col < k) ? decomposedTrlwe.rlwes[lvl].a[col] : decomposedTrlwe.rlwes[lvl].b;
            for (auto col2 = 0; col2 < k + 1; col2++) {
                for (auto lvl2 = 0; lvl2 < level2; lvl2++) {
                    auto &out = (col2 < k) ? decomposedTrlweRes.rlwes[lvl2].a[col2] : decomposedTrlweRes.rlwes[lvl2].b;
                    auto &curr2 = (col2 < k) ? trlgsw.trgsws[lvl2].trlweSamples[lvl][col].a[col2]
                                             : trlgsw.trgsws[lvl2].trlweSamples[lvl][col].b;
                    polynomialMulAccNaiveT32(out, curr, curr2);
                }
            }
        }
    }
    recomposeTrlwe(output, decomposedTrlweRes, param);
}

//todo: optimize?
void trlgswExternalProductNtt14(Trlwe& output, const TrlgswDft14& trlgswDft14Input, Trlwe& trlweInput, const YatfheParameters& param) {
    const auto k = trlweInput.k;
    const auto level1 = param.l;
    const auto level2 = param.l2;
    vector<TrlweDft14> trlweDftRes14(level2, TrlweDft14(param.k, param.N));

    DecomposedTrlwe decomposedTrlwe {param};
    DecomposedTrlweDft14 decomposedTrlweDft14 {param, level1};
    DecomposedTrlwe decomposedTrlweIntt {param, level2};

    gadgetDecomposeTrlwe(decomposedTrlwe, trlweInput, param);

    for (auto i = 0; i < decomposedTrlwe.l; i++) {
        applyNttForAB(decomposedTrlweDft14.rlweDfts[i], decomposedTrlwe.rlwes[i]);
    }

////#pragma omp parallel for collapse(2) private(out)
//    for (auto lvl2 = 0; lvl2 < level2; lvl2++) {
//        for (auto lvl = 0; lvl < level1; lvl++) {
//            for (auto col = 0; col < k + 1; col++) {
//                auto &curr = (col < k) ? decomposedTrlweDft14.rlweDfts[lvl].a[col] : decomposedTrlweDft14.rlweDfts[lvl].b;
//                for (auto col2 = 0; col2 < k + 1; col2++) {
//                    auto &out = (col2 < k) ? trlweDftRes14[lvl2].a[col2] : trlweDftRes14[lvl2].b;
//                    auto &curr2 = (col2 < k) ? trlgswDft14Input.trgswDfts[lvl2].trlweDftSamples[lvl][col].a[col2]
//                                             : trlgswDft14Input.trgswDfts[lvl2].trlweDftSamples[lvl][col].b;
//                    calModularInnerProductNtt(out, curr, curr2);
//                }
//            }
//        }
//        applyInttForAB(decomposedTrlweIntt.rlwes[lvl2], trlweDftRes14[lvl2]);
//    }
    for (auto lvl = 0; lvl < level1; lvl++) {
        for (auto col = 0; col < k + 1; col++) {
            auto &curr = (col < k) ? decomposedTrlweDft14.rlweDfts[lvl].a[col] : decomposedTrlweDft14.rlweDfts[lvl].b;
            for (auto col2 = 0; col2 < k + 1; col2++) {
                for (auto lvl2 = 0; lvl2 < level2; lvl2++) {
                    auto &out = (col2 < k) ? trlweDftRes14[lvl2].a[col2] : trlweDftRes14[lvl2].b;
                    auto &curr2 = (col2 < k) ? trlgswDft14Input.trgswDfts[lvl2].trlweDftSamples[lvl][col].a[col2]
                                             : trlgswDft14Input.trgswDfts[lvl2].trlweDftSamples[lvl][col].b;
                    calModularInnerProductNtt(out, curr, curr2);
                }
            }
        }
    }
    for (auto lvl2 = 0; lvl2 < level2; lvl2++) {
        applyInttForAB(decomposedTrlweIntt.rlwes[lvl2], trlweDftRes14[lvl2]);
    }
    recomposeTrlwe(output, decomposedTrlweIntt, param);
}
