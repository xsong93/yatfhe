//
// Created by Xintong Song on 2024/4/29.
//
#include "gtest/gtest.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yatfhe/tlwe.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/bootstrapping.h"
#include "yatfhe/numeric_functions.h"
#include "yautil/tool.h"
#include "yautil/initializer.h"

TEST(Cmux, Cmux) {
    YatfheParameters param {};
    param.N = 1024;
    param.k = 2;
    yatfheInit(param);
    int ti = 0;
    while (ti++ < 10) {
        cout << "iter: " << ti << endl;

        // key gen
        TrgswKey trgswKey {param};
        TrlweKey& trlweKey = trgswKey.trlweKey;
        trlweKeyGen(trlweKey);

        // trgsw enc
        Trgsw trgsw {param};
        TrgswDft trgswDft {param};
        Integer mu1 = genIntUniformDist(0, 1);
        trgswEncrypt(trgsw, param, trgswKey, mu1);
        printf( "trgsw dec: %d.\n", trgswDecrypt(trgsw, param, trgswKey));

        // data gen
        Trlwe in2 {param.k, param.N};
        IntPolynomial plain {param.N}; // Z/pZ
        TorusPolynomial plainT {param.N};
        for (auto i = 0; i < plain.N; i++) {
            plain.coeffs[i] = genIntUniformDist(-param.torusBase / 2, param.torusBase / 2 - 1);
            plainT.coeffs[i] = modSwitchToTorus32(plain.coeffs[i], param.torusBase);
        }
        symEncTrlweMultiSample(in2, trlweKey, plainT.coeffs);
        printTrlweAB(in2, "input");

        // pre dec
        IntPolynomial decIn {param.N};
        symDecTrlweToInt(decIn, in2, trlweKey, param.torusBase);
        printArray(decIn.coeffs, "decIn");

        // rotator gen
        int a = genIntUniformDist(0, param.N * 2);;
        cout << "a: " << a << endl;

        // test data gen
        IntPolynomial rotInP {param.N};
        Trlwe rotIn {param.k, param.N};
        if (mu1 == 1) {
            trlweRotate(rotIn, in2, a);
        } else {
            copyTrlwe(rotIn, in2, true, true);
        }
        symDecTrlweToInt(rotInP, rotIn, trlweKey, param.torusBase);
        printArray(rotInP.coeffs, "expect");

        // cmux
        Trlwe out {param.k, param.N};
        controlMux(out, in2, a, trgsw, param);
        printTrlweAB(out, "out");

        // dec
        IntPolynomial decP {param.N};
        symDecTrlweToInt(decP, out, trlweKey, param.torusBase);
        printArray(decP.coeffs, "real");

        //verify
        for (auto i = 0; i < decP.N; i++) {
            ASSERT_EQ(rotInP.coeffs[i], decP.coeffs[i]);
        }
    }
    printBanner("Cmux");
}