//
// Created by ic on 24-4-30.
//
#include "gtest/gtest.h"
#include "yatfhe/ntt64.h"
#include "yatfhe/polynomial.h"
#include "yatfhe/numeric_functions.h"
#include "yautil/time_counter.h"
#include <iostream>
#include <gmp.h>
#include "yautil/tool.h"
#include "yatfhe/ntt.h"
using namespace std;


//std::vector<NttType> extractValues64(const std::string& input) {
//    std::vector<NttType> result;
//    std::istringstream iss(input);
//    std::string token;
//
//    while (std::getline(iss, token, ' ')) {
//        size_t pos = token.find(':');
//        if (pos != std::string::npos) {
//            std::string valueStr = token.substr(pos + 1);
//            int value = std::stoi(valueStr);
//            result.push_back(value);
//        }
//    }
//
//    return result;
//}
//std::vector<int32_t> extractValues32(const std::string& input) {
//    std::vector<int32_t> result;
//    std::istringstream iss(input);
//    std::string token;
//
//    while (std::getline(iss, token, ' ')) {
//        size_t pos = token.find(':');
//        if (pos != std::string::npos) {
//            std::string valueStr = token.substr(pos + 1);
//            int value = std::stoi(valueStr);
//            result.push_back(value);
//        }
//    }
//
//    return result;
//}

TEST(ntt64_test, ntt64_test){
    int N = 4096;
    int depth = clog2(N);
    TwParam nwc_tw(depth);
    TwParam nwc_itw(depth);
    TwRom tw_rom(N);
    genTW_ROM64(tw_rom);
    genNWCparam64(nwc_tw, N, tw_rom, STR_NTT);
    genNWCparam64(nwc_itw, N, tw_rom, STR_INTT);
    IntPolynomial a{N}, b{N}, res{N}, ref{N}, test{N};
    Ntt64Polynomial a_ntt{N}, b_ntt{N}, mul_ntt{N};
    for (int i = 0; i < N; i++) {
        a.coeffs[i] = i - (N>>1);
        b.coeffs[i] = i + (1<<24);
    }
    polynomialMulNaiveModQ(ref, a, b, 1l<<32);
    doNTT64(a_ntt,a,nwc_tw);
    doNTT64(b_ntt,b,nwc_tw);
    for (int i = 0; i < N; i++) {
        mul_ntt.coeffs[i] = modMULT64(a_ntt.coeffs[i], b_ntt.coeffs[i]);
    }
    doINTT64(res,mul_ntt,nwc_itw);
    for (int i = 0; i < N; i++) {
        EXPECT_EQ(res.coeffs[i], ref.coeffs[i]);
    }
    cout<<"breakpoint"<<endl;
    printBanner("NTT_NWC_64");
}

TEST(ntt64_test, single_test) {
    int N = 1024;
    int depth = clog2(N);
    TwParam nwc_tw(depth);
    TwParam nwc_itw(depth);
    TwRom tw_rom(N);
    genTW_ROM64(tw_rom);
    std::string str_ntt = "NWC-DIT-NR-NNT";
    std::string str_intt = "NWC-DIF-RN-INNT";
    genNWCparam64(nwc_tw,N,tw_rom,str_ntt);
    genNWCparam64(nwc_itw,N,tw_rom,str_intt);
    Ntt64Polynomial a{N}, a_ntt{N}, res{N};
    for (int i = 0; i < N; i++) {
        a.coeffs[i] = i;
    }
    DIT_NR64(a_ntt, a, nwc_tw);
    DIF_RN64(res, a_ntt, nwc_itw);
    cout<<"breakpoint"<<endl;
}


TEST(MODMULT_TEST, test1){
    uint64_t x = 0, y = 0;
    uint64_t res1 = 0, res2 = 0;

    int N = 4096;
    int lvl = clog2(N);
    for (int i = 0 ; i < N; i++ ) {
        x += ((MOD64+1)>>lvl);
        y = x + 1;
        res1 = modMULT64(x,y);
        res2 = modMul(x,y);
        EXPECT_EQ(res1, res2);
    }
    printBanner("modmult test");

}
TEST(MULT_TEST,rand_test){
    uint64_t x = 0, y = 0;
    uint64_t res1 = 0, res2 = 0, res3 = 0;
    int N = 100000000;
    for (int i = 0; i < N; i++) {
        x = genUInt64UniformDist(0,(MOD64 - 1));
        y = genUInt64UniformDist(0,(MOD64 - 1));
        res1 = modMULT64(x,y);
        res2 = fastmm(x,y);
        res3 = fastmm_opt(x,y);
        if (res1 != res3) {
            cout<<"ref = "<<res1<<" res = "<<res3<<endl;
        }
        EXPECT_EQ(res1, res3);
    }
}
TEST(MULT_TEST,big_rand_test){
    uint64_t x = 0, y = 0;
    uint64_t res1 = 0, res2 = 0, res3 = 0;
    int N = 10000000;
    for (int i = 0; i < N; i++) {
        x = genUInt64UniformDist((MOD64 + 1)>>1,(MOD64 - 1));
        y = genUInt64UniformDist((MOD64 + 1)>>1,(MOD64 - 1));
        res1 = modMULT64(x,y);
        res2 = fastmm(x,y);
        res3 = fastmm_opt(x,y);
        EXPECT_EQ(res1, res3);
    }
}
TEST(MULT_TEST,little_rand_test){
    uint64_t x = 0, y = 0;
    uint64_t res1 = 0, res2 = 0;
    int N = 10000000;
    for (int i = 0; i < N; i++) {
        x = genUInt64UniformDist(0,N);
        y = genUInt64UniformDist(0,N);
        res1 = modMULT64(x,y);
        res2 = fastmm(x,y);
        EXPECT_EQ(res1, res2);
    }
}

TEST(MULT_TEST,single_test){
    uint64_t x = 0, y = 0;
    uint64_t res1 = 0, res2 = 0, res3 = 0;

    x = 18446744069397807105ul;
    y = 1099511627520ul;
    res1 = modMULT64(x,y);
    res2 = fastmm(x,y);
    res3 = fastmm_opt(x,y);
    cout<<"Ref = "<<res1<<" Res = "<<res3<<endl;
    EXPECT_EQ(res1, res3);
}

TEST(MULT_TEST,range_test){
    uint64_t x = 0, y = 0;
    uint64_t res1 = 0, res2 = 0;
    for (int i = 1; i < 1024; i++) {
        x = MOD64 - 1;
        y = MOD64 - i;
        res1 = modMULT64(x,y);
        res2 = fastmm(x,y);
        if (res1 != res2){
            cout<<"error at: "<<i<<endl;
            cout<<"Expect: "<<res1<<", Got: "<<res2<<endl;
        }
    }
}

//TEST(INTT,DEBUG){
//    int N = 32;
//    string a = "[0:5260404028879888055 1:2080200617370937984 2:8635497518897110094 3:17074629399922632159 4:532751275413342706 5:94730927289256948 6:11800853980674729369 7:14181761097009271855 8:6012948908770388828 9:15183787253499871307 10:17178040517662139140 11:16490084431612065324 12:1281149636277345433 13:1027905122032817414 14:11259222992026333351 15:7307492464117241631 16:13480420096253804778 17:15589735580482405852 18:827757031063907394 19:8930242279963162792 20:17041204284390487620 21:10401226624546424232 22:6757777188599747093 23:647659515394691263 24:4322206117477043811 25:15732020634513371785 26:1089572772745239303 27:16700679773495330012 28:190484478297455371 29:7627482078881090742 30:16351998012850196391 31:2986772681141039532 ]";
//    string ref = "[0:374207359 1:909129280 2:136638919 3:-1938377611 4:1948900086 5:-1066868068 6:27031941 7:1292742822 8:-471076183 9:308975578 10:1460128347 11:1915833895 12:1315060084 13:1896488892 14:-1630838405 15:-1132467645 16:703239536 17:380821871 18:308425268 19:-1100670068 20:-836101959 21:821753131 22:537169643 23:-532641213 24:-848097661 25:1221703882 26:-1761980014 27:-2092300082 28:-1353902721 29:-232174243 30:1200564444 31:-855495649 ]";
//    std::vector<NttType> aPrime = extractValues64(a);
//    std::vector<int32_t> ref_vec = extractValues32(ref);
//
//    LagrangePolynomial test_vec{N};
//    for (int i = 0; i < N; i++) {
//        test_vec.coeffs[i] = aPrime[i];
//    }
//    TorusPolynomial resIntt{N};
//    applyIntt(resIntt, test_vec);
//    for (int i = 0; i < N; i++) {
//        if (resIntt.coeffs[i] != ref_vec[i]) {
//            cout<<"error at:"<<i<<endl;
//            cout<<"ref = "<<ref_vec[i]<<" res = "<<resIntt.coeffs[i]<<endl;
//        }
//    }
//
//}

