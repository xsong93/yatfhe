//
// Created by ic on 24-4-11.
//
#include "gtest/gtest.h"
#include "yatfhe/myNtt.h"
#include "yatfhe/polynomial.h"
#include "yatfhe/numeric_functions.h"
#include "yautil/time_counter.h"
#include <iostream>
#include <gmp.h>
#include "yautil/tool.h"
using namespace std;

TEST(MYNTT_TEST,test_init) {
    int N = 1024;
    Ntt32_TW TW(N>>1);
    Ntt32_iTW iTW(N>>1);
    std::cout<<MOD<<std::endl;
    genTW(TW);
    geniTW(iTW, TW);
    print_myNtt(TW,iTW);
    NttPolynomial a{N},b {N}, res{N};
    for (int i = 0; i < N; i++) {
        a.coeffs[i] = i;
    }
    std::cout<<"a:N = "<<a.N<<"; b:N = "<<b.N<<std::endl;
    doNTT32(b,a,TW);
    printNttPoly(b);
    doINTT32(res, b, iTW);
    printNttPoly(res);
}

TEST(MYNTT_TEST, modSUB_test){
    Ntt32 a = 0, b = 1;
    std::cout<<modSUB(a,b)<<std::endl;
}
TEST(MYNTT_TEST, modADD_test){
    Ntt32 a = 4293918719, b = 74203740;
    std::cout<<modADD(a,b)<<std::endl;
}

TEST(MYNTT_TEST, mod_inv_test) {
    int N = 4;
    Ntt32_TW TW(N>>1);
    Ntt32_iTW iTW(N>>1);
    std::cout<<MOD<<std::endl;
    genTW(TW);
    geniTW(iTW, TW);
    print_myNtt(TW,iTW);
}

TEST(MYNTT_TEST, modMULT_test) {
    Ntt32 a = 4293918719, b = 37101870;
    std::cout<< modMULT(a,b)<<std::endl;
}

TEST(MYNTT_TEST,construct_test){
    int N = 16;
    ROM rom(N);
    cout<<"ntt_rom : {N, TW_N, PHI_N} = {"<<rom.ntt_rom.N<<", "<<rom.ntt_rom.tw_N<<", "<<rom.ntt_rom.phi_N<<"}"<<endl;
    genROM(rom);
    printROM(rom);
}

TEST(MYNTT_TEST, NWC_NTT_TEST) {
    int N = 64;
    cout<<"MOD = "<<MOD<<endl;
    ROM rom(N);
    DIF_ROM dif_rom(clog2(N));
    genROM(rom);
    genDIF_ROM(dif_rom);
//    printROM(rom);
    NttPolynomial a{N}, b{N}, res{N};
    for (int i = 0; i < N ; i++) {
        a.coeffs[i] = N - i ;
    }
    printNttPoly(a);
//    NWC_NTT32(b,a,rom.ntt_rom);
    DIF_NR(b,a,rom.ntt_rom);
    printNttPoly(b);
//    NWC_INTT32(res, b, rom.intt_rom);
    DIF_RN(res,b,dif_rom.intt_tw);
    printNttPoly(res);
}

TEST(MYNTT_TEST, try_NTT_TEST) {
    int N = 8;
    cout<<"MOD = "<<MOD<<endl;
    ROM rom(N);
    genROM(rom);
    printROM(rom);
    NttPolynomial a{N}, b{N}, res{N};
    for (int i = 0; i < N; i++) {
        a.coeffs[i] = i;
    }
    Ntt32_TW TW(N>>1);
    Ntt32_iTW iTW(N>>1);
    genTW(TW);
    geniTW(iTW, TW);
    print_myNtt(TW,iTW);
}

TEST(MYNTT_TEST, bit_rev_test) {
    int N = 1024;
    std::vector<Ntt32> vec;
    std::cout<<vec.size()<<std::endl;
    for (int i = 0; i < N; i++) {
        vec.push_back(Ntt32(i));
    }
    std::cout<<"vec size is "<<vec.size()<<std::endl;
    std::cout<<"Vector elements:"<<std::endl;
    for (int j = 0; j < N; j++) {
        std::cout<<vec[j]<<" ";
    }
    std::cout<<std::endl;
    bit_rev(vec);
    std::cout<<"Vector reversed elements:"<<std::endl;
    for (int k = 0; k < N; k++) {
        std::cout << vec[k] << " ";
    }
    std::cout<<std::endl;


}

TEST(MYNTT_TEST, dif_rom_test) {
    int N = 32;
    N = clog2(N);
    cout<<"MOD = "<<MOD<<endl;
    DIF_ROM rom(N);
    genDIF_ROM(rom);
    std::cout<<"size of uint "<< sizeof(uint)<<std::endl;
    std::cout<<"size of uint32 "<< sizeof(uint32_t)<<std::endl;
    int32_t a = -3;
    std::cout<<"-3 >> 1 = "<<(a>>1)<<std::endl;
}

TEST(MYNTT_TEST, nwc_rom_test) {
    int N = 32;
    int depth = clog2(N);
    TW_PARAM nwc_tw(depth);
    TW_ROM tw_rom(N);
    genTW_ROM(tw_rom);
    std::string str = "NWC-DIF-RN-INNT";
    genNWCparam(nwc_tw,N,tw_rom,str);

}

TEST(MYNTT_TEST, nwc_ntt_test) {
    int N = 8;
    int depth = clog2(N);
    TW_PARAM nwc_tw(depth);
    TW_PARAM nwc_itw(depth);
    TW_ROM tw_rom(N);
    genTW_ROM(tw_rom);
    std::string str_ntt = "NWC-DIT-NR-NNT";
    std::string str_intt = "NWC-DIF-RN-INNT";
    genNWCparam(nwc_tw,N,tw_rom,str_ntt);
    genNWCparam(nwc_itw,N,tw_rom,str_intt);
    NttPolynomial a{N}, b{N}, res{N};
    for (int i = 0; i < N ; i++) {
        a.coeffs[i] = i ;
    }
    printNttPoly(a);
//    NWC_NTT32(b,a,rom.ntt_rom);
    DIT_NR(b,a,nwc_tw);
    printNttPoly(b);
//    NWC_INTT32(res, b, rom.intt_rom);
    DIF_RN(res,b,nwc_itw);
    printNttPoly(res);


}