//
// Created by Xintong Song on 2024/5/23.
//
#include "yautil/initializer.h"
#include "yatfhe/ntt14.h"
#include "yatfhe/ntt24.h"
#include "yatfhe/ntt32.h"
#include "yatfhe/ntt64.h"
#include "yatfhe/ntt_hexl.h"
#include "yautil/time_counter.h"
#include "yatfhe/numeric.h"
#include "yautil/multi_threading.h"

int64_t TORUS_Q;
int64_t LWE_Q;
Integer MESSAGE_P;
Integer INT_MAX_VALUE;
Integer INT_MIN_VALUE;
Integer TORUS_MAX;
Integer TORUS_MIN;
Integer LWE_MAX;
Integer LWE_MIN;
NttType NTT_MAX;
NttType NTT_MIN;

void calGadgetVectorW(YatfheParameters& param) {
    auto qCRT = param.qCRT;
    for (size_t i = 0; i < param.dh; i++) {
        auto qHi = param.qh[i];
        param.w[i] = qCRT / qHi * (modInverse(qCRT / qHi, qHi));
    }
//    printf("w1:%ld, w2:%ld\n", param.w[0], param.w[1]);
}

void calGadgetVectorZ(YatfheParameters& param) {
    auto qCRT = param.qCRT;
    for (size_t i = 0; i < param.d; i++) {
        auto qj = param.qd[i];
        auto qjTilde = qCRT / qj;
        param.z[i] = qjTilde * (modInverse(qjTilde, qj));
    }
//    printf("z1:%ld, z2:%ld, z3:%ld, z4:%ld\n", param.z[0], param.z[1], param.z[2], param.z[3]);
}

void initYatfhe(YatfheParameters& param) {
    TORUS_Q = param.q;
    LWE_Q = param.qLwe;
    MESSAGE_P = param.torusBase;
    INT_MAX_VALUE = static_cast<Integer>((TORUS_Q - 1) >> 1);
    INT_MIN_VALUE = static_cast<Integer>(-(TORUS_Q >> 1));
    TORUS_MAX = INT_MAX_VALUE;
    TORUS_MIN = INT_MIN_VALUE;
    LWE_MAX = static_cast<Integer>((LWE_Q - 1) >> 1);
    LWE_MIN = static_cast<Integer>(-(LWE_Q >> 1));
    NTT_MAX = param.qNtt - 1;
    NTT_MIN = 0;
    NttHexl::initNttHexl(param.N, param.qNtt);
    NttHexl::initNttRotMap(param.N);
    NttHexl::initNttRotMinusOneMap(param.N);
    NttHexl::initNttGadgetRecompMap(param.torusBits, param.radixBits, param.l, param.N);
    NttNative64::initGlobalParamsNtt(param.N);
    NttNative32::initGlobalParamsNtt(param.N);
    NttNative24::initGlobalParamsNtt(param.N);
    NttNative14::initGlobalParamsNtt(param.N);
    ThreadPool::initThreadPool();
    for (size_t d = 0; d < param.dl; d++) {
        auto& dh = param.dh;
        param.taoU[dh + d] = static_cast<int>((modInverse(param.qLow / param.ql[d], param.ql[d])));
        param.taoUInv[dh + d] = static_cast<int>((longModP(param.qLow / param.ql[d], param.ql[d])));
    }
    calGadgetVectorW(param);
    calGadgetVectorZ(param);
    COUNT_TIME("init timer", std::cout << std::endl;)
}