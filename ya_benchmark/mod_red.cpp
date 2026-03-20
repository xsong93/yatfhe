//
// Created for anonymous review.
//

#include "yatfhe/numeric.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yautil/initializer.h"
#include "yautil/time_counter.h"

uint32_t mod_inverse(uint32_t R, uint32_t m) {
    int32_t t = 0, new_t = 1;
    int32_t r = m, new_r = R;

    while (new_r != 0) {
        int32_t quotient = r / new_r;
        std::swap(t, new_t -= quotient * t);
        std::swap(r, new_r -= quotient * r);
    }

    if (r > 1) return 0; // No modular inverse exists
    return (t < 0) ? t + m : t;
}

int main() {
    YatfheParameters param {};
    initYatfhe(param);

    Ntt64 a = 1234567890123456789ULL;
    Ntt64 b = 9876543210987654321ULL;
    Ntt16 a1 = 59923;
    Ntt16 b1 = 65535;
    Ntt14 a14 = 11280;
    Ntt14 b14 = 11111;
    Ntt24 a24 = 16777216;
    Ntt24 b24 = 16777216;


    cout << longModP(INT64_MIN, TORUS_Q) << endl;
    cout << barrettReduceT32(INT64_MIN) << endl;
    cout << longModP((1l << 32), TORUS_Q) << endl;
    cout << montgomoryReduceT32(INT64_MIN) << endl;
    COUNT_TIME("longModP", for(size_t i = 0; i < 10000; i++){longModP(INT64_MIN, TORUS_Q);})
    COUNT_TIME("barrettReduceT32", for(size_t i = 0; i < 10000; i++){barrettReduceT32(INT64_MIN);})
}