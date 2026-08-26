//
// Created by xsong93 on 08/26/2026.
//

#include <chrono>
#include <vector>
#include <iostream>
#include <iomanip>

#include "yatfhe/gadget_decomposition.h"
#include "yatfhe/ntt_hexl.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/yatfhe_parameters.h"
#include "yautil/initializer.h"

using namespace std::chrono;
using namespace NttHexl;

int main() {
    YatfheParameters p{};
    initYatfhe(p);

    const int N = p.N, L = p.l, la = p.lApprox, reps = 2000;

    TorusPolynomial a{N};
    for (int i = 0; i < N; i++) {
        a.coeffs[i] = static_cast<Torus>(i * 2654435761u);
    }

    NttPolynomial A{N}, B{N}, C{N};
    applyNtt(A, a);
    applyNtt(B, a);

    auto t0 = steady_clock::now();
    for (int r = 0; r < reps; r++) {
        applyNtt(A, a);
    }
    auto t1 = steady_clock::now();
    double ntt = duration<double, std::micro>(t1 - t0).count() / reps;

    t0 = steady_clock::now();
    for (int r = 0; r < reps; r++) {
        calModularInnerProductNtt(C, A, B);
    }
    t1 = steady_clock::now();
    double pw = duration<double, std::micro>(t1 - t0).count() / reps;

    // The scalar gadget decomposition SS performs: N coefficients x L digits
    DecomposedData d{L};
    t0 = steady_clock::now();
    for (int r = 0; r < reps; r++) {
        for (int j = 0; j < N; j++) {
            signedGadgetDecomposition(d, a.coeffs[j], p);
        }
    }
    t1 = steady_clock::now();
    double gd = duration<double, std::micro>(t1 - t0).count() / reps;

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "N=" << N << "  l=" << L << "  lApprox=" << la << "   (k=1)\n\n";
    std::cout << "  NTT (N=" << N << ")                     "
              << std::setw(7) << ntt << " us\n";
    std::cout << "  pointwise mult-accumulate      "
              << std::setw(7) << pw << " us\n";
    std::cout << "  scalar gadget decomp, N x l=" << L << "  "
              << std::setw(7) << gd
              << " us   <- Sec III-B counts this as O(N) adds\n";

    std::cout << "\nREBUILT MODEL FROM MEASURED CONSTANTS (per key component):\n";
    double ep = la * 2 * ntt + la * 4 * pw + 2 * ntt;  // level*(k+1) NTT + level*(k+1)^2 pw + INTT
    double ss = la * (gd + L * ntt + 2 * L * pw + ntt);  // level SS calls, each: GD + L NTT + 2L pw + INTT

    std::cout << "  T_EP  = lApprox*2*NTT + lApprox*4*pw + 2*NTT            = "
              << std::setw(6) << ep << " us  (measured 17.92)\n";
    std::cout << "  T_SS  = lApprox*(GD + l*NTT + 2l*pw + INTT), 4-way par  = "
              << std::setw(6) << ss / la * 1.0 << " us  (measured 44.07)\n";
    std::cout << "     (per-thread, i.e. /lApprox since the " << la
              << " calls run in parallel)\n";
    std::cout << "\n  gadget decomposition alone is "
              << std::setprecision(0)
              << 100 * gd / (gd + L * ntt + 2 * L * pw + ntt)
              << "% of one SS call\n";

    return 0;
}