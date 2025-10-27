//
// Created by xintong on 10/27/25.
//

#include "yautil/ya_serializer.h"

// Serialize NttPolynomial
void serialize(const NttPolynomial& p, std::ostream& os) {
    writePOD(os, p.N);
    os.write(reinterpret_cast<const char*>(p.coeffs.data()), p.N * sizeof(NttType));
}

void deserialize(NttPolynomial& p, std::istream& is) {
    readPOD(is, p.N);
    p.coeffs.resize(p.N);
    is.read(reinterpret_cast<char*>(p.coeffs.data()), p.N * sizeof(NttType));
}

// Serialize TorusPolynomial
void serialize(const TorusPolynomial& p, std::ostream& os) {
    writePOD(os, p.N);
    os.write(reinterpret_cast<const char*>(p.coeffs.data()), p.N * sizeof(Torus));
}

void deserialize(TorusPolynomial& p, std::istream& is) {
    readPOD(is, p.N);
    p.coeffs.resize(p.N);
    is.read(reinterpret_cast<char*>(p.coeffs.data()), p.N * sizeof(Torus));
}

// Serialize DecompPolynomial
void serialize(const DecompPolynomial& p, std::ostream& os) {
    writePOD(os, p.N);
    os.write(reinterpret_cast<const char*>(p.coeffs.data()), p.N * sizeof(Decomp));
}

void deserialize(DecompPolynomial& p, std::istream& is) {
    readPOD(is, p.N);
    p.coeffs.resize(p.N);
    is.read(reinterpret_cast<char*>(p.coeffs.data()), p.N * sizeof(Decomp));
}

// Serialize Trlwe
void serialize(const Trlwe& t, std::ostream& os) {
    writePOD(os, t.k);
    serialize(t.b, os);
    writePOD(os, static_cast<int>(t.a.size()));
    for (const auto& poly : t.a) serialize(poly, os);
}

void deserialize(Trlwe& t, std::istream& is) {
    readPOD(is, t.k);
    deserialize(t.b, is);
    int a_size;
    readPOD(is, a_size);
    t.a.resize(a_size);
    for (auto& poly : t.a) deserialize(poly, is);
}

// Serialize TrlweDft
void serialize(const TrlweDft& t, std::ostream& os) {
    writePOD(os, t.k);
    writePOD(os, t.a_initialized);
    serialize(t.b, os);
    if (t.a_initialized) {
        writePOD(os, static_cast<int>(t.a.size()));
        for (const auto& poly : t.a) serialize(poly, os);
    }
}

void deserialize(TrlweDft& t, std::istream& is) {
    readPOD(is, t.k);
    readPOD(is, t.a_initialized);
    deserialize(t.b, is);
    if (t.a_initialized) {
        int a_size;
        readPOD(is, a_size);
        t.a.resize(a_size);
        for (auto& poly : t.a) deserialize(poly, is);
    }
}

// Serialize TrgswMPDft
void serialize(const TrgswMPDft& t, std::ostream& os) {
    writePOD(os, t.l);
    writePOD(os, t.k);
    writePOD(os, t.isHalf);
    // Serialize c (vector<vector<TrlweDft>>)
    if (!t.isHalf) {
        writePOD(os, static_cast<int>(t.c.size()));
        for (const auto &inner: t.c) {
            writePOD(os, static_cast<int>(inner.size()));
            for (const auto &trlwe: inner) serialize(trlwe, os);
        }
    }
    // Serialize cPrime (vector<TrlweDft>)
    writePOD(os, static_cast<int>(t.cPrime.size()));
    for (const auto& trlwe : t.cPrime) serialize(trlwe, os);
}

void deserialize(TrgswMPDft& t, std::istream& is) {
    readPOD(is, t.l);
    readPOD(is, t.k);
    readPOD(is, t.isHalf);
    // Deserialize c
    int outer_size;
    readPOD(is, outer_size);
    t.c.resize(outer_size);
    for (auto& inner : t.c) {
        int inner_size;
        readPOD(is, inner_size);
        inner.resize(inner_size);
        for (auto& trlwe : inner) deserialize(trlwe, is);
    }
    // Deserialize cPrime
    int cPrime_size;
    readPOD(is, cPrime_size);
    t.cPrime.resize(cPrime_size);
    for (auto& trlwe : t.cPrime) deserialize(trlwe, is);
}