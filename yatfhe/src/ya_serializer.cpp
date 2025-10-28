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

    // Only deserialize c if !isHalf (mirror serialization)
    if (!t.isHalf) {
        int outer_size;
        readPOD(is, outer_size);
        t.c.resize(outer_size);
        for (auto& inner : t.c) {
            int inner_size;
            readPOD(is, inner_size);
            inner.resize(inner_size);
            for (auto& trlwe : inner) deserialize(trlwe, is);
        }
    } else {
        t.c.clear();  // Ensure c is empty if isHalf
    }

    // Deserialize cPrime
    int cPrime_size;
    readPOD(is, cPrime_size);
    t.cPrime.resize(cPrime_size);
    for (auto& trlwe : t.cPrime) deserialize(trlwe, is);
}

void deserializeBskMP(BootstrappingKeyMP& bsk, const std::string& filename, int n) {
    std::ifstream inFile(filename, std::ios::binary);
    bsk.bskDft.resize(n);  // Initialize outer vector (size = n)

    for (int i = 0; i < n; ++i) {
        bsk.bskDft[i].resize(1);       // Initialize inner vector (size = 1)
        deserialize(bsk.bskDft[i][0], inFile);  // Read TrgswMPDft in order
    }

    // Read remaining fields (n, group, isHalf) if they were serialized
    inFile.read(reinterpret_cast<char*>(&bsk.n), sizeof(int));
    inFile.read(reinterpret_cast<char*>(&bsk.group), sizeof(int));
    inFile.read(reinterpret_cast<char*>(&bsk.isHalf), sizeof(bool));
}

void deserializeBskMPOpt(BootstrappingKeyMPOpt& bskOpt, const std::string& filename, int n) {
    std::ifstream inFile(filename, std::ios::binary);
    if (!inFile) throw std::runtime_error("Failed to open file");

    bskOpt.bskFirst.resize(1);      // bskFirst is a vector of size 1
    bskOpt.bskDft.resize(n - 1);     // bskDft has n-1 elements

    // Step 1: Deserialize bskFirst[0]
    deserialize(bskOpt.bskFirst[0], inFile);

    // Step 2: Deserialize bskDft[i][0] for i = 0 to n-2
    for (int i = 0; i < n - 1; ++i) {
        bskOpt.bskDft[i].resize(1);  // Each bskDft[i] is a vector of size 1
        deserialize(bskOpt.bskDft[i][0], inFile);
    }

    // Step 3: Deserialize remaining fields (n, group, isHalf, initialized)
    inFile.read(reinterpret_cast<char*>(&bskOpt.n), sizeof(int));
    inFile.read(reinterpret_cast<char*>(&bskOpt.group), sizeof(int));
    inFile.read(reinterpret_cast<char*>(&bskOpt.isHalf), sizeof(bool));
    inFile.read(reinterpret_cast<char*>(&bskOpt.initialized), sizeof(bool));
}

void deserializeBskLazy(BootstrappingKeyMPLazy& bskLazy, const std::string& filename, int n) {
    std::ifstream inFile(filename, std::ios::binary);
    if (!inFile) throw std::runtime_error("Failed to open file");

    bskLazy.bskFirst.resize(1);       // bskFirst is a vector of size 1
    bskLazy.bskTrim.resize(n - 1);    // bskTrim has n-1 elements
    bskLazy.bskDecompA.resize(n - 1); // bskDecompA has n-1 elements

    // Step 1: Deserialize bskFirst[0]
    deserialize(bskLazy.bskFirst[0], inFile);

    // Step 2: Deserialize bskTrim[i][0] for i = 0 to n-2
    for (int i = 0; i < n - 1; ++i) {
        bskLazy.bskTrim[i].resize(1);  // Each bskTrim[i] is a vector of size 1
        deserialize(bskLazy.bskTrim[i][0], inFile);
    }

    // Step 3: Deserialize bskDecompA[i][0] for i = 0 to n-2
    for (int i = 0; i < n - 1; ++i) {
        bskLazy.bskDecompA[i].resize(1);  // Each bskDecompA[i] is a vector of size 1
        deserializeNestedVector(bskLazy.bskDecompA[i][0], inFile);
    }

    // Step 4: Deserialize metadata (n, level, group, initialized)
    inFile.read(reinterpret_cast<char*>(&bskLazy.n), sizeof(int));
    inFile.read(reinterpret_cast<char*>(&bskLazy.level), sizeof(int));
    inFile.read(reinterpret_cast<char*>(&bskLazy.group), sizeof(int));
    inFile.read(reinterpret_cast<char*>(&bskLazy.initialized), sizeof(bool));
}

void deserializeBskLazyPipe(BootstrappingKeyMPLazy& bskLazy, const std::string& filename, int n) {
    std::ifstream inFile(filename, std::ios::binary);
    if (!inFile) throw std::runtime_error("Failed to open file");

    bskLazy.bskFirst.resize(1);       // bskFirst is a vector of size 1
    bskLazy.bskTrim.resize(n - 1);    // bskTrim has n-1 elements
    bskLazy.bskDecompA.resize(n - 2); // bskDecompA has n-1 elements

    // Step 1: Deserialize bskFirst[0], bskFull[0]
    deserialize(bskLazy.bskFirst[0], inFile);
    deserialize(bskLazy.bskTrim[0][0], inFile);

    // Step 2: Deserialize bskTrim[i][0] for i = 1 to n-3,
    // and bskDecompA[i][0] for i = 1 to n-4
    for (int i = 1; i < n - 2; ++i) {
        bskLazy.bskTrim[i].resize(1);  // Each bskTrim[i] is a vector of size 1
        deserialize(bskLazy.bskTrim[i][0], inFile);
        if (i < n - 3) {
            bskLazy.bskDecompA[i].resize(1);   // Each bskDecompA[i] is a vector of size 1
            deserializeNestedVector(bskLazy.bskDecompA[i][0], inFile);
        }
    }

    // Step 4: Deserialize metadata (n, level, group, initialized)
    inFile.read(reinterpret_cast<char*>(&bskLazy.n), sizeof(int));
    inFile.read(reinterpret_cast<char*>(&bskLazy.level), sizeof(int));
    inFile.read(reinterpret_cast<char*>(&bskLazy.group), sizeof(int));
    inFile.read(reinterpret_cast<char*>(&bskLazy.initialized), sizeof(bool));
}