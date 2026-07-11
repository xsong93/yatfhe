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
    writePOD(os, t.N);
    serialize(t.b, os);
    writePOD(os, static_cast<int>(t.a.size()));
    for (const auto& poly : t.a) serialize(poly, os);
}

void deserialize(Trlwe& t, std::istream& is) {
    readPOD(is, t.k);
    readPOD(is, t.N);
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
    } else {
        // Ensure a stale value from a previous use of this object (e.g. a ping-ponged
        // buffer being deserialized into repeatedly) can't survive: without this, a
        // later resize(k, ...) onto an already-correctly-sized `a` is a no-op, leaving
        // old NTT data in place for the next accumulation step to add onto.
        t.a.clear();
    }
}

// Serialize TrlevDft
void serialize(const TrlevDft& t, std::ostream& os) {
    writePOD(os, t.l);
    writePOD(os, static_cast<int>(t.trlweDfts.size()));
    for (const auto& trlweDft : t.trlweDfts) serialize(trlweDft, os);
}

void deserialize(TrlevDft& t, std::istream& is) {
    readPOD(is, t.l);
    int trlweDfts_size;
    readPOD(is, trlweDfts_size);
    t.trlweDfts.resize(trlweDfts_size);
    for (auto& trlweDft : t.trlweDfts) deserialize(trlweDft, is);
}

// Serialize TrgswMP
void serialize(const TrgswMP& t, std::ostream& os) {
    writePOD(os, t.l);
    writePOD(os, t.k);
    writePOD(os, t.isHalf);
    // Serialize c (vector<vector<Trlwe>>)
    if (!t.isHalf) {
        writePOD(os, static_cast<int>(t.c.size()));
        for (const auto &inner: t.c) {
            writePOD(os, static_cast<int>(inner.size()));
            for (const auto &trlwe: inner) serialize(trlwe, os);
        }
    }
    // Serialize cPrime (vector<Trlwe>)
    writePOD(os, static_cast<int>(t.cPrime.size()));
    for (const auto& trlwe : t.cPrime) serialize(trlwe, os);
}

void deserialize(TrgswMP& t, std::istream& is) {
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

void serializeBskWWL24(const BootstrappingKeyWWL24& t, const std::string& filename) {
    std::ofstream os(filename, std::ios::binary | std::ios::trunc);
    for (auto i = 0; i < t.n; i++) {
        serialize(t.bskDft[i], os);
    }
    serialize(t.s2Dft, os);
    writePOD(os, t.n);
    writePOD(os, t.group);
    os.close();
}

void deserializeBskWWL24(BootstrappingKeyWWL24& bsk, const std::string& filename, const int n) {
    std::ifstream inFile(filename, std::ios::binary);
    if (!inFile) throw std::runtime_error("Failed to open file");

    bsk.bskDft.resize(n);  // Initialize outer vector (size = n)

    for (int i = 0; i < n; i++) {
        deserialize(bsk.bskDft[i], inFile);  // Read TrgswMPDft in order
    }
    deserialize(bsk.s2Dft, inFile);

    // Read remaining fields (n, group, isHalf) if they were serialized
    readPOD(inFile, bsk.n);
    readPOD(inFile, bsk.group);
    inFile.close();
}

void serializeBskMP(const BootstrappingKeyMP& t, const std::string& filename) {
    std::ofstream os(filename, std::ios::binary | std::ios::trunc);
    for (auto i = 0; i < t.n; i++) {
        serialize(t.bskDft[i][0], os);
    }
    writePOD(os, t.n);
    writePOD(os, t.group);
    writePOD(os, t.isHalf);
    os.close();
}

void deserializeBskMP(BootstrappingKeyMP& bsk, const std::string& filename, const int n) {
    std::ifstream inFile(filename, std::ios::binary);
    if (!inFile) throw std::runtime_error("Failed to open file");

    bsk.bskDft.resize(n);  // Initialize outer vector (size = n)

    for (int i = 0; i < n; i++) {
        bsk.bskDft[i].resize(1);       // Initialize inner vector (size = 1)
        deserialize(bsk.bskDft[i][0], inFile);  // Read TrgswMPDft in order
    }

    // Read remaining fields (n, group, isHalf) if they were serialized
    readPOD(inFile, bsk.n);
    readPOD(inFile, bsk.group);
    readPOD(inFile, bsk.isHalf);
    inFile.close();
}

void serializeBskMPOpt(const BootstrappingKeyMPOpt& t, const std::string& filename) {
    std::ofstream os(filename, std::ios::binary | std::ios::trunc);
    serialize(t.bskFirst[0], os);
    for (auto i = 0; i < t.n - 1; i++) {
        serialize(t.bskDft[i][0], os);
    }
    writePOD(os, t.n);
    writePOD(os, t.group);
    os.close();
}

void deserializeBskMPOpt(BootstrappingKeyMPOpt& bskOpt, const std::string& filename, const int n) {
    std::ifstream inFile(filename, std::ios::binary);
    if (!inFile) throw std::runtime_error("Failed to open file");

    bskOpt.bskFirst.resize(1);      // bskFirst is a vector of size 1
    bskOpt.bskDft.resize(n - 1);     // bskDft has n-1 elements

    // Step 1: Deserialize bskFirst[0]
    deserialize(bskOpt.bskFirst[0], inFile);

    // Step 2: Deserialize bskDft[i][0] for i = 0 to n-2
    for (int i = 0; i < n - 1; i++) {
        bskOpt.bskDft[i].resize(1);  // Each bskDft[i] is a vector of size 1
        deserialize(bskOpt.bskDft[i][0], inFile);
    }

    // Step 3: Deserialize remaining fields (n, group)
    readPOD(inFile, bskOpt.n);
    readPOD(inFile, bskOpt.group);
    inFile.close();
}

void serializeBskMPLazy(const BootstrappingKeyMPLazy& t, const std::string& filename) {
    std::ofstream os(filename, std::ios::binary | std::ios::trunc);
    serialize(t.bskFirst[0], os);
    for (auto i = 0; i < t.n - 1; i++) {
        serialize(t.bskTrim[i][0], os);
    }
    for (auto i = 0; i < t.n - 1; i++) {
        for (auto lvl = 0; lvl < t.level; lvl++) {
            serializeNestedVector(t.bskDecompA[t.decompIndex(i, lvl)], os);
        }
    }
    writePOD(os, t.n);
    writePOD(os, t.level);
    writePOD(os, t.group);
    writePOD(os, t.initialized);
    os.close();
}

void deserializeBskMPLazy(BootstrappingKeyMPLazy& bskLazy, const std::string& filename, const int n) {
    std::ifstream inFile(filename, std::ios::binary);
    if (!inFile) throw std::runtime_error("Failed to open file");

    bskLazy.bskFirst.resize(1);       // bskFirst is a vector of size 1
    bskLazy.bskTrim.resize(n - 1);    // bskTrim has n-1 elements

    // Step 1: Deserialize bskFirst[0]
    deserialize(bskLazy.bskFirst[0], inFile);

    // Step 2: Deserialize bskTrim[i][0] for i = 0 to n-2
    for (int i = 0; i < n - 1; ++i) {
        bskLazy.bskTrim[i].resize(1);  // Each bskTrim[i] is a vector of size 1
        deserialize(bskLazy.bskTrim[i][0], inFile);
    }

    // level/group aren't known until read from metadata at the end of the file, but each
    // TrgswMPDft already carries its own level count, and group is always 1 in this
    // (non-TERNARY) format, so we can recover what decompIndex() needs here.
    bskLazy.level = bskLazy.bskTrim[0][0].l;
    bskLazy.group = 1;
    bskLazy.bskDecompA.resize(static_cast<size_t>(n - 1) * bskLazy.level); // bskDecompA has (n-1)*level cells

    // Step 3: Deserialize bskDecompA[i][lvl] for i = 0 to n-2, lvl = 0 to level-1
    for (int i = 0; i < n - 1; ++i) {
        for (int lvl = 0; lvl < bskLazy.level; ++lvl) {
            deserializeNestedVector(bskLazy.bskDecompA[bskLazy.decompIndex(i, lvl)], inFile);
        }
    }

    // Step 4: Deserialize metadata (n, level, group, initialized)
    readPOD(inFile, bskLazy.n);
    readPOD(inFile, bskLazy.level);
    readPOD(inFile, bskLazy.group);
    readPOD(inFile, bskLazy.initialized);
    inFile.close();
}

void serializeBskLazyPipe(const BootstrappingKeyMPLazyPipe& t, const std::string& filename) {
    std::ofstream os(filename, std::ios::binary | std::ios::trunc);

    writePOD(os, t.level);

    // Step 1: Write bskFirst[0] and s2 (written outside loop)
    serialize(t.bskFirst[0], os);
    serialize(t.s2Dft, os);

    // Step 2: Write alternating pattern
    for (int i = 0; i < t.n - 1; ++i) {

        // Write bsk
        for (auto lvl = 0; lvl < t.level; lvl++) {
            serializeNestedVector(t.bskDecompA[t.decompIndex(i, lvl)], os);
            // serialize(t.bskB[t.decompIndex(i, lvl)], os);
            serializeNestedVector(t.bskDecompB[t.decompIndex(i, lvl)], os);
        }
    }
    os.close();
}

void deserializeBskLazyPipe(BootstrappingKeyMPLazyPipe& bskLazy, const std::string& filename, const int n) {
    std::ifstream inFile(filename, std::ios::binary);
    if (!inFile) throw std::runtime_error("Failed to open file");

    readPOD(inFile, bskLazy.level);

    // Resize vectors to match serialization structure
    bskLazy.bskFirst.resize(1);

    // Step 1: Read bskFirst[0] and s2 (written outside loop)
    deserialize(bskLazy.bskFirst[0], inFile);
    deserialize(bskLazy.s2Dft, inFile);

    bskLazy.group = 1;
    bskLazy.bskDecompA.resize(static_cast<size_t>(n - 1) * bskLazy.level);
    // bskLazy.bskB.resize(static_cast<size_t>(n - 1) * bskLazy.level);
    bskLazy.bskDecompB.resize(static_cast<size_t>(n - 1) * bskLazy.level);

    // Step 2: Read alternating pattern
    for (int i = 0; i < n - 1; ++i) {

        // Read bsk
        for (auto lvl = 0; lvl < bskLazy.level; lvl++) {
            deserializeNestedVector(bskLazy.bskDecompA[bskLazy.decompIndex(i, lvl)], inFile);
            // deserialize(bskLazy.bskB[bskLazy.decompIndex(i, lvl)], inFile);
            deserializeNestedVector(bskLazy.bskDecompB[bskLazy.decompIndex(i, lvl)], inFile);
        }
    }

    inFile.close();
}

void serializeBskLazyPipeAlt(const BootstrappingKeyMPLazyPipeAlt& t, const std::string& filename) {
    std::ofstream os(filename, std::ios::binary | std::ios::trunc);

    // Step 1: Write bskFirst[0] and s2Dft
    serialize(t.bskFirst[0], os);
    serialize(t.s2Dft, os);

    // Step 2: Write bskPrime[i][0]
    for (int i = 0; i < t.n - 1; ++i) {
        serialize(t.bskPrime[i], os);
    }

    // Step 3: Write meta data
    writePOD(os, t.level);
    os.close();
}

void deserializeBskLazyPipeAlt(BootstrappingKeyMPLazyPipeAlt& bskLazy, const std::string& filename, const int n) {
    std::ifstream inFile(filename, std::ios::binary);
    if (!inFile) throw std::runtime_error("Failed to open file");

    // Resize vectors to match serialization structure
    bskLazy.bskFirst.resize(1);
    bskLazy.bskPrime.resize(n - 1);

    // Step 1: Read bskFirst[0] and s2Dft
    deserialize(bskLazy.bskFirst[0], inFile);
    deserialize(bskLazy.s2Dft, inFile);

    // Step 2: Read bskPrime[i]
    for (int i = 0; i < n - 1; ++i) {
        deserialize(bskLazy.bskPrime[i], inFile);
    }

    // Step 3: Read meta data
    readPOD(inFile, bskLazy.level);
    inFile.close();
}