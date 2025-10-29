//
// Created by xintong on 10/27/25.
//

#ifndef BASE_YA_SERIALIZER_H
#define BASE_YA_SERIALIZER_H

#include <fstream>
#include <vector>
#include <cstdint>
#include "yatfhe/polynomial.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/trgsw.h"
#include "yatfhe/bootstrapping.h"

void serialize(const NttPolynomial& p, std::ostream& os);

void deserialize(NttPolynomial& p, std::istream& is);

void serialize(const TorusPolynomial& p, std::ostream& os);

void serialize(const DecompPolynomial& p, std::ostream& os);

void deserialize(DecompPolynomial& p, std::istream& is);

void serialize(const Trlwe& t, std::ostream& os);

void deserialize(Trlwe& t, std::istream& is);

void deserialize(TorusPolynomial& p, std::istream& is);

void serialize(const TrlweDft& t, std::ostream& os);

void deserialize(TrlweDft& t, std::istream& is);

void serialize(const TrgswMPDft& t, std::ostream& os);

void deserialize(TrgswMPDft& t, std::istream& is);

void serializeBskMP(const BootstrappingKeyMP& t, const std::string& filename);

void serializeBskMPOpt(const BootstrappingKeyMPOpt& t, const std::string& filename);

void deserializeBskMP(BootstrappingKeyMP& bsk, const std::string& filename, int n);

void deserializeBskMPOpt(BootstrappingKeyMPOpt& bskOpt, const std::string& filename, int n);

void serializeBskMPLazy(const BootstrappingKeyMPLazy& t, const std::string& filename);

void deserializeBskMPLazy(BootstrappingKeyMPLazy& bskLazy, const std::string& filename, int n);

void deserializeBskLazyPipe(BootstrappingKeyMPLazyPipe& bskLazy, const std::string& filename, int n);

// Helper to write/read POD types and vectors
template <typename T>
void writePOD(std::ostream& os, const T& val) {
    os.write(reinterpret_cast<const char*>(&val), sizeof(T));
}

template <typename T>
void readPOD(std::istream& is, T& val) {
    is.read(reinterpret_cast<char*>(&val), sizeof(T));
}

template <typename T>
void serializeNestedVector(const std::vector<T>& vec, std::ostream& os) {
    // Write size of current dimension
    size_t size = vec.size();
    os.write(reinterpret_cast<const char*>(&size), sizeof(size_t));

    // Recursively serialize elements
    for (const auto& elem : vec) {
        serializeNestedVector(elem, os);  // Recursion unwinds at Int16Polynomial
    }
}

// Template specialization for DecompPolynomial (base case)
template <>
inline void serializeNestedVector<DecompPolynomial>(const std::vector<DecompPolynomial>& vec, std::ostream& os) {
    size_t size = vec.size();
    os.write(reinterpret_cast<const char*>(&size), sizeof(size_t));
    for (const auto& poly : vec) {
        serialize(poly, os);  // Calls Step 1
    }
}

template <typename T>
void deserializeNestedVector(std::vector<T>& vec, std::istream& is) {
    // Read size of current dimension
    size_t size;
    is.read(reinterpret_cast<char*>(&size), sizeof(size_t));
    vec.resize(size);

    // Recursively deserialize elements
    for (auto& elem : vec) {
        deserializeNestedVector(elem, is);
    }
}

// Template specialization for DecompPolynomial (base case)
template <>
inline void deserializeNestedVector<DecompPolynomial>(std::vector<DecompPolynomial>& vec, std::istream& is) {
    size_t size;
    is.read(reinterpret_cast<char*>(&size), sizeof(size_t));
    vec.resize(size);
    for (auto& poly : vec) {
        deserialize(poly, is);
    }
}

#endif //BASE_YA_SERIALIZER_H
