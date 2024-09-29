//
// Created by Xintong Song on 2024/2/18.
//

#ifndef HLS_YATFHE_TOOL_H
#define HLS_YATFHE_TOOL_H

#include <iostream>
#include "yatfhe/tlwe.h"
#include "yatfhe/trlwe.h"
#include "yatfhe/trgsw.h"
#include "yautil/ansi_color.h"

using namespace std;

template <typename T>
void printElement(int index, const T& value) {
    if constexpr (std::is_integral<T>::value && std::is_signed<T>::value && !std::is_same<T, char>::value && !std::is_same<T, bool>::value) {
        // Handle signed integral types (int, int64_t, etc.)
        printf("%d:%s%lld%s ", index, ANSI_COLOR_YELLOW, static_cast<long long>(value), ANSI_COLOR_RESET);
    }
    else if constexpr (std::is_integral<T>::value && std::is_unsigned<T>::value && !std::is_same<T, char>::value) {
        // Handle unsigned integral types (uint32_t, uint64_t, etc.)
        printf("%d:%s%llu%s ", index, ANSI_COLOR_YELLOW, static_cast<unsigned long long>(value), ANSI_COLOR_RESET);
    }
    else if constexpr (std::is_floating_point<T>::value) {
        // Handle floating point types (float, double, etc.)
        printf("%d:%s%f%s ", index, ANSI_COLOR_YELLOW, value, ANSI_COLOR_RESET);
    }
    else if constexpr (std::is_same<T, char>::value) {
        // Handle char type
        printf("%d:%s'%c'%s ", index, ANSI_COLOR_YELLOW, value, ANSI_COLOR_RESET);
    }
    else if constexpr (std::is_same<T, bool>::value) {
        // Handle bool type
        printf("%d:%s%s%s ", index, ANSI_COLOR_YELLOW, value ? "true" : "false", ANSI_COLOR_RESET);
    }
    else if constexpr (std::is_same<T, std::string>::value) {
        // Handle std::string
        printf("%d:%s%s%s ", index, ANSI_COLOR_YELLOW, value.c_str(), ANSI_COLOR_RESET);
    }
    else {
        // Fallback for any other types
        std::cout << index << ":" << ANSI_COLOR_YELLOW << value << ANSI_COLOR_RESET << " ";
    }
}

template <typename T>
void printTlweAB(const T& in, const string& msg) {
    cout << msg << ": a: [";
    for (int i = 0; i < in.n; i++) {
        cout << i << ":" << ANSI_COLOR_YELLOW << in.a[i] << ANSI_COLOR_RESET << " ";
    }
    cout <<"]" << endl << "b: [" << ANSI_COLOR_YELLOW << in.b << ANSI_COLOR_RESET << "]" << endl << endl;
}

void printRlweAB(const Rlwe& in, const string& msg);

void printTrlweAB(const Trlwe& in, const string& msg);

void printTrlweDftAB(const TrlweDft& in, const string& msg);

void printTrgsw(const Trgsw& in , const string& msg);

void printDecomposedTrlweAB(const DecomposedTrlwe& in, const string& msg);

void printDecomposedTrlweNttAB(const DecomposedTrlweDft& in, const string& msg);

template <typename T>
void printArray(const vector<T>& in, const string& msg) {
    cout << msg <<": [";
    for (int i = 0; i < in.size(); i++) {
//        cout << i << ":" << ANSI_COLOR_YELLOW << in[i] << ANSI_COLOR_RESET <<" ";
//        printf("%d:%s%d%s ", i, ANSI_COLOR_YELLOW, in[i], ANSI_COLOR_RESET);
        printElement(i, in[i]);
    }
    cout <<"]" <<endl << endl;
}

template <typename T, size_t N>
void printArray(const T (&in)[N], const string& msg) {
    cout << msg <<": [";
    for (size_t i = 0; i < N; i++) {
//        cout << i << ":" << ANSI_COLOR_YELLOW << in[i] << ANSI_COLOR_RESET <<" ";
//        printf("%d:%s%d%s ", i, ANSI_COLOR_YELLOW, in[i], ANSI_COLOR_RESET);
        printElement(i, in[i]);
    }
    cout <<"]" <<endl << endl;
}

void printPolyMat(const vector<vector<IntPolynomial>>& in, const string& msg);

template <typename T>
void printPolyVec(const vector<T>& in, const string& msg) {
    cout << msg << ": [";
    for (int i = 0; i < in.size(); i++) {
        for (int j = 0; j < in[0].N; j++) {
            cout << i << "," << j << ":" << ANSI_COLOR_YELLOW <<  in[i].coeffs[j] << ANSI_COLOR_RESET  <<" ";
        }
    }
    cout <<"]" <<endl << endl;
}

void printBanner(const string& msg);

#endif //HLS_YATFHE_TOOL_H
