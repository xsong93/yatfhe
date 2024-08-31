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
        printf("%d:%s%d%s ", i, ANSI_COLOR_YELLOW, in[i], ANSI_COLOR_RESET);
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
