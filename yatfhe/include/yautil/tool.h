//
// Created by Xintong Song on 2024/2/18.
//

#ifndef HLS_YATFHE_TOOL_H
#define HLS_YATFHE_TOOL_H

#include <iostream>
#include "yatfhe/tlwe.h"
#include "yatfhe/trlwe.h"

using namespace std;

void printTlweAB(const Tlwe& in, const string& msg);

void printTrlweAB(const Trlwe& in, const string& msg);

template <typename T>
void printArray(const vector<T>& in, const string& msg) {
    cout << msg <<": [";
    for (int i = 0; i < in.size(); i++) {
        cout << i << ":" << in[i] <<" ";
    }
    cout <<"]" <<endl << endl;
}

void printPolyMat(const vector<vector<IntPolynomial>>& in, const string& msg);

template <typename T>
void printPolyVec(const vector<T>& in, const string& msg) {
    cout << msg << ": [";
    for (int i = 0; i < in.size(); i++) {
        for (int j = 0; j < in[0].N; j++) {
            cout << i << "," << j << ":" << in[i].coeffs[j] <<" ";
        }
    }
    cout <<"]" <<endl << endl;
}

#endif //HLS_YATFHE_TOOL_H
