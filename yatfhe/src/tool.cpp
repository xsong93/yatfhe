//
// Created by Xintong Song on 2024/2/18.
//
#include <iostream>
#include <vector>
#include "yautil/tool.h"
#include "yatfhe/trlwe.h"

using namespace std;

void printTlweAB(const Tlwe& in, const string& msg) {
    cout << msg<< ": a: [";
    for (int i = 0; i < in.n; i++) {
        cout << i << ":" << in.a[i] <<" ";
    }
    cout <<"]" << endl << "b: [" << in.b << "]" << endl << endl;
}

void printTrlweAB(const Trlwe& in, const string& msg) {
    cout << msg<< ": a: [";
    for (int i = 0; i < in.k; i++) {
        for (int j = 0; j < in.b.N; j++) {
            cout << i << "," << j << ":" << in.a[i].coeffs[j] <<" ";
        }
    }
    cout <<"]" << endl << "b: [";
    for (int j = 0; j < in.b.N; j++) {
        cout << j << ":" << in.b.coeffs[j] <<" ";
    }
    cout <<"]" << endl << endl;
}

void printPolyMat(const vector<vector<IntPolynomial>>& in, const string& msg) {
    cout << msg <<": [";
    for (int i = 0; i < in.size(); i++) {
        for (int j = 0; j < in[i].size(); j++) {
            for (int k = 0; k < in[i][j].N; k++) {
                cout << i << "," << j << "," << k << ":" << in[i][j].coeffs[k] << " ";
            }
        }
    }
    cout <<"]" <<endl << endl;
}

void printBanner(const string& msg) {
    string l = ">>>>>>>>>>>>>>>>>>>>>>>> ";
    string r = " test passed! <<<<<<<<<<<<<<<<<<<<<<<<";
    auto length = l.size() + r.size() + msg.size();
    for (auto i = 0 ; i < length; i++) {
        std::cout << "-";
    }
    std::cout << std::endl;
    std::cout << l + msg + r << std::endl;
    for (auto i = 0 ; i < length; i++) {
        std::cout << "-";
    }
    std::cout << std::endl;
}