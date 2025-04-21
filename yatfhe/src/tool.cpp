//
// Created by Xintong Song on 2024/2/18.
//
#include <iostream>
#include <vector>
#include "yautil/tool.h"
#include "yatfhe/trlwe.h"

using namespace std;

void printRlweAB(const Rlwe& in, const string& msg) {
    cout << ANSI_COLOR_CYAN << msg << ANSI_COLOR_RESET << ": a: ";
    for (int i = 0; i < in.k; i++) {
        cout << "[";
        for (int j = 0; j < in.b.N; j++) {
            printf("%d,%d:%s%d%s  " , i, j, ANSI_COLOR_YELLOW, in.a[i].coeffs[j], ANSI_COLOR_RESET);
        }
        cout <<"] ";
    }
    cout << endl << "b: [";
    for (int j = 0; j < in.b.N; j++) {
        printf("%d:%s%d%s  ", j, ANSI_COLOR_YELLOW, in.b.coeffs[j], ANSI_COLOR_RESET);
    }
    cout <<"]" << endl << endl;
}

void printTrlweDftAB(const TrlweDft& in, const string& msg) {
    cout << ANSI_COLOR_CYAN << msg << ANSI_COLOR_RESET << ": a: ";
    for (int i = 0; i < in.k; i++) {
        cout << "[";
        for (int j = 0; j < in.b.N; j++) {
            printf("%d,%d:%s%lu%s  " , i, j, ANSI_COLOR_YELLOW, in.a[i].coeffs[j], ANSI_COLOR_RESET);
        }
        cout <<"] ";
    }
    cout << endl << "b: [";
    for (int j = 0; j < in.b.N; j++) {
        printf("%d:%s%lu%s  ", j, ANSI_COLOR_YELLOW, in.b.coeffs[j], ANSI_COLOR_RESET);
    }
    cout <<"]" << endl << endl;
}

void printDecomposedTrlweAB(const DecomposedTrlwe& in, const string& msg) {
    cout << ANSI_COLOR_CYAN << msg << ": " << ANSI_COLOR_RESET;
    for (auto l = 0; l < in.l; l++) {
        printf("level %d: a: ", l);
        for (int i = 0; i < in.trlwes[l].k; i++) {
            cout << "[";
            for (int j = 0; j < in.trlwes[l].b.N; j++) {
                printf("%d,%d:%s%d%s  ", i, j, ANSI_COLOR_YELLOW, in.trlwes[l].a[i].coeffs[j], ANSI_COLOR_RESET);
            }
            cout << "] ";
        }
        cout << endl << "b: [";
        for (int j = 0; j < in.trlwes[l].b.N; j++) {
            printf("%d:%s%d%s  ", j, ANSI_COLOR_YELLOW, in.trlwes[l].b.coeffs[j], ANSI_COLOR_RESET);

        }
        cout << "]" << endl << endl;
    }
}

void printDecomposedTrlweNttAB(const DecomposedTrlweDft& in, const string& msg) {
    cout << ANSI_COLOR_CYAN << msg << ": " << ANSI_COLOR_RESET;
    for (auto l = 0; l < in.l; l++) {
        printf("level %d: a: ", l);
        for (int i = 0; i < in.rlweDfts[l].k; i++) {
            cout << "[";
            for (int j = 0; j < in.rlweDfts[l].b.N; j++) {
                printf("%d,%d:%s%lu%s  ", i, j, ANSI_COLOR_YELLOW, in.rlweDfts[l].a[i].coeffs[j], ANSI_COLOR_RESET);
            }
            cout << "] ";
        }
        cout << endl << "b: [";
        for (int j = 0; j < in.rlweDfts[l].b.N; j++) {
            printf("%d:%s%lu%s  ", j, ANSI_COLOR_YELLOW, in.rlweDfts[l].b.coeffs[j], ANSI_COLOR_RESET);

        }
        cout << "]" << endl << endl;
    }
}

void printPolyMat(const vector<vector<IntPolynomial>>& in, const string& msg) {
    cout << ANSI_COLOR_CYAN << msg << ANSI_COLOR_RESET <<": [";
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
        std::cout << ANSI_COLOR_GREEN << "-";
    }
    std::cout << std::endl << ANSI_COLOR_GREEN << l + msg + r << std::endl;
    for (auto i = 0 ; i < length; i++) {
        std::cout << ANSI_COLOR_GREEN << "-";
    }
    std::cout << std::endl;
}