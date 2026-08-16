//
// Created by Xintong Song on 2024/2/18.
//
#include <iostream>
#include <vector>
#include "yautil/tool.h"

#include <thread>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>

#include "yatfhe/trlwe.h"

using namespace std;

void printRlweAB(const Rlwe& in, const string& msg) {
#ifdef PRINTER_ON
    cout << ANSI_COLOR_CYAN << msg << ANSI_COLOR_RESET << ": a: ";
    for (int i = 0; i < in.k; i++) {
        cout << "[";
        for (int j = 0; j < in.b.N; j++) {
            printf("%d,%d:%s%ld%s  " , i, j, ANSI_COLOR_YELLOW, (int64_t)in.a[i].coeffs[j], ANSI_COLOR_RESET);
        }
        cout <<"] ";
    }
    cout << endl << "b: [";
    for (int j = 0; j < in.b.N; j++) {
        printf("%d:%s%ld%s  ", j, ANSI_COLOR_YELLOW, (int64_t)in.b.coeffs[j], ANSI_COLOR_RESET);
    }
    cout <<"]" << endl << endl;
#endif
}

void printTrlweDftAB(const TrlweDft& in, const string& msg) {
#ifdef PRINTER_ON
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
#endif
}

void printDecomposedTrlweAB(const DecomposedTrlwe& in, const string& msg) {
#ifdef PRINTER_ON
    cout << ANSI_COLOR_CYAN << msg << ": " << ANSI_COLOR_RESET;
    for (auto l = 0; l < in.l; l++) {
        printf("level %d: a: ", l);
        for (int i = 0; i < in.trlwes[l].k; i++) {
            cout << "[";
            for (int j = 0; j < in.trlwes[l].b.N; j++) {
                printf("%d,%d:%s%ld%s  ", i, j, ANSI_COLOR_YELLOW, (int64_t)in.trlwes[l].a[i].coeffs[j], ANSI_COLOR_RESET);
            }
            cout << "] ";
        }
        cout << endl << "b: [";
        for (int j = 0; j < in.trlwes[l].b.N; j++) {
            printf("%d:%s%ld%s  ", j, ANSI_COLOR_YELLOW, (int64_t)in.trlwes[l].b.coeffs[j], ANSI_COLOR_RESET);

        }
        cout << "]" << endl << endl;
    }
#endif
}

void printDecomposedTrlweNttAB(const DecomposedTrlweDft& in, const string& msg) {
#ifdef PRINTER_ON
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
#endif
}

void printPolyMat(const vector<vector<IntPolynomial>>& in, const string& msg) {
#ifdef PRINTER_ON
    cout << ANSI_COLOR_CYAN << msg << ANSI_COLOR_RESET <<": [";
    for (int i = 0; i < in.size(); i++) {
        for (int j = 0; j < in[i].size(); j++) {
            for (int k = 0; k < in[i][j].N; k++) {
                cout << i << "," << j << "," << k << ":" << in[i][j].coeffs[k] << " ";
            }
        }
    }
    cout <<"]" <<endl << endl;
#endif
}

void printBanner(const string& msg) {
#ifdef PRINTER_ON
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
#endif
}

size_t clearFileCache(const std::string& filename) {
    const int fd = ::open(filename.c_str(), O_RDONLY);
    if (fd < 0) {
        return 0;   // no such file
    }
    struct stat st{};
    if (fstat(fd, &st) != 0 || st.st_size <= 0) {
        close(fd);
        return 0;
    }
    const auto len = static_cast<size_t>(st.st_size);
    // POSIX_FADV_DONTNEED drops clean pages only, so flush first
    fdatasync(fd);
    posix_fadvise(fd, 0, 0, POSIX_FADV_DONTNEED);

    size_t resident = 0;
    const auto pageSize = static_cast<size_t>(sysconf(_SC_PAGESIZE));
    void* map = mmap(nullptr, len, PROT_READ, MAP_SHARED, fd, 0);
    if (map != MAP_FAILED) {
        vector<unsigned char> present((len + pageSize - 1) / pageSize, 0);
        if (mincore(map, len, present.data()) == 0) {
            for (const unsigned char page : present) {
                resident += (page & 1u);
            }
        }
        munmap(map, len);
    }
    close(fd);
    if (resident != 0) {
        cerr << "clearFileCache(" << filename << "): " << resident
             << " pages still resident, so the next read is not cold" << endl;
    }
    return resident;
}

void clearFileCache() {
    // sudo -n so a missing password fails immediately instead of blocking on a prompt.
    if (system("sync; echo 3 | sudo -n tee /proc/sys/vm/drop_caches > /dev/null 2>&1") != 0) {
        static bool warned = false;
        if (!warned) {
            cerr << "clearFileCache(): dropping the page cache needs root, so nothing was dropped." << endl;
            warned = true;
        }
        return;
    }
    std::this_thread::sleep_for(std::chrono::milliseconds(10)); // Let it settle
}