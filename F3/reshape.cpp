#include <iostream>
#include <string>
#include <filesystem>
#include <algorithm>

#include "F3Vec.h"
#include "Symmetry336.h"
#include "Scheme.h"
#include "Utils.h"

// Smirnov's symmetries again
using SymA = SignFlipSymmetry<4, 0x0000000000050502ULL, 0x0000000000060106ULL>;
using SymB = SignFlipSymmetry<4, 0x0000000000281728ULL, 0x0000000000393906ULL>;
using SymC = SignFlipSymmetry<4, 0x0000010601060606ULL, 0x0000020202050502ULL>;

using MyScheme = Scheme<SymA, SymB, SymC>;

int main(int argc, char* argv[]) {
    if (argc < 8) {
        std::cout << "Usage: " << argv[0] << " <n> <m> <p> <input_file.txt> <new_n> <new_m> <new_p>" << std::endl;
        return 1;
    }

    int n = std::stoi(argv[1]);
    int m = std::stoi(argv[2]);
    int p = std::stoi(argv[3]);
    std::string filename = argv[4];
    
    int new_n = std::stoi(argv[5]);
    int new_m = std::stoi(argv[6]);
    int new_p = std::stoi(argv[7]);

    MyScheme scheme;
    Utils::loadRaw(scheme, filename);

    for (int i = 0; i < new_n; ++i) {
        for (int j = 0; j < new_m; ++j) {
            for (int k = 0; k < new_p; ++k) {
                if (i < n && j < m && k < p) continue;
                uint64_t a = 1ULL << (8 * i + j);
                uint64_t b = 1ULL << (8 * j + k);
                uint64_t c = 1ULL << (8 * k + i);
                scheme.tensors.emplace_back(F3Vec(a, 0), F3Vec(b, 0), F3Vec(c, 0));
            }
        }
    }

    uint64_t a_mask = 0;
    uint64_t b_mask = 0;
    uint64_t c_mask = 0;

    for (int i=0; i < new_n; ++i) {
        for (int j=0; j < new_m; ++j) {
            a_mask |= (1ULL << (8 * i + j));
        }
    }

    for (int j=0; j < new_m; ++j) {
        for (int k=0; k < new_p; ++k) {
            b_mask |= (1ULL << (8 * j + k));
        }
    }

    for (int k=0; k < new_p; ++k) {
        for (int i=0; i < new_n; ++i) {
            c_mask |= (1ULL << (8 * k + i));
        }
    }

    for (auto& t : scheme.tensors) {
        t.a.plus &= a_mask;
        t.a.minus &= a_mask;
        t.b.plus &= b_mask;
        t.b.minus &= b_mask;
        t.c.plus &= c_mask;
        t.c.minus &= c_mask;
    }

    scheme.tensors.erase(
        std::remove_if(scheme.tensors.begin(), scheme.tensors.end(), [](const auto &t) { return t.isZero(); }),
        scheme.tensors.end()
    );

    std::string dirStr = "solutions/" + std::to_string(new_n) + "," + std::to_string(new_m) + "," + std::to_string(new_p);
    std::string rankDirStr = dirStr + "/rank" + std::to_string(scheme.tensors.size());
    std::filesystem::create_directories(rankDirStr);
    std::string hashStr = Utils::genSchemeHash(scheme);
    std::string rawOut = rankDirStr + "/" + hashStr + ".txt";
    std::string expOut = rankDirStr + "/" + hashStr + ".exp";

    Utils::saveRaw(scheme, rawOut);
    Utils::saveReadable(scheme, expOut);

    std::cout<< rawOut << "," << scheme.tensors.size() << "\n";

    return 0;
}