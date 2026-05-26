//
// Created by isaacw on 5/26/26.
//

#include <iostream>
#include <string>
#include <vector>

#include "F3Vec.h"
#include "Symmetry336.h"
#include "RankOneTensor.h"
#include "Scheme.h"
#include "Utils.h"

// These are the exact symmetries Smirnov's <3,3,6> has
using SymA = SignFlipSymmetry<4, 0x0000000000050502ULL, 0x0000000000060106ULL>;
using SymB = SignFlipSymmetry<4, 0x0000000000281728ULL, 0x0000000000393906ULL>;
using SymC = SignFlipSymmetry<4, 0x0000010601060606ULL, 0x0000020202050502ULL>;

using MyScheme = Scheme<SymA, SymB, SymC>;

int main(int argc, char** argv) {
    if (argc < 5) {
        std::cout << "Usage: " << argv[0] << " <n> <m> <p> <file_to_verify.txt>" << std::endl;
        return 1;
    }

    int n = std::stoi(argv[1]);
    int m = std::stoi(argv[2]);
    int p = std::stoi(argv[3]);
    std::string filename = argv[4];

    MyScheme targetScheme;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            for (int k = 0; k < p; k++) {
                uint64_t a = 1ULL << (8 * i + j);
                uint64_t b = 1ULL << (8 * j + k);
                uint64_t c = 1ULL << (8 * k + i);
                targetScheme.tensors.emplace_back(F3Vec(a,0), F3Vec(b,0), F3Vec(c,0));
            }
        }
    }
    MyScheme testScheme;
    try {
        Utils::loadRaw(testScheme, filename);
    } catch (const std::exception& e) {
        std::cerr << "Failed to load file: " << e.what() << std::endl;
        return 1;
    }
    if (Utils::verifyDecomposition(testScheme, targetScheme)) {
        std::cout << "CORRECT" << std::endl;
    } else {
        std::cout << "INCORRECT" << std::endl;
    }
    return 0;
}