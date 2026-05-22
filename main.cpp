#include <iostream>
#include <string>
#include <random>
#include <filesystem>

#include "F3Vec.h"
#include "Symmetry336.h"
#include "RankOneTensor.h"
#include "Scheme.h"
#include "FastRNG.h"
#include "Utils.h"

// These are the exact symmetries Smirnov's <3,3,6> has
using SymA = SignFlipSymmetry<4, 0x0000000000050502ULL, 0x0000000000060106ULL>;
using SymB = SignFlipSymmetry<4, 0x0000000000281728ULL, 0x0000000000393906ULL>;
using SymC = SignFlipSymmetry<4, 0x0000010601060606ULL, 0x0000020202050502ULL>;

using MyScheme = Scheme<SymA, SymB, SymC>;

void randomWalk(MyScheme& scheme, const MyScheme& targetScheme, int maxIterations, int n, int m, int p) {
    std::random_device rd;
    FastRNG rng(rd(), rd());
    std::string dirStr = "solutions/" + std::to_string(n) + "," + std::to_string(m) + "," + std::to_string(p); // this is where we will save the results
    std::filesystem::create_directories(dirStr); // in case it doesn't exist already
    for (int step = 0; step < maxIterations; ++step) {
        scheme.generateFlipCandidates();
        if (scheme.flipCandidates.empty()) {
            break; // there are no flips left to be done
        }
        uint64_t randIdx = rng.randomInt(scheme.flipCandidates.size());
        scheme.flip(scheme.flipCandidates[randIdx], rng);
        for (int i = scheme.tensors.size() - 1; i >= 0; --i) {
            if (scheme.tensors[i].isZero()) {
                scheme.tensors[i] = scheme.tensors.back();
                scheme.tensors.pop_back();
            }
        }
    }
    if (Utils::verifyDecomposition(scheme,targetScheme)) {
        std::string filename = Utils::genSchemeHash(scheme);
        Utils::saveRaw(scheme, dirStr + "/" + filename+".txt");
        Utils::saveReadable(scheme, dirStr + "/" + filename+".exp");
        std::cout << filename << ", " << scheme.tensors.size() << std::endl;
    } else {
        std::cout << "ERROR: input and output don't match";
    }
}

int main(int argc, char* argv[]) {
    if (argc < 5) {
        std::cout << "Usage: " << argv[0] << "<dim1> <dim2> <dim3> <path_length> [filename]" << std::endl;
    }
    int n = std::stoi(argv[1]);
    int m = std::stoi(argv[2]);
    int p = std::stoi(argv[3]);

    int path_length = std::stoi(argv[4]);

    MyScheme targetScheme;

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            for (int k = 0; k < p; ++k) {
                uint64_t a_mask = 1ULL << (8 * i + j);
                uint64_t b_mask = 1ULL << (8 * j + k);
                uint64_t c_mask = 1ULL << (8 * k + i);
                targetScheme.tensors.emplace_back(F3Vec(a_mask, 0), F3Vec(b_mask, 0), F3Vec(c_mask, 0));
            }
        }
    }

    MyScheme activeScheme;

    // Check if a file was passed as an argument
    if (argc > 5) {
        std::string filename = argv[4];
        std::cout << "Loading state from: " << filename << "\n";
        Utils::loadRaw(activeScheme, filename);
    } else {
        std::cout << "Generating standard algorithm target as starting state.\n";
        // Deep copy the target into our active workspace
        activeScheme = targetScheme;
    }

    randomWalk(activeScheme, targetScheme, path_length, n, m, p);
    return 0;
}