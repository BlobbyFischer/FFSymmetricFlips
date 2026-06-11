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
using SymA = SignFlipSymmetry<4, 0x202020202fdfd02, 0x60606060606f906>;
using SymB = SignFlipSymmetry<4, 0x282828282828d728, 0x606060606f9f906>;
using SymC = SignFlipSymmetry<4, 0x606f906f9060606, 0x202020202fdfd02>;

using MyScheme = Scheme<SymA, SymB, SymC>;

void randomWalk(MyScheme& scheme, int maxIterations, int n, int m, int p) {
    std::random_device rd;
    FastRNG rng(rd(), rd()); // only use the slow randomness here to generate random seeds for the fast prng
    std::string dirStr = "solutions/" + std::to_string(n) + "," + std::to_string(m) + "," + std::to_string(p); // this is where we will save the results
    std::filesystem::create_directories(dirStr); // in case it doesn't exist already
    scheme.generateFlipCandidates(); // ensure this is correct
    uint64_t best_rank = scheme.tensors.size(); // perhaps this is too big?
    for (int step = 0; step < maxIterations; ++step) {
        if (scheme.flipCandidates.empty() || (rng.randomInt(3000) == 0 && scheme.tensors.size() - best_rank < 4)) {
            for (uint64_t i = 0; i < rng.randomInt(3) + 1; ++i) {
                // loop so that we can increase by more at a time, hopefully avoiding deeper local minima
                scheme.plus(rng);
            }
            scheme.generateFlipCandidates();
        }
        uint64_t randIdx = rng.randomInt(scheme.flipCandidates.size());
        FlipCandidate chosenFlip = scheme.flipCandidates[randIdx];
        scheme.flip(chosenFlip, rng);
        uint16_t idx1 = chosenFlip.index1;
        uint16_t idx2 = chosenFlip.index2;
        if (scheme.tensors[idx1].isZero() || scheme.tensors[idx2].isZero()) {
            // this should happen rarely, so we can afford to be a little slow here
            for (int i = scheme.tensors.size() - 1; i >= 0; --i) {
                if (scheme.tensors[i].isZero()) {
                    scheme.tensors[i] = scheme.tensors.back();
                    scheme.tensors.pop_back();
                }
            }
            scheme.generateFlipCandidates();
            // did we find a new best?
            if (scheme.tensors.size() < best_rank) {
                best_rank = scheme.tensors.size();
                // removed verification. Can verify using a separate function
                std::string filename = Utils::genSchemeHash(scheme);
                std::filesystem::create_directories(dirStr + "/rank" + std::to_string(best_rank));
                Utils::saveRaw(scheme, dirStr + "/rank" + std::to_string(best_rank) + "/" + filename +".txt");
                Utils::saveReadable(scheme, dirStr + "/rank" + std::to_string(best_rank) + "/" + filename+".exp");
                std::cout << filename << ", " << scheme.tensors.size() << std::endl;
            }
        }
        else { // need to update flipCandidates
            // remove all outdated flips involving idx1 and idx2
            for (int i = scheme.flipCandidates.size() - 1; i >= 0; --i) {
                if (scheme.flipCandidates[i].index1 == idx1 ||
                    scheme.flipCandidates[i].index2 == idx1 ||
                    scheme.flipCandidates[i].index1 == idx2 ||
                    scheme.flipCandidates[i].index2 == idx2) {
                    scheme.flipCandidates[i] = scheme.flipCandidates.back();
                    scheme.flipCandidates.pop_back();
                }
            }
            // add all updated flips involving idx1 and idx2
            scheme.generateFlipCandidatesForTensor(idx1);
            scheme.generateFlipCandidatesForTensor(idx2,idx1);
        }
    }
}

int main(int argc, char* argv[]) {
    if (argc < 5) {
        std::cout << "Usage: " << argv[0] << " <dim1> <dim2> <dim3> <path_length> [filename]" << std::endl;
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
        std::string filename = argv[5];
        std::cout << "Loading state from: " << filename << "\n";
        Utils::loadRaw(activeScheme, filename);
    } else {
        std::cout << "Generating standard algorithm target as starting state.\n";
        // Deep copy the target into our active workspace
        activeScheme = targetScheme;
    }

    randomWalk(activeScheme, path_length, n, m, p);
    return 0;
}