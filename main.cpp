#include <iostream>
#include <string>

#include "F3Vec.h"
#include "Symmetry336.h"
#include "RankOneTensor.h"
#include "Scheme.h"
#include "FastRNG.h"
#include "utils.h"

using SymA = SignFlipSymmetry<4, 0xFF00FF00FF00FF00ULL, 0xAAAAAAAAAAAAAAAAULL>;
using SymB = SignFlipSymmetry<4, 0x0000000000000000ULL, 0xFF00FF00FF00FF00ULL>;
using SymC = SignFlipSymmetry<4, 0xAAAAAAAAAAAAAAAAULL, 0x0000000000000000ULL>;

using MyScheme = Scheme<SymA, SymB, SymC>;

void runRandomWalk(MyScheme& scheme, const MyScheme& targetScheme, int maxIterations) {
    FastRNG rng(1337, 4242);

    std::cout << "Starting Random Walk. Initial rank: " << scheme.tensors.size() << "\n";

    for (int step = 0; step < maxIterations; ++step) {
        scheme.generateFlipCandidates();

        if (scheme.flipCandidates.empty()) {
            std::cout << "Walk terminated: Local Minimum.\n";
            break;
        }

        uint64_t randomIndex = rng.randomInt(scheme.flipCandidates.size());
        scheme.flip(scheme.flipCandidates[randomIndex], rng);

        bool rankReduced = false;
        for (int i = scheme.tensors.size() - 1; i >= 0; --i) {
            if (scheme.tensors[i].isZero()) {
                scheme.tensors[i] = scheme.tensors.back();
                scheme.tensors.pop_back();
                rankReduced = true;
            }
        }

        if (rankReduced) {
            std::cout << "Step " << step << " | Rank reduced! New rank: " << scheme.tensors.size() << "\n";

            // Output the state whenever we hit a reduction!
            Utils::saveRaw(scheme, "current_state.txt");
            Utils::saveReadable(scheme, "current_readable.txt");

            /*if (scheme.tensors.size() <= 7) {
                std::cout << "Target Rank found! Stopping.\n";
                break;
            }*/
        }
    }

    std::cout << "\n--- Verification ---\n";
    if (Utils::verifyDecomposition(scheme, targetScheme)) {
        std::cout << "[SUCCESS] The output scheme mathematically equals the target tensor!\n";
    } else {
        std::cout << "[FAILED] The output scheme is broken. Math does not match.\n";
    }
}

int main(int argc, char* argv[]) {
    MyScheme targetScheme;

    // Build the mathematical target <2,2,2> for verification purposes
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            for (int k = 0; k < 2; ++k) {
                uint64_t a_mask = 1ULL << (8 * i + j);
                uint64_t b_mask = 1ULL << (8 * j + k);
                uint64_t c_mask = 1ULL << (8 * k + i);
                targetScheme.tensors.emplace_back(F3Vec(a_mask, 0), F3Vec(b_mask, 0), F3Vec(c_mask, 0));
            }
        }
    }

    MyScheme activeScheme;

    // Check if a file was passed as an argument
    if (argc > 1) {
        std::string filename = argv[1];
        std::cout << "Loading state from: " << filename << "\n";
        Utils::loadRaw(activeScheme, filename);
    } else {
        std::cout << "No file provided. Generating <2,2,2> target as starting state.\n";
        // Deep copy the target into our active workspace
        activeScheme = targetScheme;
    }

    runRandomWalk(activeScheme, targetScheme, 100000);

    return 0;
}