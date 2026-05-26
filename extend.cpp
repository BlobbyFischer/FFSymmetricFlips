#include <iostream>
#include <string>
#include <filesystem>

#include "F3Vec.h"
#include "Symmetry336.h"
#include "RankOneTensor.h"
#include "Scheme.h"
#include "Utils.h"

// Use the exact same symmetries you are currently working with
using SymA = SignFlipSymmetry<4, 0x0000000000050502ULL, 0x0000000000060106ULL>;
using SymB = SignFlipSymmetry<4, 0x0000000000281728ULL, 0x0000000000393906ULL>;
using SymC = SignFlipSymmetry<4, 0x0000010601060606ULL, 0x0000020202050502ULL>;

using MyScheme = Scheme<SymA, SymB, SymC>;

int main(int argc, char* argv[]) {
    if (argc < 8) {
        std::cout << "Usage: " << argv[0] << " <input_file.txt> <n> <m> <p> <new_n> <new_m> <new_p>\n";
        return 1;
    }

    std::string filename = argv[1];
    int n = std::stoi(argv[2]);
    int m = std::stoi(argv[3]);
    int p = std::stoi(argv[4]);
    
    int new_n = std::stoi(argv[5]);
    int new_m = std::stoi(argv[6]);
    int new_p = std::stoi(argv[7]);

    MyScheme scheme;
    std::cout << "Loading state from: " << filename << "\n";
    Utils::loadRaw(scheme, filename);

    std::cout << "Original rank: " << scheme.tensors.size() << "\n";

    // Pad with the standard algorithm orbits
    for (int i = 0; i < new_n; ++i) {
        for (int j = 0; j < new_m; ++j) {
            for (int k = 0; k < new_p; ++k) {
                // If it's inside the original n x m x p box, skip it (it's already in the scheme)
                if (i < n && j < m && k < p) continue;

                // Otherwise, add the standard basis tensor for this coordinate
                uint64_t a_mask = 1ULL << (8 * i + j);
                uint64_t b_mask = 1ULL << (8 * j + k);
                uint64_t c_mask = 1ULL << (8 * k + i);
                scheme.tensors.emplace_back(F3Vec(a_mask, 0), F3Vec(b_mask, 0), F3Vec(c_mask, 0));
            }
        }
    }

    std::cout << "New padded rank: " << scheme.tensors.size() << "\n";

    // Create the appropriate target directory
    std::string dirStr = "solutions/" + std::to_string(new_n) + "," + std::to_string(new_m) + "," + std::to_string(new_p);
    std::string rankDirStr = dirStr + "/rank_" + std::to_string(scheme.tensors.size());
    std::filesystem::create_directories(rankDirStr);

    // Hash and save
    std::string hashStr = Utils::genSchemeHash(scheme);
    std::string rawOut = rankDirStr + "/" + hashStr + ".txt";
    std::string expOut = rankDirStr + "/" + hashStr + ".exp";

    Utils::saveRaw(scheme, rawOut);
    Utils::saveReadable(scheme, expOut);

    std::cout << "Saved padded scheme to:\n  " << rawOut << "\n";

    return 0;
}