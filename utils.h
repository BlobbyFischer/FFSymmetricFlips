#pragma once
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include "Scheme.h"

namespace Utils {

    // --- VERIFICATION LOGIC ---

    // Projects the entire Scheme (including all orbits) into a dense 64x64x64 grid
    template<typename SchemeType>
    std::vector<int8_t> buildDenseTensor(const SchemeType& scheme) {
        std::vector<int8_t> dense(64 * 64 * 64, 0);
        
        for (const auto& t : scheme.tensors) {
            for (const auto& ot : t.generateOrbit()) {
                for (int i = 0; i < 64; ++i) {
                    int va = ot.a.get(i);
                    if (!va) continue;
                    va = (va == 2) ? -1 : 1; // Map GF(3) '2' to standard '-1'
                    
                    for (int j = 0; j < 64; ++j) {
                        int vb = ot.b.get(j);
                        if (!vb) continue;
                        vb = (vb == 2) ? -1 : 1;
                        
                        for (int k = 0; k < 64; ++k) {
                            int vc = ot.c.get(k);
                            if (!vc) continue;
                            vc = (vc == 2) ? -1 : 1;
                            
                            int idx = (i * 64 + j) * 64 + k;
                            dense[idx] += va * vb * vc;
                        }
                    }
                }
            }
        }
        
        // Finalize GF(3) arithmetic
        for (auto& val : dense) {
            val %= 3;
            if (val < 0) val += 3;
            if (val == 2) val = -1; // Keep as -1, 0, 1 for easy debugging
        }
        return dense;
    }

    template<typename SchemeType>
    bool verifyDecomposition(const SchemeType& current, const SchemeType& target) {
        return buildDenseTensor(current) == buildDenseTensor(target);
    }

    // --- RAW DATA I/O ---

    template<typename SchemeType>
    void saveRaw(const SchemeType& s, const std::string& filename) {
        std::ofstream out(filename);
        // Write out the raw hex masks for perfect serialization
        for (const auto& t : s.tensors) {
            out << std::hex << t.a.plus << " " << t.a.minus << " "
                << t.b.plus << " " << t.b.minus << " "
                << t.c.plus << " " << t.c.minus << "\n";
        }
    }

    template<typename SchemeType>
    void loadRaw(SchemeType& s, const std::string& filename) {
        std::ifstream in(filename);
        if (!in) throw std::runtime_error("Could not open file: " + filename);
        
        s.tensors.clear();
        uint64_t ap, am, bp, bm, cp, cm;
        // Parse the hex masks directly back into the bit-sliced architecture
        while (in >> std::hex >> ap >> am >> bp >> bm >> cp >> cm) {
            s.tensors.emplace_back(F3Vec(ap, am), F3Vec(bp, bm), F3Vec(cp, cm));
        }
    }

    // --- HUMAN READABLE EXPORT ---

    // Helper to print a single vector in formal basis notation (e.g., +a11 -a12)
    void printReadableVec(const F3Vec& v, std::ofstream& out, char basisName) {
        bool first = true;
        out << "(";
        for (int i = 0; i < 64; ++i) {
            int val = v.get(i);
            if (val) {
                if (!first) out << " ";

                // Calculate 1-indexed matrix coordinates (1 to 8)
                int row = (i / 8) + 1;
                int col = (i % 8) + 1;

                out << (val == 1 ? "+" : "-") << basisName << row << col;
                first = false;
            }
        }

        if (first) out << "0"; // Fallback if the vector is completely empty
        out << ")";
    }

    template<typename SchemeType>
    void saveReadable(const SchemeType& s, const std::string& filename) {
        std::ofstream out(filename);
        out << "Scheme Rank: " << s.tensors.size() << "\n";
        out << "====================================\n";

        for (size_t i = 0; i < s.tensors.size(); ++i) {
            const auto& t = s.tensors[i];
            out << "Orbit " << i + 1 << ":\n";
            auto orbit = t.generateOrbit();

            // Print each element of the orbit as a full tensor product
            for (size_t g = 0; g < orbit.size(); ++g) {
                out << "  ";
                printReadableVec(orbit[g].a, out, 'a');
                out << " * ";
                printReadableVec(orbit[g].b, out, 'b');
                out << " * ";
                printReadableVec(orbit[g].c, out, 'c');
                out << "\n";
            }
            out << "------------------------------------\n";
        }
    }
}