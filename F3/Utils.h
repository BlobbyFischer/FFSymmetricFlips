#pragma once
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include "Scheme.h"

namespace Utils {

    // Projects the entire Scheme into a 64x64x64 tensor. Useful for checking correctness
    template<typename SchemeType>
    std::vector<int8_t> buildDenseTensor(const SchemeType& scheme) {
        std::vector<int8_t> dense(64 * 64 * 64, 0);
        for (const auto& t : scheme.tensors) {
            for (const auto& ot : t.generateOrbit()) {
                for (int i = 0; i < 64; ++i) {
                    int va = ot.a.get(i);
                    if (!va) continue;
                    va = (va == 2) ? -1 : 1;
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
        for (auto& val : dense) {
            val %= 3;
            if (val < 0) val += 3;
            if (val == 2) val = -1; // Keep in {-1, 0, 1}
        }
        return dense;
    }

    template<typename SchemeType>
    bool verifyDecomposition(const SchemeType& current, const SchemeType& target) {
        return buildDenseTensor(current) == buildDenseTensor(target);
    }

    template<typename SchemeType>
    void saveRaw(const SchemeType& s, const std::string& filename) {
        std::ofstream out(filename);
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

    // Print a vector in human readable format
    inline void readableVec(const F3Vec &v, std::ofstream &out, char basisName) {
        bool first = true;
        out << "(";
        for (int i = 0; i < 64; ++i) {
            if (int val = v.get(i)) {

                // Calculate 1-indexed matrix coordinates (1 to 8)
                const int row = (i / 8) + 1;
                const int col = (i % 8) + 1;

                out << (val == 1 ? (first ? " " : "+") : "-") << basisName << row << col;
                first = false;
            }
        }

        if (first) out << "0"; // This vector is empty. We should never need this...
        out << ")";
    }

    template<typename SchemeType>
    void saveReadable(const SchemeType& s, const std::string& filename) {
        std::ofstream out(filename);

        for (size_t i = 0; i < s.tensors.size(); ++i) {
            const auto& t = s.tensors[i];
            auto orbit = t.generateOrbit();
            for (size_t g = 0; g < orbit.size(); ++g) {
                out << "";
                readableVec(orbit[g].a, out, 'a');
                out << "*";
                readableVec(orbit[g].b, out, 'b');
                out << "*";
                readableVec(orbit[g].c, out, 'c');
                out << "\n";
            }
            out << "\n"; // Space out the orbits a bit. Technically not needed...
        }
    }

    // a hash for the scheme which is orbit invariant.
    template<typename SchemeType>
    std::string genSchemeHash(const SchemeType& scheme) {
        uint64_t hash = 0;
        for (const auto& t : scheme.tensors) {
            uint64_t sa = t.a.support();
            uint64_t sb = t.b.support();
            uint64_t sc = t.c.support();
            uint64_t tensorHash = sa;
            tensorHash = (tensorHash ^ (sb >> 32) ^ sb) * 0x9E3779B185EBCA87ULL;
            tensorHash = (tensorHash ^ (sc >> 32) ^ sc) * 0xC2B2AE3D27D4EB4FULL;
            hash += tensorHash;
        }
        std::stringstream ss;
        ss << std::hex << std::setfill('0') << std::setw(16) << hash;
        return ss.str();
    }
}