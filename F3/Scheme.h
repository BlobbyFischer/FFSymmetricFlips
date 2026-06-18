//
// Created by isaacw on 5/19/26.
//

#ifndef FFSYMMETRICFLIPS_SCHEME_H
#define FFSYMMETRICFLIPS_SCHEME_H
#include <cstdint>
#include <vector>

#include "FastRNG.h"
#include "RankOneTensor.h"

struct FlipCandidate {
    uint16_t index1;
    uint16_t index2;
    uint8_t axis; // 0 for match at a, 1 for match at b, 2 for match at c
    uint8_t g; // applying g to index2 gives a match at axis
    bool isNegEqual; // do we only have equality up to a sign?
};

template <typename Sym1, typename Sym2, typename Sym3>
class Scheme {
public:
    using TensorType = RankOneTensor<Sym1, Sym2, Sym3>;
    std::vector<TensorType> tensors;
    std::vector<FlipCandidate> flipCandidates;

    void generateFlipCandidates() {
        flipCandidates.clear();
        for (uint16_t t1_index=0; t1_index<tensors.size(); t1_index++) {
            const auto& t1 = tensors[t1_index];
            auto t1_orbit = t1.generateOrbit();
            uint64_t supp_t1_a = t1.a.support();
            uint64_t supp_t1_b = t1.b.support();
            uint64_t supp_t1_c = t1.c.support();

            for (uint16_t i = t1_index + 1; i < tensors.size(); i++) {
                const auto& t2 = tensors[i];
                bool checkA = (supp_t1_a == t2.a.support());
                bool checkB = (supp_t1_b == t2.b.support());
                bool checkC = (supp_t1_c == t2.c.support());
                if (!checkA && !checkB && !checkC) continue;
                if (checkA) {
                    for (uint8_t g=0; g < Sym1::ORBIT_SIZE; g++) {
                        const auto& mapped_t1 = t1_orbit[g];
                        if (mapped_t1.a == t2.a) {
                            flipCandidates.push_back({i,t1_index,0,g,false});
                        } else if (mapped_t1.a.negEqual(t2.a)) {
                            flipCandidates.push_back({i,t1_index,0,g,true});
                        }
                    }
                }if (checkB) {
                    for (uint8_t g=0; g < Sym1::ORBIT_SIZE; g++) {
                        const auto& mapped_t1 = t1_orbit[g];
                        if (mapped_t1.b == t2.b) {
                            flipCandidates.push_back({i,t1_index,1,g,false});
                        } else if (mapped_t1.b.negEqual(t2.b)) {
                            flipCandidates.push_back({i,t1_index,1,g,true});
                        }
                    }
                }if (checkC) {
                    for (uint8_t g=0; g < Sym1::ORBIT_SIZE; g++) {
                        const auto& mapped_t1 = t1_orbit[g];
                        if (mapped_t1.c == t2.c) {
                            flipCandidates.push_back({i,t1_index,2,g,false});
                        } else if (mapped_t1.c.negEqual(t2.c)) {
                            flipCandidates.push_back({i,t1_index,2,g,true});
                        }
                    }
                }
            }
        }
    }

    void generateFlipCandidatesForTensor(uint16_t targetIndex, uint16_t skipIndex = 0xFFFF) {
        const auto& t1 = tensors[targetIndex];
        auto t1_orbit = t1.generateOrbit();
        uint64_t supp_t1_a = t1.a.support();
        uint64_t supp_t1_b = t1.b.support();
        uint64_t supp_t1_c = t1.c.support();

        for (uint16_t i = 0; i < tensors.size(); i++) {
            if (i == targetIndex || i == skipIndex) continue;
            const auto& t2 = tensors[i];
            bool checkA = (supp_t1_a == t2.a.support());
            bool checkB = (supp_t1_b == t2.b.support());
            bool checkC = (supp_t1_c == t2.c.support());
            if (!checkA && !checkB && !checkC) continue;
            if (checkA) {
                for (uint8_t g=0; g < Sym1::ORBIT_SIZE; g++) {
                    const auto& mapped_t1 = t1_orbit[g];
                    if (mapped_t1.a == t2.a) {
                        flipCandidates.push_back({i,targetIndex,0,g,false});
                    } else if (mapped_t1.a.negEqual(t2.a)) {
                        flipCandidates.push_back({i,targetIndex,0,g,true});
                    }
                }
            }if (checkB) {
                for (uint8_t g=0; g < Sym1::ORBIT_SIZE; g++) {
                    const auto& mapped_t1 = t1_orbit[g];
                    if (mapped_t1.b == t2.b) {
                        flipCandidates.push_back({i,targetIndex,1,g,false});
                    } else if (mapped_t1.b.negEqual(t2.b)) {
                        flipCandidates.push_back({i,targetIndex,1,g,true});
                    }
                }
            }if (checkC) {
                for (uint8_t g=0; g < Sym1::ORBIT_SIZE; g++) {
                    const auto& mapped_t1 = t1_orbit[g];
                    if (mapped_t1.c == t2.c) {
                        flipCandidates.push_back({i,targetIndex,2,g,false});
                    } else if (mapped_t1.c.negEqual(t2.c)) {
                        flipCandidates.push_back({i,targetIndex,2,g,true});
                    }
                }
            }
        }
    }

    void flip(const FlipCandidate& f, FastRNG& rng) {
        const bool swap = rng.coinToss();
        uint16_t idx1 = swap ? f.index2 : f.index1;
        uint16_t idx2 = swap ? f.index1 : f.index2;
        auto& t1  = tensors[idx1];
        auto& t2 = tensors[idx2];

        F3Vec t2_a_mapped = Sym1::apply(t2.a, f.g);
        F3Vec t2_b_mapped = Sym2::apply(t2.b, f.g);
        F3Vec t2_c_mapped = Sym3::apply(t2.c, f.g);

        switch (f.axis) {
            case 0: {
                if (f.isNegEqual) t2_b_mapped = F3Vec(t2_b_mapped.minus, t2_b_mapped.plus);
                if (rng.coinToss()) {
                    t2.c = t2_c_mapped - t1.c;
                    t2.b = t2_b_mapped;
                    t2.a = t1.a;
                    t1.b += t2_b_mapped;
                } else {
                    t2.c = t2_c_mapped + t1.c;
                    t2.b = t2_b_mapped;
                    t2.a = t1.a;
                    t1.b -= t2_b_mapped;
                }
                break;
            }
            case 1: {
                if (f.isNegEqual) t2_c_mapped = F3Vec(t2_c_mapped.minus, t2_c_mapped.plus);
                if (rng.coinToss()) {
                    t2.a = t2_a_mapped - t1.a;
                    t2.c = t2_c_mapped;
                    t2.b = t1.b;
                    t1.c += t2_c_mapped;
                } else {
                    t2.a = t2_a_mapped + t1.a;
                    t2.c = t2_c_mapped;
                    t2.b = t1.b;
                    t1.c -= t2_c_mapped;
                }
                break;
            }
            case 2: {
                if (f.isNegEqual) t2_a_mapped = F3Vec(t2_a_mapped.minus, t2_a_mapped.plus);
                if (rng.coinToss()) {
                    t2.b = t2_b_mapped - t1.b;
                    t2.a = t2_a_mapped;
                    t2.c = t1.c;
                    t1.a += t2_a_mapped;
                } else {
                    t2.b = t2_b_mapped + t1.b;
                    t2.a = t2_a_mapped;
                    t2.c = t1.c;
                    t1.a -= t2_a_mapped;
                }
                break;
            }
        }
    }
    void plus(FastRNG& rng) {
        if (tensors.size() < 2) return;
        // choose the two indices
        uint16_t idx1 = rng.randomInt(tensors.size());
        uint16_t idx2 = rng.randomInt(tensors.size());
        while (idx1 == idx2) idx2 = rng.randomInt(tensors.size()); // make sure they are actually different
        uint8_t g = rng.randomInt(Sym1::ORBIT_SIZE); // choose a random group element
        F3Vec a1 = tensors[idx1].a;
        F3Vec b1 = tensors[idx1].b;
        F3Vec c1 = tensors[idx1].c;
        F3Vec a2 = Sym1::apply(tensors[idx2].a, g);
        F3Vec b2 = Sym2::apply(tensors[idx2].b, g);
        F3Vec c2 = Sym3::apply(tensors[idx2].c, g);
        switch (rng.randomInt(3)) {
            case 0:
                tensors[idx1] = TensorType(a1,b1+b2,c1);
                tensors[idx2] = TensorType(a1,b2,c2-c1);
                tensors.emplace_back(a2-a1,b2,c2);
                break;
            case 1:
                tensors[idx1] = TensorType(a1,b1,c1+c2);
                tensors[idx2] = TensorType(a2-a1,b1,c2);
                tensors.emplace_back(a2,b2-b1,c2);
                break;
            case 2:
                tensors[idx1] = TensorType(a1+a2,b1,c1);
                tensors[idx2] = TensorType(a2,b2-b1,c1);
                tensors.emplace_back(a2,b2,c2-c1);
                break;
        }
    }

    // Now we want to try to add "hidden flip" functionality

    bool isOrbitZero(const F3Vec& a, const F3Vec& b, const F3Vec& c) const { // perhaps there is a better way to do this
        F3Vec a1 = Sym1::apply(a, 1); F3Vec b1 = Sym2::apply(b, 1); F3Vec c1 = Sym3::apply(c, 1);
        F3Vec a2 = Sym1::apply(a, 2); F3Vec b2 = Sym2::apply(b, 2); F3Vec c2 = Sym3::apply(c, 2);
        F3Vec a3 = Sym1::apply(a, 3); F3Vec b3 = Sym2::apply(b, 3); F3Vec c3 = Sym3::apply(c, 3);
        uint64_t supp_a = a.support();
        while (supp_a) {
            int i = __builtin_ctzll(supp_a);
            supp_a &= supp_a - 1; // This just clears the lowest non-zero bit
            int val_a0 = a.get(i); int val_a1 = a1.get(i); int val_a2 = a2.get(i); int val_a3 = a3.get(i);
            uint64_t supp_b = b.support();
            while (supp_b) {
                int j = __builtin_ctzll(supp_b);
                supp_b &= supp_b - 1;
                int val_b0 = b.get(j); int val_b1 = b1.get(j); int val_b2 = b2.get(j); int val_b3 = b3.get(j);
                uint64_t supp_c = c.support();
                while (supp_c) {
                    int k = __builtin_ctzll(supp_c);
                    supp_c &= supp_c - 1;
                    int val_c0 = c.get(k);
                    int sum = val_a0 * val_b0 * val_c0 +
                              val_a1 * val_b1 * c1.get(k) +
                              val_a2 * val_b2 * c2.get(k) +
                              val_a3 * val_b3 * c3.get(k);

                    if (sum % 3 != 0) return false;
                }
            }
        }
        return true;
    }

    bool findAndApplyMutatedFlip(FastRNG& rng) {
        struct HiddenFlip {
            uint16_t idx1; uint16_t idx2;
            uint8_t axis; uint8_t g;
            bool isNegEqual;
            F3Vec delta;
        };
        std::vector<HiddenFlip> hiddenFlips;
        for (uint16_t i = 0; i < tensors.size(); ++i) {
            const auto& t1 = tensors[i];
            for (uint16_t j = 0; j < tensors.size(); ++j) {
                if (i == j) continue;
                const auto& t2 = tensors[j];
                for (uint8_t g = 0; g < Sym1::ORBIT_SIZE; ++g) {
                    F3Vec g_t2_a = Sym1::apply(t2.a, g);
                    F3Vec neg_g_t2_a = F3Vec(g_t2_a.minus, g_t2_a.plus); // Safe unary minus
                    F3Vec delta_a_pos = g_t2_a - t1.a;
                    F3Vec delta_a_neg = neg_g_t2_a - t1.a;
                    if (delta_a_pos.support() != 0 && isOrbitZero(delta_a_pos, t1.b, t1.c)) {
                        hiddenFlips.push_back({i, j, 0, g, false, delta_a_pos});
                    } else if (delta_a_neg.support() != 0 && isOrbitZero(delta_a_neg, t1.b, t1.c)) {
                        hiddenFlips.push_back({i, j, 0, g, true, delta_a_neg});
                    }
                    F3Vec g_t2_b = Sym2::apply(t2.b, g);
                    F3Vec neg_g_t2_b = F3Vec(g_t2_b.minus, g_t2_b.plus);
                    F3Vec delta_b_pos = g_t2_b - t1.b;
                    F3Vec delta_b_neg = neg_g_t2_b - t1.b;
                    if (delta_b_pos.support() != 0 && isOrbitZero(t1.a, delta_b_pos, t1.c)) {
                        hiddenFlips.push_back({i, j, 1, g, false, delta_b_pos});
                    } else if (delta_b_neg.support() != 0 && isOrbitZero(t1.a, delta_b_neg, t1.c)) {
                        hiddenFlips.push_back({i, j, 1, g, true, delta_b_neg});
                    }
                    F3Vec g_t2_c = Sym3::apply(t2.c, g);
                    F3Vec neg_g_t2_c = F3Vec(g_t2_c.minus, g_t2_c.plus);
                    F3Vec delta_c_pos = g_t2_c - t1.c;
                    F3Vec delta_c_neg = neg_g_t2_c - t1.c;
                    if (delta_c_pos.support() != 0 && isOrbitZero(t1.a, t1.b, delta_c_pos)) {
                        hiddenFlips.push_back({i, j, 2, g, false, delta_c_pos});
                    } else if (delta_c_neg.support() != 0 && isOrbitZero(t1.a, t1.b, delta_c_neg)) {
                        hiddenFlips.push_back({i, j, 2, g, true, delta_c_neg});
                    }
                }
            }
        }
        if (hiddenFlips.empty()) return false;
        auto chosen = hiddenFlips[rng.randomInt(hiddenFlips.size())];
        if (chosen.axis == 0) tensors[chosen.idx1].a += chosen.delta;
        else if (chosen.axis == 1) tensors[chosen.idx1].b += chosen.delta;
        else if (chosen.axis == 2) tensors[chosen.idx1].c += chosen.delta;
        FlipCandidate f = {chosen.idx1, chosen.idx2, chosen.axis, chosen.g, chosen.isNegEqual};
        this->flip(f, rng);
        return true;
    }
};


#endif //FFSYMMETRICFLIPS_SCHEME_H
