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

    // CONSTRUCTORS

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

    bool flip(const FlipCandidate& f, FastRNG& rng) {
        const bool swap = rng.coinToss();
        uint16_t idx1 = swap ? f.index2 : f.index1;
        uint16_t idx2 = swap ? f.index1 : f.index2;
        auto& t1  = tensors[idx1];
        auto& t2 = tensors[idx2];

        ZTVec t2_a_mapped = Sym1::apply(t2.a, f.g);
        ZTVec t2_b_mapped = Sym2::apply(t2.b, f.g);
        ZTVec t2_c_mapped = Sym3::apply(t2.c, f.g);

        switch (f.axis) {
            case 0: {
                if (f.isNegEqual) t2_b_mapped = ZTVec(t2_b_mapped.minus, t2_b_mapped.plus);
                if (rng.coinToss()) {
                    if (!t2_c_mapped.canSubtract(t1.c) || !t1.b.canAdd(t2_b_mapped)) return false;
                    t2.c = t2_c_mapped - t1.c;
                    t2.b = t2_b_mapped;
                    t2.a = t1.a;
                    t1.b += t2_b_mapped;
                } else {
                    if (!t2_c_mapped.canAdd(t1.c) || !t1.b.canSubtract(t2_b_mapped)) return false;
                    t2.c = t2_c_mapped + t1.c;
                    t2.b = t2_b_mapped;
                    t2.a = t1.a;
                    t1.b -= t2_b_mapped;
                }
                break;
            }
            case 1: {
                if (f.isNegEqual) t2_c_mapped = ZTVec(t2_c_mapped.minus, t2_c_mapped.plus);
                if (rng.coinToss()) {
                    if (!t2_a_mapped.canSubtract(t1.a) || !t1.c.canAdd(t2_c_mapped)) return false;
                    t2.a = t2_a_mapped - t1.a;
                    t2.c = t2_c_mapped;
                    t2.b = t1.b;
                    t1.c += t2_c_mapped;
                } else {
                    if (!t2_a_mapped.canAdd(t1.a) || !t1.c.canSubtract(t2_c_mapped)) return false;
                    t2.a = t2_a_mapped + t1.a;
                    t2.c = t2_c_mapped;
                    t2.b = t1.b;
                    t1.c -= t2_c_mapped;
                }
                break;
            }
            case 2: {
                if (f.isNegEqual) t2_a_mapped = ZTVec(t2_a_mapped.minus, t2_a_mapped.plus);
                if (rng.coinToss()) {
                    if (!t2_b_mapped.canSubtract(t1.b) || !t1.a.canAdd(t2_a_mapped)) return false;
                    t2.b = t2_b_mapped - t1.b;
                    t2.a = t2_a_mapped;
                    t2.c = t1.c;
                    t1.a += t2_a_mapped;
                } else {
                    if (!t2_b_mapped.canAdd(t1.b) || !t1.a.canSubtract(t2_a_mapped)) return false;
                    t2.b = t2_b_mapped + t1.b;
                    t2.a = t2_a_mapped;
                    t2.c = t1.c;
                    t1.a -= t2_a_mapped;
                }
                break;
            }
        }
        return true;
    }
    void plus(FastRNG& rng) {
        if (tensors.size() < 2) return;
        // choose the two indices
        uint16_t idx1 = rng.randomInt(tensors.size());
        uint16_t idx2 = rng.randomInt(tensors.size());
        while (idx1 == idx2) idx2 = rng.randomInt(tensors.size()); // make sure they are actually different
        uint8_t g = rng.randomInt(Sym1::ORBIT_SIZE); // choose a random group element
        ZTVec a1 = tensors[idx1].a;
        ZTVec b1 = tensors[idx1].b;
        ZTVec c1 = tensors[idx1].c;
        ZTVec a2 = Sym1::apply(tensors[idx2].a, g);
        ZTVec b2 = Sym2::apply(tensors[idx2].b, g);
        ZTVec c2 = Sym3::apply(tensors[idx2].c, g);
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
};


#endif //FFSYMMETRICFLIPS_SCHEME_H
