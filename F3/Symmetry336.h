//
// Created by isaacw on 5/19/26.
//

#pragma once
#include <cstdint>

#include "F3Vec.h"

#ifndef FFSYMMETRICFLIPS_SYMMETRY336_H
#define FFSYMMETRICFLIPS_SYMMETRY336_H

// The template parameters define the distinct masks for this specific symmetry component
template <int OrbitSize, uint64_t Gen1Mask, uint64_t Gen2Mask>
struct SignFlipSymmetry {
    static constexpr int ORBIT_SIZE = OrbitSize;

    static constexpr uint64_t MASK_11 = Gen1Mask & Gen2Mask;
    static constexpr uint64_t MASK_10 = Gen1Mask & ~Gen2Mask;
    static constexpr uint64_t MASK_01 = ~Gen1Mask & Gen2Mask;
    static constexpr uint64_t MASK_00 = ~Gen1Mask & ~Gen2Mask;

    static inline uint64_t get_footprint(uint64_t support) { // returns which of the four signature classes the support is in
        uint8_t output = 0;
        if (support & MASK_00) output |= 1;
        if (support & MASK_01) output |= 2;
        if (support & MASK_10) output |= 4;
        if (support & MASK_11) output |= 8;
        return output;
    }

    static inline uint64_t get_kernel_mask(uint8_t fp_other1, uint8_t fp_other2) { // we know the footprints of the other two terms, then what is the kernel of this term? i.e fix b and c, then what are the a such that a*b*c is 0 (when symmetrized)        uint8_t forbidden = 0;
        uint8_t forbidden = 0;
        for (int i=0; i<4; i++) {
            if ((fp_other1 >> i) & 1) {
                for (int j = 0; j < 4; j++) {
                    if ((fp_other2 >> j) & 1) {
                        forbidden |= (1 << (i ^ j)); // this just tells us which symmetries we *can't* use if we want to be in the kernel
                    }
                }
            }
        }
        uint64_t kernel = 0;
        if (!(forbidden & 1)) kernel |= MASK_00;
        if (!(forbidden & 2)) kernel |= MASK_01;
        if (!(forbidden & 4)) kernel |= MASK_10;
        if (!(forbidden & 8)) kernel |= MASK_11;
        return kernel;
    }

    [[nodiscard]] static F3Vec apply_mask(const F3Vec& v, uint64_t mask) {
        return F3Vec((v.plus & ~mask) | (v.minus & mask), (v.plus & mask) | (v.minus & ~mask));
    }

    [[nodiscard]] static F3Vec apply(const F3Vec& v, int g) {
        switch(g) {
            case 0: return v;
            case 1: return apply_mask(v, Gen1Mask);
            case 2: return apply_mask(v, Gen2Mask);
            case 3: return apply_mask(v, Gen1Mask ^ Gen2Mask);
            default: return v;
        }
    }
};

#endif //FFSYMMETRICFLIPS_SYMMETRY336_H
