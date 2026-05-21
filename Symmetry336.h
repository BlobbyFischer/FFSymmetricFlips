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
