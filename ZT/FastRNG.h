//
// Created by isaacw on 5/20/26.
//

#ifndef FFSYMMETRICFLIPS_FASTRNG_H
#define FFSYMMETRICFLIPS_FASTRNG_H


#pragma once
#include <cstdint>

class FastRNG {
private:
    uint64_t s[2];

    static inline uint64_t rotl(const uint64_t x, int k) {
        return (x << k) | (x >> (64 - k));
    }

public:
    FastRNG(uint64_t seed1, uint64_t seed2) {
        s[0] = seed1;
        s[1] = seed2;
    }

    uint64_t next() {
        const uint64_t s0 = s[0];
        uint64_t s1 = s[1];
        const uint64_t result = s0 + s1;

        s1 ^= s0;
        s[0] = rotl(s0, 24) ^ s1 ^ (s1 << 16);
        s[1] = rotl(s1, 37);

        return result;
    }

    bool coinToss() {
        return (next() & 1ULL) == 1ULL;
    }

    uint64_t randomInt(uint64_t max_val) {
        return next() % max_val;
    }
};

#endif //FFSYMMETRICFLIPS_FASTRNG_H
