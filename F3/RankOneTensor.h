//
// Created by isaacw on 5/19/26.
//

#ifndef FFSYMMETRICFLIPS_RANKONETENSOR_H
#define FFSYMMETRICFLIPS_RANKONETENSOR_H

#include <array>

#include "F3Vec.h"

template <typename Sym1, typename Sym2, typename Sym3> // for now this will be F3Vec and then the symmetries for the a, b, c components
class RankOneTensor {
public:
    F3Vec a, b, c;

    RankOneTensor() = default;

    RankOneTensor(const F3Vec& a, const F3Vec& b, const F3Vec& c) : a(a), b(b), c(c) {}

    [[nodiscard]] std::array<RankOneTensor, Sym1::ORBIT_SIZE> generateOrbit() const {
        std::array<RankOneTensor, Sym1::ORBIT_SIZE> orbit;
        for (int g=0; g<Sym1::ORBIT_SIZE; g++) {
            orbit[g].a = Sym1::apply(a, g);
            orbit[g].b = Sym2::apply(b, g);
            orbit[g].c = Sym3::apply(c, g);
        }
        return orbit;
    }

    [[nodiscard]] bool isZero() const {
        if (a.support() == 0) return true;
        if (b.support() == 0) return true;
        if (c.support() == 0) return true;
        return false;
    }

    [[nodiscard]] bool isGhost() const { // checks if the orbit is zero, not just one of the rank one tensors is zero
        uint8_t fp_a = Sym1::get_footprint(a.support());
        uint8_t fp_b = Sym2::get_footprint(b.support());
        uint8_t fp_c = Sym3::get_footprint(c.support());
        if ((a.support() & ~Sym1::get_kernel_mask(fp_b,fp_c)) == 0) return true;
        if ((b.support() & ~Sym2::get_kernel_mask(fp_a,fp_c)) == 0) return true;
        if ((c.support() & ~Sym3::get_kernel_mask(fp_a,fp_b)) == 0) return true;
        return false;
    }
};


#endif //FFSYMMETRICFLIPS_RANKONETENSOR_H
