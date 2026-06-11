//
// Created by isaacw on 5/19/26.
//

#ifndef FFSYMMETRICFLIPS_RANKONETENSOR_H
#define FFSYMMETRICFLIPS_RANKONETENSOR_H

#include <array>

#include "ZTVec.h"

template <typename Sym1, typename Sym2, typename Sym3> // for now this will be ZTVec and then the symmetries for the a, b, c components
class RankOneTensor {
public:
    ZTVec a, b, c;

    RankOneTensor() = default;

    RankOneTensor(const ZTVec& a, const ZTVec& b, const ZTVec& c) : a(a), b(b), c(c) {}

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
        return a.support() == 0 || b.support() == 0 || c.support() == 0;
    }
};


#endif //FFSYMMETRICFLIPS_RANKONETENSOR_H
