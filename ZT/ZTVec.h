//
// Created by isaacw on 5/19/26.
//

#ifndef FFSYMMETRICFLIPS_ZTVec_H
#define FFSYMMETRICFLIPS_ZTVec_H
#include <cstdint>
#include <stdexcept>


class ZTVec {
public:
    uint64_t plus; // is the number +1
    uint64_t minus; // is the number -1

    // get rid of numbers that are +1-1
    void cleanup() {
        uint64_t overlap = plus & minus;
        plus ^= overlap;
        minus ^= overlap;
    }

    ZTVec(): plus(0), minus(0) {}
    ZTVec(uint64_t const plus, uint64_t const minus): plus(plus), minus(minus) {
        cleanup();
    }

    [[nodiscard]] bool canAdd(const ZTVec &other) const {
        return (plus & other.plus) == 0 && (minus & other.minus) == 0;
    }

    ZTVec operator+(const ZTVec &other) const { // assumes we have already checked we can add
        uint64_t p = (plus ^ other.plus) | (minus & other.minus);
        uint64_t q = (minus ^ other.minus) | (plus & other.plus);
        return {p & ~q, q & ~p};
    }

    ZTVec operator+=(const ZTVec &other) {
        *this = *this + other;
        return *this;
    }

    [[nodiscard]] bool canSubtract(const ZTVec &other) const {
        return (plus & other.minus) == 0 && (minus & other.plus) == 0;
    }

    ZTVec operator-(const ZTVec &other) const { // assumes you have checked if you can subtract
        ZTVec negOther(other.minus, other.plus);
        return *this + negOther;
    }

    ZTVec operator-=(const ZTVec &other) {
        *this = *this - other;
        return *this;
    }

    bool operator==(const ZTVec &other) const {
        return plus == other.plus && minus == other.minus;
    }

    [[nodiscard]] bool negEqual(const ZTVec &other) const {
        return plus == other.minus && minus == other.plus;
    }

    bool operator!=(const ZTVec &other) const {
        return plus != other.plus || minus != other.minus;
    }

    [[nodiscard]] uint8_t get(int index) const {
        if (index < 0 || index >= 64) throw std::out_of_range("index out of range");
        if ((plus >> index) & 1ULL) return 1;
        if ((minus >> index) & 1ULL) return 2;
        return 0;
    }

    void set(int index, uint8_t value) {
        if (index < 0 || index >= 64) throw std::out_of_range("index out of range");
        value %= 3;
        uint64_t mask = ~(1ULL << index);
        plus &= mask;
        minus &= mask;
        if (value == 1) plus |= (1ULL << index);
        else if (value == 2) minus |= (1ULL << index);
    }

    [[nodiscard]] uint64_t support() const {
        return plus | minus;
    }
};


#endif //FFSYMMETRICFLIPS_ZTVec_H
