#ifndef DIV_H
#define DIV_H

#include "util.h"
#include <array>

inline int lgl_flr(uint64_t x) {
  assert(x >= 1); // clz < 1 undefined
  return 64 - __builtin_clzll(x);
}

class DivLut {
  using Numerator = unsigned long;
  using Denominator = int;
  constexpr static int lg_TableSize = 8;
  constexpr static int TableSize = 1 << lg_TableSize;

public:
  class Divisor {
    Numerator p;

    public:
    constexpr Divisor() : p(0) {}
    constexpr Divisor(Numerator d) : // ceiling 2^k/d
      p(1 == d ? (~0) : ((1UL << 63) - 1 + d) / d * 2) { } // k-64 = 0
    friend Numerator operator/(__uint128_t n, const Divisor& d) { return (n * d.p) >> 64; }
    Divisor operator<<(int n) const {
      Divisor d = *this;
      d.p >>= n;
      return d;
    }
    Divisor operator/(Numerator n) const {
      Divisor d = *this;
      d.p *= n;
      return d;
    }
  };

private:
  std::array<Divisor, TableSize> t;

public:
  constexpr DivLut() { for (auto i = 1; i < TableSize; i++) t[i] = i; }
  Divisor operator[](const Numerator d) const {
    const Denominator k = lgl_flr(d) - lg_TableSize;
    auto r = t[d >> k] << k;
    return d > TableSize? r : t[d];
  }
  DivLut& operator<<=(int n) { for (auto& d : t) d = d << n; return *this;}
  DivLut& operator/=(Numerator n) { for (auto& d : t) d = d / n; return *this;}
};
#endif
