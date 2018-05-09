#ifndef DIV_H
#define DIV_H

#include "util.h"
#include <array>

class DivLut {
  using Denominator = int;
  constexpr static int lg_TableSize = 8;
  constexpr static int TableSize = 1 << lg_TableSize;

public:
  using Numerator = unsigned long;
  class Divisor;

  class Gen {
    friend class Divisor;
    using u128 = __uint128_t;
    u128 p;
    int lg_q;

    public:
    constexpr Gen(Numerator d) :
      p((((u128)1 << 127) - 1 + d) / d),
      lg_q(127) {}
    constexpr Gen(u128 p, int lg_q) : p(p), lg_q(lg_q) {}
    Gen operator/(uint64_t n) const {
      return Gen((p >> ceil_lg((u128)n)) * n, lg_q - ceil_lg((u128)n));
    }
  };
  class Divisor {
    friend class Gen;
    Numerator p;

    public:
    constexpr Divisor() : p(0) {}
    constexpr Divisor(Numerator d) : // ceiling 2^k/d
      p(1 == d ? (~0) : ((1UL << 63) - 1 + d) / d * 2) { } // k-64 = 0
    constexpr Divisor(Gen g) : p((uint64_t)(
        ceil_lg(g.p) - g.lg_q == 64 ? ~0UL
        : g.lg_q < 64 ? g.p << (64 - g.lg_q)
        : g.p >> (g.lg_q - 64))) {}

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
};
#endif
