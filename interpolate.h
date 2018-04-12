#ifndef INTERPOLATE_H
#define INTERPOLATE_H

#include "div.h"
#include "lin.h"
#include "padded_vector.h"

#include <vector>
#include <algorithm>
#include <array>
#include <iostream>
#include <iterator>
#include <cmath>
#include <tuple>
#include <cmath>

#if IACA == 1
#include <iacaMarks.h>
#else
#define IACA_START 
#define IACA_END 
#endif

template <int record_bytes>
class IBase {
public:
  using Vector = PaddedVector<record_bytes>;
  using Index = int64_t;

  static constexpr int Recurse = -1;
  static constexpr bool Precompute = true;
  static constexpr bool Intercept = true;

  template <bool fold=false>
  struct Lut {
    // maybe want a.size() since we truncate
    Lut(const Vector& a) : A(a), lgScale(std::max(0L, lg(A.size() - 1UL) + lg((uint64_t)A.back()) - 64L)) {
      if (fold) divisors /= (A.size() - 1);
      d_range_width = DivLut::Gen(A.back() - A[0]) / (A.size() - 1);
    }

    const Vector& A;
    int lgScale;
    DivLut::Divisor d_range_width;
    DivLut divisors;

    Index operator()(const Key x, const Index left, const Index right) {
      if (fold) return left + (Key)((x - A[left]) >> lgScale) / divisors[(A[right] - A[left])];
      return left + (uint64_t)(((x - A[left]) >> lgScale) * (right-left)) /
        divisors[(A[right] - A[left]) >> lgScale];
    }
    Index operator()(const Key x) { return (uint64_t)(x-A[0]) / d_range_width; }
    Index operator()(const Key x, const Index mid) {
      return mid + (Key)(x - A[mid]) / d_range_width;
    }
  };
  template <bool precompute=true>
    struct Float {
      Float(const Vector& a) : A(a), f_aL(A[0]),
      f_width_range( (double)(A.size() - 1) / (double)(A.back() - A[0])) {}

      const Vector& A;
      const double f_aL;
      const double f_width_range;

      Index operator()(const Key x, const Index left, const Index right) {
        return left + ((double)x - (double)(A[left])) /
          (double)(A[right] - A[left]) * (double)(right-left);
      }

      Index operator()(const double x, const Index mid) {
        return mid + (x - (double)A[mid]) * f_width_range;
      }

      Index operator()(const Key x) {
        return precompute? (Index)(((double)x - f_aL) * f_width_range) :
          (*this)(x, 0, A.size()-1);
      }
    };

  struct Exp {
    Exp(const Vector& a) : A(a), y_1(A[A.size()>>1]
        ), diff_y_01(A[0] - y_1
          ),  a_0(diff_y_01 == (y_1 - A[A.size()-1]) ? 0.99999999999999 :
            diff_y_01 / (double)(y_1 - A[A.size()-1])
            ), lg_a_0(log(a_0)
              ), diff_scale(A[0] - a_0 * y_1
                ), d(A.size()>>1)
    {}

    const Vector& A;
    const double y_1, diff_y_01, a_0, lg_a_0, diff_scale, d;
/*
def exp_next(points, y_star):
    x, y = zip(*points)
    y = [i - y_star for i in y]
    d = [x[1] - x[0]]
    assert((d[0] - x[2] + x[1])**2<4)
    a = (y[0] - y[1]) / (y[1] - y[2])
    b = (y[0] - y[1]) / (y[0] - a * y[1])
    return x[1] + d[0]  * log(b) / log(a)
*/
    Index operator()(const Key x, const Index x_0, const Index x_2) {
      const Index x_1 = (x_0 + x_2) >> 1;
      auto d = x_1 - x_0;
      double y_0 = A[x_0] - x, y_1 = A[x_1] - x, y_2 = A[x_2] - x;
      auto a = (y_0 - y_1) / (y_1 - y_2),
           b = (y_0 - y_1) / (y_0 - a * y_1);
      //std::cerr << a << ' ' << log2(a) << ' ' << b << ' ' << log2(b) << '\n';
      //std::cerr << d * log(b) / log(a) << ' ' << d * floor(log2(b)) / floor(log2(a)) << '\n';
      return x_1 + d * log(b) / log(a);
    }

    Index operator()(const Key x) {
      auto b_0 = diff_y_01 / (diff_scale + a_0 * x - x); 
      return d + d * log(b_0) / lg_a_0;
    }
  };

  struct Hyp {
    Hyp(const Vector& a) : A(a), y_1(A[A.size()>>1]
        ), diff_y_01(A[0] - y_1
          ),  a_0(diff_y_01 == (y_1 - A[A.size()-1]) ? 0.99999999999999 :
            diff_y_01 / (double)(y_1 - A[A.size()-1])
            ), diff_scale(A[0] - a_0 * A[A.size()-1]
              ), d(A.size()>>1
                ), d_a((1.0+a_0) * d)
    {}

    const Vector& A;
    const double y_1, diff_y_01, a_0, diff_scale, d, d_a;
/*
def hyp_next(points, y_star):
    x, y = zip(*points)
    y = [i - y_star for i in y]
    d = [x[1] - x[0]]
    assert((d[0] - x[2] + x[1])**2<4)
    a = (y[0] - y[1]) / (y[1] - y[2])
    b = (y[0] - y[1]) / (y[0] - a * y[1])
    return x[1] + d[0] * y[1] * (a + 1) / (y[0] - a * y[2])
*/
    Index operator()(const Key x, const Index x_0, const Index x_2) {
      const Index x_1 = (x_0 + x_2) >> 1;
      auto d = x_1 - x_0;
      double y_0 = A[x_0] - x, y_1 = A[x_1] - x, y_2 = A[x_2] - x;
      auto a = (y_0 - y_1) / (y_1 - y_2);
      return x_1 + d * y_1 * (a + 1.0) / (y_0 - a * y_2);
    }

    Index operator()(const Key x) {
      return d + d_a * (y_1 - x) / (diff_scale - x * (a_0 + 1.0));
    }
  };

  struct Hyp3 {
    Hyp3(const Vector& a) : A(a), y_1(A[A.size()>>1]
        ), diff_y_01(A[0] - y_1
          ),  a_0(diff_y_01 == (y_1 - A[A.size()-1]) ? 0.99999999999999 :
            diff_y_01 / (double)(y_1 - A[A.size()-1])
            ), diff_scale(A[0] - a_0 * A[A.size()-1]
              ), d(A.size()>>1
                ), d_a((1.0+a_0) * d)
    {}
    const Vector& A;
    const double y_1, diff_y_01, a_0, diff_scale, d, d_a;
    /*
     * def hyp3(points, y_star):
         x, y = zip(*[(x, y - y_star) for x, y in points])
         error = (x[1] - x[2]) * (x[1] - x[0]) * y[1] * (y[2] - y[0]) / (
             (x[1] - x[2]) * (y[0] - y[1]) * y[2] + (x[1] - x[0]) * (y[1] - y[2]) * y[0])
         return x[1] + error
     */
    Index operator()(const Key y, const Index x_0, const Index x_1,
        const Index x_2) const {
      double y_0 = A[x_0] - y, y_1 = A[x_1] - y, y_2 = A[x_2] - y,
        error = y_1 * (x_1 - x_2) * (x_1 - x_0) * (y_2 - y_0) /
          (y_2 * (x_1 - x_2) * (y_0 - y_1) + y_0 * (x_1 - x_0) * (y_1 - y_2));
      return x_1 + (Index)error;
    }

    Index operator()(const Key x) {
      return d + d_a * (y_1 - x) / (diff_scale - x * (a_0 + 1.0));
    }
  };

  struct IntDiv {
    IntDiv(const Vector& a) : A(a),
    i_range_width((A.back() - A[0]) / (A.size() - 1)) {}

    const Vector& A;
    Key i_range_width;

    Index operator()(const Key x, const Index left, const Index right) {
      return left + (x-A[left]) / ((A[right]-A[left]) / (right-left));
    }
    Index operator()(const Key x, const Index mid) {
      return mid + (x - A[mid]) / i_range_width;
    }
    Index operator()(const Key x) {
      return (x - A[0]) / i_range_width;
    }
  };

protected:
  const Vector& A;

  IBase(const Vector& v) : A(v) {}
};

template <int record_bytes>
class Interpolation3 : public IBase<record_bytes> {
  using Super = IBase<record_bytes>;
  using Vector = typename Super::Vector;
  using typename Super::Index;
  using Super::A;
  using Interpolate = typename Super::Hyp3;
  using Linear = LinearUnroll<Vector>;
  static constexpr int nIter = Super::Recurse;
  static constexpr int guardOff=64;
  static constexpr bool min_width = false;

  Interpolate interpolate;

  __attribute__((always_inline))
  Key linear_search(const Key x, Index y) const {
    if (A[y] >= x) {
      return A[Linear::reverse(A, y, x)];
    } else {
      return A[Linear::forward(A, y + 1, x)];
    }
  }

  public:
  Interpolation3(const Vector& v) : Super(v), interpolate(A) { 
    assert(A.size() >= 1);
  }

  __attribute__((always_inline))
  Key operator()(const Key x) {
    Index left = 0, right = A.size() - 1, next_1 = A.size()>>1,
          next_2 = interpolate(x);
    for (int i = 1; nIter < 0 || i < nIter; i++) {
      if (next_2 - next_1 < guardOff && next_2 - next_1 > -guardOff) return linear_search(x, next_2);
      assert(next_1 >= left); assert(next_1 <= right); assert(next_2 >= left);
      assert(next_2 <= right); assert(next_1 != next_2);
      if (next_1 < next_2) {
        assert(A[next_1] <= x); // f(x) <= f(x') ==> x <= x'
        left = next_1;
      } else {
        assert(A[next_1] >= x); // f(x) >= f(x') ==> x >= x'
        right = next_1;
      }
      if (next_2+guardOff >= right) {
        auto r = A[Linear::reverse(A, right, x)];
        IACA_END
        return r;
      } else if (next_2-guardOff <= left) {
        auto r = A[Linear::forward(A, left, x)];
        IACA_END
        return r;
      }
      next_1 = next_2;

      assert(left<right); assert(left >= 0); assert(right < A.size());
      assert(next_1 != left); assert(next_1 != right);

      next_2 = interpolate(x, left, next_1, right);

      assert(next_2 >= left); assert(next_2 <= right);
    }
    return linear_search(x, next_2);
  }
};

template <int record_bytes
         ,class Interpolate = typename IBase<record_bytes>::template Lut<>
         ,int nIter = IBase<record_bytes>::Recurse
         ,int guardOff=32
         ,bool min_width = false
         >
class Interpolation : public IBase<record_bytes> {
  using Super = IBase<record_bytes>;
  using Vector = typename Super::Vector;
  using typename Super::Index;
  using Super::A;
  using Linear = LinearUnroll<Vector>;

  Interpolate interpolate;

  public:
  Interpolation(const Vector& v) : Super(v), interpolate(A) { }

  __attribute__((always_inline))
  Key operator()(const Key x) {
    assert(A.size() >= 1);
    Index left = 0, right = A.size() - 1, next = interpolate(x);
    for (int i = 1; nIter < 0 ? true : i < nIter; i++) {
      IACA_START
      if (A[next] < x) left = next+1;
      else if (A[next] > x) right = next-1;
      else return A[next];
      if (min_width) {
        if (right - left <= guardOff) return A[Linear::reverse(A, right, x)];
      } else if (left == right) return A[left];

      assert(left<right);
      assert(left >= 0); assert(right < A.size());
      next = interpolate(x, left, right);
      assert(next > -32); assert(next < A.size()+32);

      if (guardOff >= 0 && next+guardOff >= right) {
        auto r = A[Linear::reverse(A, right, x)];
        IACA_END
        return r;
      } else if (guardOff >= 0 && next-guardOff <= left) {
        auto r = A[Linear::forward(A, left, x)];
        IACA_END
        return r;
      }
      assert(next >= left); assert(next <= right);
    }
    // linear search base case
    if (A[next] >= x) {
      return A[Linear::reverse(A, next, x)];
    } else {
      return A[Linear::forward(A, next + 1, x)];
    }
  }
};

template <int record_bytes
         ,int nIter = IBase<record_bytes>::Recurse
         ,class Interpolate = typename IBase<record_bytes>::template Lut<>
         ,int guardOff=32
         >
class InterpolationSlope : public IBase<record_bytes> {
  using Super = IBase<record_bytes>;
  using Vector = typename Super::Vector;
  using typename Super::Index;
  using Super::A;
  using Super::Recurse;
  using Linear = LinearUnroll<Vector>;

  Interpolate interpolate;

  public:
  InterpolationSlope(const Vector& v) : Super(v), interpolate(A) { }

  // TODO replace with flatten?
  __attribute__((always_inline))
  Key operator()(const Key x) {
    assert(A.size() >= 1);
    // set bounds and do first interpolation
    Index left = 0, right = A.size() - 1, next = interpolate(x);
    for (int i = 1; (nIter < 0 ? true : i < nIter); i++) {
      IACA_START

      // update bounds and check for match
      if (A[next] < x) left = next+1;
      else if (A[next] > x) right = next-1;
      else return A[next];
      if (left == right) return A[left];

      // next interpolation
      assert(left<right);
      assert(left >= 0); assert(right < A.size());
      next = interpolate(x, next);

      assert(next > -A.get_pad()); assert(next < A.size()+A.get_pad());
      if (nIter == Recurse) { 
        // apply guards
        if (guardOff == -1) next = std::min(std::max(left, next), right);
        else {
          if (next+guardOff >= right) return A[Linear::reverse(A, right, x)];
          else if (next-guardOff <= left) return A[Linear::forward(A, left, x)];
        }
        assert(next >= left); assert(next <= right);
      }
    }
    // linear search base case
    if (A[next] >= x) {
      return A[Linear::reverse(A, next, x)];
    } else {
      return A[Linear::forward(A, next + 1, x)];
    }

    return 0; 
  }
};

/*
 * i_naive : Naive-IS
 * i_opt : Optimized-IS
 * i_seq : Interpolation-Sequential
 * i_recompute : Don't re-use slope
 * i_no_guard : guard = 0
 * i_fp : use FP division
 * i_idiv : use int division
 */
#define i_naive(RECORD) Interpolation<RECORD, typename IBase<RECORD>::template Float<false>, IBase<RECORD>::Recurse, -1>
#define i_opt(RECORD) InterpolationSlope<RECORD>
#define i_seq(RECORD) InterpolationSlope<RECORD, 1>
#define i_recompute(RECORD) Interpolation<RECORD>
#define i_no_guard(RECORD) InterpolationSlope<RECORD, IBase<RECORD>::Recurse, typename IBase<RECORD>::template Lut<>, -1>
#define i_fp(RECORD) InterpolationSlope<RECORD, IBase<RECORD>::Recurse, typename IBase<RECORD>::template Float<>>
#define i_idiv(RECORD) InterpolationSlope<RECORD, IBase<RECORD>::Recurse, typename IBase<RECORD>::IntDiv>
#define i_hyp(RECORD) Interpolation3<RECORD>
#endif
