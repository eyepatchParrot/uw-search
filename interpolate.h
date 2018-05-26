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

template <int record_bytes> class IBase {
public:
  using Vector = PaddedVector<record_bytes>;
  using Index = int64_t;

  static constexpr int Recurse = -1;
  static constexpr bool Precompute = true;
  static constexpr bool Intercept = true;

  template <bool appx = true> struct Float {
    Float(const Vector &a)
        : A(a), slope(FixedPoint::Gen(A.size() - 1) / (A.back() - A[0])),
          f_aL(A[0]), f_width_range((double)((uint64_t)A.size() - 1) /
                                    (double)(A.back() - A[0])) {}

    const Vector &A;
    const FixedPoint slope;
    const double f_aL;
    const double f_width_range;

    Index operator()(const Key x, const Index left, const Index right) {
      return left + ((double)x - (double)(A[left])) /
                        (double)(A[right] - A[left]) * (double)(right - left);
    }

    Index operator()(const Key x, const Index mid) {
      return appx ? (x < A[mid] ? mid - slope * (uint64_t)(A[mid] - x)
                                : mid + slope * (uint64_t)(x - A[mid]))
                  : mid + (Index)(((double)x - (double)A[mid]) * f_width_range);
    }

    Index operator()(const Key x) {
      return appx ? slope * (uint64_t)(x - A[0])
                  : (Index)(((double)x - f_aL) * f_width_range);
    }
  };

  template <bool t = false, bool appx = false> struct Hyp3 {
    Hyp3(const Vector &a)
        : A(a), d((uint64_t)A.size() >> 1), y_1(A[d]), diff_y_01(A[0] - y_1),
          a_0(diff_y_01 == (y_1 - A.back())
                  ? 0.99999999999999
                  : diff_y_01 / (double)(y_1 - A.back())),
          diff_scale(A[0] - a_0 * A.back()), d_a((1.0 + a_0) * d) {}
    const Vector &A;
    const Index d;
    const Key y_1;
    const double diff_y_01, a_0, diff_scale, d_a;
    /*
     * def hyp3(points, y_star):
         x, y = zip(*[(x, y - y_star) for x, y in points])
         error = (x[1] - x[2]) * (x[1] - x[0]) * y[1] * (y[2] - y[0]) / (
             (x[1] - x[2]) * (y[0] - y[1]) * y[2] + (x[1] - x[0]) * (y[1] -
     y[2]) * y[0])
         return x[1] + error
     */
    Index operator()(const Key y, const Index x_0, const Index x_1,
                     const Index x_2) const {
      double y_0 = A[x_0] - y, y_1 = A[x_1] - y, y_2 = A[x_2] - y,
             error = y_1 * (x_1 - x_2) * (x_1 - x_0) * (y_2 - y_0) /
                     (y_2 * (x_1 - x_2) * (y_0 - y_1) +
                      y_0 * (x_1 - x_0) * (y_1 - y_2));
      return x_1 + (Index)error;
    }

    Index operator()(const Key x) {
      return d + (Index)(d_a * (y_1 - x) / (diff_scale - x * (a_0 + 1.0)));
    }
  };

  struct IntDiv {
    IntDiv(const Vector &a)
        : A(a), i_range_width((A.back() - A[0]) / ((uint64_t)A.size() - 1)) {}

    const Vector &A;
    Key i_range_width;

    Index operator()(const Key x, const Index left, const Index right) {
      return left + (x - A[left]) / ((A[right] - A[left]) / (right - left));
    }
    Index operator()(const Key x, const Index mid) {
      return mid + (x - A[mid]) / i_range_width;
    }
    Index operator()(const Key x) { return (x - A[0]) / i_range_width; }
  };

protected:
  const Vector &A;

  IBase(const Vector &v) : A(v) {}
};

template <int record_bytes, int guard_off,
          class Interpolate = typename IBase<record_bytes>::template Hyp3<> >
class i_hyp : public IBase<record_bytes> {
  using Super = IBase<record_bytes>;
  using Vector = typename Super::Vector;
  using typename Super::Index;
  using Super::A;
  using Linear = LinearUnroll<Vector>;
  static constexpr int nIter = Super::Recurse;
  static constexpr bool min_width = false;

  Interpolate interpolate;

  __attribute__((always_inline)) Key linear_search(const Key x, Index y) const {
    if (A[y] >= x) {
      return A[Linear::reverse(A, y, x)];
    } else {
      return A[Linear::forward(A, y + 1, x)];
    }
  }

public:
  i_hyp(const Vector &v) : Super(v), interpolate(A) { assert(A.size() >= 1); }

  __attribute__((always_inline)) Key operator()(const Key x) {
    Index left = 0, right = A.size() - 1, next_1 = A.size() >> 1,
          next_2 = interpolate(x);
    for (int i = 1; nIter < 0 || i < nIter; i++) {
      if (next_2 - next_1 <= guard_off && next_2 - next_1 >= -guard_off)
        return linear_search(x, next_2);
      assert(next_1 >= left);
      assert(next_1 <= right);
      assert(next_2 >= left);
      assert(next_2 <= right);
      assert(next_1 != next_2);
      if (next_1 < next_2) {
        assert(A[next_1] <= x); // f(x) <= f(x') ==> x <= x'
        left = next_1;
      } else {
        assert(A[next_1] >= x); // f(x) >= f(x') ==> x >= x'
        right = next_1;
      }
      if (next_2 + guard_off >= right) {
        auto r = A[Linear::reverse(A, right, x)];
        IACA_END
        return r;
      } else if (next_2 - guard_off <= left) {
        auto r = A[Linear::forward(A, left, x)];
        IACA_END
        return r;
      }
      next_1 = next_2;

      assert(left < right);
      assert(left >= 0);
      assert(right < A.size());
      assert(next_1 != left);
      assert(next_1 != right);

      next_2 = interpolate(x, left, next_1, right);

      assert(next_2 >= left);
      assert(next_2 <= right);
    }
    return linear_search(x, next_2);
  }
};

template <int record_bytes,
          class Interpolate = typename IBase<record_bytes>::template Float<>,
          int nIter = IBase<record_bytes>::Recurse, int guard_off = 16,
          bool min_width = false>
class Interpolation : public IBase<record_bytes> {
  using Super = IBase<record_bytes>;
  using Vector = typename Super::Vector;
  using typename Super::Index;
  using Super::A;
  using Linear = LinearUnroll<Vector>;

  Interpolate interpolate;

public:
  Interpolation(const Vector &v) : Super(v), interpolate(A) {}

  __attribute__((always_inline)) Key operator()(const Key x) {
    assert(A.size() >= 1);
    Index left = 0, right = A.size() - 1, next = interpolate(x);
    for (int i = 1; nIter < 0 ? true : i < nIter; i++) {
      IACA_START
      if (A[next] < x)
        left = next + 1;
      else if (A[next] > x)
        right = next - 1;
      else
        return A[next];
      if (min_width) {
        if (right - left <= guard_off)
          return A[Linear::reverse(A, right, x)];
      } else if (left == right)
        return A[left];

      assert(left < right);
      assert(left >= 0);
      assert(right < A.size());
      next = interpolate(x, left, right);
      assert(next > -32);
      assert(next < A.size() + 32);

      if (guard_off >= 0 && next + guard_off >= right) {
        auto r = A[Linear::reverse(A, right, x)];
        IACA_END
        return r;
      } else if (guard_off >= 0 && next - guard_off <= left) {
        auto r = A[Linear::forward(A, left, x)];
        IACA_END
        return r;
      }
      assert(next >= left);
      assert(next <= right);
    }
    // linear search base case
    if (A[next] >= x) {
      return A[Linear::reverse(A, next, x)];
    } else {
      return A[Linear::forward(A, next + 1, x)];
    }
  }
};

template <int record_bytes, int nIter = IBase<record_bytes>::Recurse,
          class Interpolate = typename IBase<record_bytes>::template Float<>,
          int guard_off = 8>
class InterpolationSlope : public IBase<record_bytes> {
  using Super = IBase<record_bytes>;
  using Vector = typename Super::Vector;
  using typename Super::Index;
  using Super::A;
  using Super::Recurse;
  using Linear = LinearUnroll<Vector>;

  Interpolate interpolate;

public:
  InterpolationSlope(const Vector &v) : Super(v), interpolate(A) {}

  // TODO replace with flatten?
  __attribute__((always_inline)) Key operator()(const Key x) {
    assert(A.size() >= 1);
    // set bounds and do first interpolation
    Index left = 0, right = A.size() - 1, next = interpolate(x);
    for (int i = 1; (nIter < 0 ? true : i < nIter); i++) {
      IACA_START

      // update bounds and check for match
      if (A[next] < x)
        left = next + 1;
      else if (A[next] > x)
        right = next - 1;
      else
        return A[next];
      if (left == right)
        return A[left];

      // next interpolation
      assert(left < right);
      assert(left >= 0);
      assert(right < A.size());
      next = interpolate(x, next);

      if (nIter == Recurse) {
        // apply guards
        if (guard_off == -1)
          next = std::min(std::max(left, next), right);
        else {
          if (next + guard_off >= right)
            return A[Linear::reverse(A, right, x)];
          else if (next - guard_off <= left)
            return A[Linear::forward(A, left, x)];
        }
        assert(next >= left);
        assert(next <= right);
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
template <int record_bytes>
using i_naive =
    Interpolation<record_bytes,
                  typename IBase<record_bytes>::template Float<false>,
                  IBase<record_bytes>::Recurse, -1>;
template <int RECORD, int GUARD>
using i_opt =
    InterpolationSlope<RECORD, IBase<RECORD>::Recurse,
                       typename IBase<RECORD>::template Float<>, GUARD>;
template <int RECORD> using i_seq = InterpolationSlope<RECORD, 1>;
template <int RECORD, int GUARD>
using i_no_reuse =
    Interpolation<RECORD, typename IBase<RECORD>::template Float<>,
                  IBase<RECORD>::Recurse, GUARD>;
template <int RECORD>
using i_no_guard =
    InterpolationSlope<RECORD, IBase<RECORD>::Recurse,
                       typename IBase<RECORD>::template Float<>, -1>;
template <int RECORD>
using i_fp = InterpolationSlope<RECORD, IBase<RECORD>::Recurse,
                                typename IBase<RECORD>::template Float<false> >;
template <int RECORD>
using i_idiv = InterpolationSlope<RECORD, IBase<RECORD>::Recurse,
                                  typename IBase<RECORD>::IntDiv>;

#endif
