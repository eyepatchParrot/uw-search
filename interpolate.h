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

class IBase {
public:
  using Index = int64_t;

  static constexpr int Recurse = -1;
  static constexpr bool Precompute = true;
  static constexpr bool Intercept = true;

  template <bool fold=false>
  struct Lut {
    // maybe want a.size() since we truncate
    Lut(const PaddedVector<>& a) : A(a), lgScale(std::max(0L, lg(A.size() - 1UL) + lg((uint64_t)A.back()) - 64L)) {
      if (fold) divisors /= (A.size() - 1);
      d_range_width = DivLut::Gen(A.back() - A[0]) / (A.size() - 1);
    }

    const PaddedVector<>& A;
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
  template <bool precompute=false>
    struct Float {
      Float(const PaddedVector<>& a) : A(a), f_aL(A[0]),
      f_width_range( (double)(A.size() - 1) / (double)(A.back() - A[0])) {}

      const PaddedVector<>& A;
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
    Exp(const PaddedVector<>& a) : A(a), y_1(A[A.size()>>1]
        ), diff_y_01(A[0] - y_1
          ),  a_0(diff_y_01 / (double)(y_1 - A[A.size()-1])
            ), lg_a_0(log(a_0)
              ), diff_scale(A[0] - a_0 * y_1
                ), d(A.size()>>1)
    {}

    const PaddedVector<>& A;
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



  struct IntDiv {
    IntDiv(const PaddedVector<>& a) : A(a),
    i_range_width((A.back() - A[0]) / (A.size() - 1)) {}

    const PaddedVector<>& A;
    Key i_range_width;

    Index operator()(const Key x, const Index left, const Index right) {
      return left + (x-A[left]) / ((A[right]-A[left]) / (right-left));
    }
    Index operator()(const Key x) {
      return (x - A[0]) / i_range_width;
    }
  };
protected:
  const PaddedVector<>& A;

  IBase(const PaddedVector<>& v) : A(v) {}

};

template <class Interpolate = IBase::Lut<>
         ,int nIter = 1
         ,class Linear = LinearUnroll<>
         ,int guardOff=0
         >
class Interpolation : public IBase {
  Interpolate interpolate;
  const Key* a0;

  __attribute__((always_inline))
  auto is(const Key x) {
    Index left = 0;
    Index right = A.size() - 1;
    assert(A.size() >= 1);

    Index mid = interpolate(x);
#ifndef NDEBUG
//    std::cout << mid << ' ';
#endif
    for (int i = 1; (nIter < 0 ? true : i < nIter); i++) {
      IACA_START
      if (A[mid] < x) left = mid+1;
      else if (A[mid] > x) right = mid-1;
      else return A[mid];
      if (left == right) return A[left];

      assert(left<right);
      assert(left >= 0); assert(right < A.size());
      mid = interpolate(x, left, right);
      assert(mid > -32); assert(mid < A.size()+32);
#ifndef NDEBUG
//      std::cout << "mid " <<mid << ' ';
//      auto fp = IBase::Float<IBase::Precompute>(A);
//      auto d2 = (int)fp(x, left, right) - mid;
//      if (d2*d2 > 1) {
//        printf("%lu, %lu, %lu = %lu ~ %ld\n", x, left, right, mid, d2);
//      }
#endif

      if (nIter < 0) { 
        IACA_END
        if (mid+guardOff >= right) {
          auto r = A[Linear::reverse(&A[0], right, x)];
//          err += i;
          return r;
        } else if (mid-guardOff <= left) {
          auto r = A[Linear::forward(&A[0], left, x)];
          err += i;
          return r;
        }
      }
      assert(mid >= left); assert(mid <= right);
    }

    if (A[mid] >= x) {
      auto r = Linear::reverse(&A[0], mid, x);
//      err += abs(r-mid);
      return A[r];
    } else {
      auto r = Linear::forward(&A[0], mid + 1, x);
//      err += abs(r-mid);
      return A[r];
    }
  }

  public:
  Interpolation(const PaddedVector<>& v) : IBase(v), interpolate(A), a0(A.cbegin()), err(0) { }

  long err;

  __attribute__((always_inline))
  Key operator()(const Key x) { return is(x); }
};

template<class Interpolate = IBase::Lut<> >
class InterpolateSlope : public IBase {
  static constexpr int guardOff=32;

  Interpolate interpolate;
  using Linear = LinearUnroll<>;

  public:
  InterpolateSlope(const PaddedVector<>& v) : IBase(v), interpolate(A) { }
  long err;
  __attribute__((always_inline))
  Key operator()(const Key x) {
    Index left = 0;
    Index right = A.size() - 1;
    assert(A.size() >= 1);

    Index mid = interpolate(x);
#ifndef NDEBUG
    std::cout << mid << ' ';
#endif
    for (int i = 1; ; i++) {
      IACA_START
      if (A[mid] < x) left = mid+1;
      else if (A[mid] > x) right = mid-1;
      else return A[mid];
      if (left == right) return A[left];

      assert(left<right);
      assert(left >= 0); assert(right < A.size());
      mid = interpolate(x, mid);

      assert(mid > -32); assert(mid < A.size()+32);
#ifndef NDEBUG
      std::cout << mid << ' ';
      auto fp = IBase::Float<IBase::Precompute>(A);
      auto d2 = (int)fp(x, left, right) - mid;
      if (d2*d2 > 1)
        printf("%lu, %lu, %lu = %lu ~ %ld\n", x, left, right, mid, d2);
#endif
      if (mid+guardOff >= right) return A[Linear::reverse(&A[0], right, x)];
      else if (mid-guardOff <= left) return A[Linear::forward(&A[0], left, x)];
      assert(mid >= left); assert(mid <= right);
    }

    return 0; 
  }
};

using InterpolationNaive = Interpolation<IBase::Float<>,IBase::Recurse, LinearUnroll<>>;
using InterpolationPrecompute = Interpolation<IBase::Float<IBase::Precompute>,IBase::Recurse,LinearUnroll<>>;
using InterpolationRecurseGuard = Interpolation<IBase::Float<IBase::Precompute>, IBase::Recurse, LinearUnroll<>, 32>;
using i_slope = InterpolateSlope<>;
using InterpolationRecurse3 = Interpolation<IBase::Float<IBase::Precompute>, 3>;
using InterpolationRecurseLut = Interpolation<IBase::Lut<>, IBase::Recurse, LinearUnroll<>, 32>;
using InterpolationLinearFp = Interpolation<IBase::Float<IBase::Precompute>>;
using InterpolationLinear = Interpolation<>;
using i_seq_intercept = Interpolation<IBase::Lut<IBase::Intercept>>;
using i_simd = Interpolation<IBase::Lut<>, 1, LinearSIMD<>>;
using i_exp_seq = Interpolation<IBase::Exp>;
using i_exp = Interpolation<IBase::Exp, IBase::Recurse>;
//using InterpolationIDiv = Interpolation<IBase::IntDiv> ;
#endif
