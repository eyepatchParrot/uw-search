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
    Lut(const PaddedVector<>& a) : A(a), lgScale(lg(A.size() - 1)) {
      if (fold) {
        divisors /= (A.size() - 1);
        //divisors <<= lgScale;
      }
      d_range_width = (DivLut::Divisor((A.back() - A[0]) >>  lgScale) << lgScale) / (A.size() - 1);
    }

    const PaddedVector<>& A;
    int lgScale;
    DivLut::Divisor d_range_width;
    DivLut divisors;

    Index operator()(const Key x, const Index left, const Index right) {
      if (fold) return left + (Key)((x - A[left]) >> lgScale) / divisors[(A[right] - A[left])];
      return left + (Key)(((x - A[left]) >> lgScale) * (right-left)) /
        divisors[(A[right] - A[left]) >> lgScale];
    }
    Index operator()(const Key x) { return (x-A[0]) / d_range_width; }
    Index operator()(const Key x, const Index mid) {
      return mid + (x - A[mid]) / d_range_width;
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
#ifndef NDEBUG
      std::cout << mid << ' ';
      auto fp = IBase::Float<IBase::Precompute>(A);
      auto d2 = (int)fp(x, left, right) - mid;
      if (d2*d2 > 1) {
        printf("%lu, %lu, %lu = %lu ~ %ld\n", x, left, right, mid, d2);
      }
#endif

      if (nIter < 0) { 
        IACA_END
        if (mid+guardOff >= right) return A[Linear::reverse(&A[0], right, x)];
        else if (mid-guardOff <= left) return A[Linear::forward(&A[0], left, x)];
      }
      assert(mid >= left); assert(mid <= right);
    }

    if (A[mid] > x) {
    //IACA_END
      auto r = A[Linear::reverse(&A[0], mid - 1, x)];
      return r;
    } else {
    //IACA_END
      auto r = A[Linear::forward(&A[0], mid, x)];
      return r;
    }
  }

  public:
  Interpolation(const PaddedVector<>& v) : IBase(v), interpolate(A), a0(A.cbegin()) { }

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
//using InterpolationIDiv = Interpolation<IBase::IntDiv> ;
#endif
