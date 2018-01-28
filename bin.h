#ifndef BIN_H
#define BIN_H

#include <cinttypes>
#include <assert.h>

#include "lin.h"                                                                                          
#include "util.h"                                                                                         
#include "padded_vector.h"                                                                                

#include <limits>

#if IACA == 1
#include <iacaMarks.h>
#else
#define IACA_START 
#define IACA_END 
#endif

template <bool RETURN_EARLY=true
         ,bool TEST_EQ=true
         ,bool FOR=false
         ,bool OVERFLOW_MATH=false
         ,int MIN_EQ_SZ=1
         ,typename Index = unsigned
         ,class Linear = LinearUnroll<>
         >
class BinaryLR {
  PaddedVector<> A;
  int lg_v;

  public:
  BinaryLR(const std::vector<Key>& _a, const std::vector<int>& indexes) : A(_a) {
    lg_v=0;
    for (auto n = A.size(); n > MIN_EQ_SZ; n -= (n/2)) lg_v++;
  }
    __attribute__((always_inline))
  Key operator()(const Key x) {
    // use pointer
    Index l = 0;
    Index r = A.size();

    for (int i =  0; FOR ? i < lg_v : r - l > 1; i++) {
      //IACA_START
      assert(l <= r);    // ordering check
      assert(l+r >= r); // overflow check
      Index m = OVERFLOW_MATH ? (l+r)/2 : l + (r-l) / 2 ;
      if (TEST_EQ) {
        if (A[m] < x) {
          l = m + 1;
        } else if (A[m] > x) {
          r = m;
        } else {
          if (RETURN_EARLY) return A[m];
          l = r = m;
        }
      } else {
        if (A[m] <= x) l = m;
        else r = m;
      }
    }
    if (MIN_EQ_SZ == 1) {
      //IACA_END
      return A[l];
    }

    Index guess = (l+r)/2;
    if (A[guess] < x) return A[Linear::forward(A.begin(), guess+1, x)];
    else return A[Linear::reverse(A.begin(), guess, x)];
  }
};

template <bool RETURN_EARLY=true
         ,bool TEST_EQ=true
         ,bool FOR=false
         ,bool POW_2=false
         ,int MIN_EQ_SZ=1
         ,typename Index = unsigned long
         ,class Linear = LinearUnroll<>
         >
class BinarySz {
  PaddedVector<> A;
  int lg_v, lg_min;

  public:
  BinarySz(const std::vector<Key>& _a, const std::vector<int>& indexes) : A(_a) {
    lg_v=lg_min=0;
    for (auto n = A.size();n > 1; n -= (n/2)) {
      lg_v++;
      if (n > MIN_EQ_SZ) lg_min++;
    }
  }

  __attribute__((always_inline))
    Key operator()(const Key x) {
      Index n = A.size();
      auto a = A.begin();
      Index left = 0L;
      if (POW_2) {
        Index mid = n - (1UL << (lg_v-1));
        left = a[mid] <= x ? mid : left;
        n -= mid;
      }
      // TODO how to set unroll = 4 if MIN_EQ_SZ > 1?, else 8?
#define LOOP \
        assert(left + n == A.size() || A[left + n] > x); \
        assert(A[left] <= x); \
        IACA_START \
        Index half = n / 2; \
        if (TEST_EQ) { \
          if (x < a[left + half]) { \
            n = half; \
          } else if (a[left + half] < x) { \
            left = left + half + 1; \
            n = n - half - 1; \
          } else { \
            if (RETURN_EARLY) return a[left + half]; \
            left += half; \
            if (FOR) i = lg_min; else n = 0; \
          } \
        } else { \
          left = a[left + half] <= x ? left + half : left; \
          if (POW_2) n = half; \
          else n -= half; \
        } 
      if (MIN_EQ_SZ == 1) {
#pragma unroll(8)
        for (int i = POW_2 ? 1 : 0; FOR ? i < lg_min : n > MIN_EQ_SZ; i++) {
          LOOP
        }
      } else {
#pragma unroll(4)
        for (int i = POW_2 ? 1 : 0; FOR ? i < lg_min : n > MIN_EQ_SZ; i++) {
          LOOP
        }
      }

      IACA_END
      if (MIN_EQ_SZ == 1) return a[left];

      Index guess = left + n/2;
      if (a[guess] < x) return a[Linear::forward(a,guess+1,x)];
      else return a[Linear::reverse(a,guess,x)];
    }
};

using BinaryLRCond = BinaryLR<false>;
using BinaryLROver = BinaryLR<false, true, false, true>;
using BinaryLRNoEq = BinaryLR<false, false, false, true>;
using BinaryLRFor = BinaryLR<false, true, true, true>;
using BinaryLRNoEqFor = BinaryLR<false, false, true, true>;
using BinaryLRLin = BinaryLR<false, false, true, true, 32>;

using BinarySzCond = BinarySz<false>;
using BinarySzNoEq = BinarySz<false, false>;
using BinarySzFor = BinarySz<false, true, true>;
using BinarySzNoEqFor = BinarySz<false, false, true>;
using BinarySzPow = BinarySz<false, false, true, true>;
using BinarySzLin = BinarySz<false, false, true, true, 32>;

using B2 = BinarySz<false, false, false, true>;
//using B0 = BinaryLinSize;
//using B1 = BinarySize<32, true, true, unsigned long, LinearUnroll<>>;

#endif
