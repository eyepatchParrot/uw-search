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

template <int record_bytes=8
         ,bool RETURN_EARLY=true
         ,bool TEST_EQ=true
         ,bool FOR=false
         ,bool POW_2=false
         ,int MIN_EQ_SZ=1
         ,typename Index = unsigned long
         >
class Binary {
  using Vector = PaddedVector<record_bytes>;
  using Linear = LinearUnroll<Vector>;

  const Vector& A;
  int lg_v, lg_min;

  public:
  Binary(const Vector& _a) : A(_a) {
    lg_v=lg_min=0;
    for (auto n = A.size();n > 1; n -= (n/2)) {
      lg_v++;
      if (n > MIN_EQ_SZ) lg_min++;
    }
  }

  __attribute__((always_inline))
    Key operator()(const Key x) {
      Index n = A.size();
      Index left = 0L;
      if (POW_2) {
        Index mid = n - (1UL << (lg_v-1));
        left = A[mid] <= x ? mid : left;
        n -= mid;
      }
      // TODO how to set unroll = 4 if MIN_EQ_SZ > 1?, else 8?
#define LOOP \
        assert(left + n == A.size() || A[left + n] > x); \
        assert(A[left] <= x); \
        IACA_START \
        Index half = n / 2; \
        if (TEST_EQ) { \
          if (x < A[left + half]) { \
            n = half; \
          } else if (A[left + half] < x) { \
            left = left + half + 1; \
            n = n - half - 1; \
          } else { \
            if (RETURN_EARLY) return A[left + half]; \
            left += half; \
            if (FOR) i = lg_min; else n = 0; \
          } \
        } else { \
          left = A[left + half] <= x ? left + half : left; \
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
      if (MIN_EQ_SZ == 1) return A[left];

      Index guess = left + n/2;
      if (A[guess] < x) return A[Linear::forward(A, guess+1,x)];
      else return A[Linear::reverse(A,guess,x)];
    }
};

#define b_lin(PAYLOAD) Binary<PAYLOAD, false, false, true, true, 32>

#endif
