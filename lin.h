#ifndef LIN_H
#define LIN_H

#include <x86intrin.h>
//#include <iacaMarks.h>

// TODO use padded vector
template <int roll=2>
class LinearSIMD {
  template <bool reverse=false>
  static int64_t aligned(const int64_t *ptr, int64_t i, const int64_t x) {
      auto vecXX = reverse? _mm256_set1_epi64x(x): _mm256_set1_epi64x(x-1);
      for (;;i = reverse?(i-16) : i+16) {
        auto sign = reverse?-1:1;
        auto av0 = _mm256_load_si256((__m256i*)(ptr + i + sign*0));
        auto av1 = _mm256_load_si256((__m256i*)(ptr + i + sign*4));
        auto av2 = _mm256_load_si256((__m256i*)(ptr + i + sign*8));
        auto av3 = _mm256_load_si256((__m256i*)(ptr + i + sign*12));
        auto cmp3 = reverse? _mm256_cmpgt_epi64(vecXX, av3) :  _mm256_cmpgt_epi64(av3, vecXX);
        auto msk3 = _mm256_movemask_epi8(cmp3);
        if (!msk3) continue;
        auto cmp0 = reverse? _mm256_cmpgt_epi64(vecXX, av0) :   _mm256_cmpgt_epi64(av0,vecXX );
        auto cmp1 = reverse? _mm256_cmpgt_epi64(vecXX, av1) :   _mm256_cmpgt_epi64(av1,vecXX );
        auto cmp2 = reverse? _mm256_cmpgt_epi64(vecXX, av2) :   _mm256_cmpgt_epi64(av2,vecXX );
        auto msk0 = _mm256_movemask_epi8(cmp0);
        auto msk1 = _mm256_movemask_epi8(cmp1);
        auto msk2 = _mm256_movemask_epi8(cmp2);
        if (msk0) return reverse? (i + 4 - _lzcnt_u32(msk0) / 8 - 0 * 4) : i + _tzcnt_u32(msk0) / 8 + 0 * 4;
        if (msk1) return reverse? (i + 4 - _lzcnt_u32(msk1) / 8 - 1 * 4) : i + _tzcnt_u32(msk1) / 8 + 1 * 4;
        if (msk2) return reverse? (i + 4 - _lzcnt_u32(msk2) / 8 - 2 * 4) : i + _tzcnt_u32(msk2) / 8 + 2 * 4;
        if (msk3) return reverse? (i + 4 - _lzcnt_u32(msk3) / 8 - 3 * 4) : i + _tzcnt_u32(msk3) / 8 + 3 * 4;
      }
  }

  template <bool reverse=false>
    static int64_t linSIMD(const int64_t* arr, const int64_t guessIx, const int64_t x) {
      auto ptr = arr;
      auto i = guessIx;
      auto misalignment = ((uintptr_t)(ptr+i) & 31)/sizeof(int64_t);
      for (int j = 0; j < 4*roll; j++)
        if (reverse? (arr[i-j] <= x) : arr[i+j] >= x) return reverse? i-j : i+j;
      i = reverse? (i-4*(roll-1) - misalignment) : i + 4*roll - misalignment;
      return aligned<reverse>(arr, i, x);
    }
public:
    static int64_t forward(const int64_t* a, const int64_t guessIx, const int64_t x) {
      return linSIMD<false>(a, guessIx, x); 
    }
    static int64_t reverse(const int64_t* a, const int64_t guessIx, const int64_t x) {
      return linSIMD<true>(a, guessIx, x);
    }
};

template <int n=8>
class LinearUnroll {
  template <bool reverse=false>
    static int64_t linUnroll(const int64_t* a, int64_t m, int64_t k) {
      for (;;m = (reverse?m-n:m+n)) {
        for (int i = 0; i < n; i++) {
          //assert(m+i < 1032); assert((m-i) > -32);
          if (reverse?(a[m-i]<=k):(a[m+i]>=k)){
            return reverse?(m-i):(m+i);
          }
        }
      }
    }
public:
    static int64_t forward(const int64_t* a, const int64_t guessIx, const int64_t x) {
      return linUnroll<false>(a, guessIx, x); 
    }
    static int64_t reverse(const int64_t* a, const int64_t guessIx, const int64_t x) {
      return linUnroll<true>(a, guessIx, x);
    }
};



#endif
