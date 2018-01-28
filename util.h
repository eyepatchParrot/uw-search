#ifndef UTIL_H
#define UTIL_H

#include <cinttypes>
#include <assert.h>

using Key = int64_t;
using SearchFn = int64_t(const Key*, int64_t, Key);

#ifdef NDEBUG
#define ERR(A, args...)
#else
#define ERR(A, args...) fprintf(stderr, A, args)
#endif

#ifndef N_RUNS
#ifndef NDEBUG
#define N_RUNS 5
#else
#define N_RUNS 5000
#endif
#endif

#ifndef N_SAMPLES
#define N_SAMPLES 10
#endif

#ifndef SUBSET_SIZE
#define SUBSET_SIZE -1
#endif

inline unsigned lg(unsigned x) {
  assert(x >= 2); // subtracting and clz < 1 is undefined.
  return 32 - __builtin_clz(x-1);
}

inline unsigned lg_flr(unsigned x) {
  assert(x >= 1);
  return 32 - __builtin_clz(x);
}

#endif
