#ifndef UTIL_H
#define UTIL_H

#include <cinttypes>
#include <assert.h>
#include <vector>
#include <string>
#include <sstream>

using Key = int64_t;
using SearchFn = int64_t(const Key*, int64_t, Key);

#ifdef NDEBUG
#define ERR(A, args...)
#else
#define ERR(A, args...) fprintf(stderr, A, args)
#endif

#ifndef N_SAMPLES
#define N_SAMPLES 10
#endif

#ifndef SUBSET_SIZE
#define SUBSET_SIZE -1
#endif

constexpr inline unsigned lg(unsigned x) {
  assert(x >= 2); // subtracting and clz < 1 is undefined.
  return 32 - __builtin_clz(x-1);
}

constexpr inline long lg(uint64_t x) {
  assert(x >= 2);
  return 64 - __builtin_clzl(x);
}

constexpr inline long lg(__uint128_t x) {
  uint64_t hi = x >> 64;
  return x < 2 ? 1
    : hi != 0 ? 64 + lg(hi) : lg((uint64_t)x);
}

constexpr inline unsigned lg_flr(unsigned x) {
  assert(x >= 1);
  return 32 - __builtin_clz(x);
}

constexpr inline int lg_flrl(uint64_t x) {
  assert(x >= 1); // clz < 1 undefined
  return 64 - __builtin_clzl(x);
}

std::vector<std::string> split(std::string s, char delim) {
  std::replace(s.begin(), s.end(), ',', ' ');
  std::vector<std::string> v;
  std::stringstream ss(s);
  for (std::stringstream ss(s); ss.good(); ) {
    std::string s;
    ss >> s;
    v.push_back(s);
  }
  return v;
}

#endif
