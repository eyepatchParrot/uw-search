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

inline unsigned lg(unsigned x) {
  assert(x >= 2); // subtracting and clz < 1 is undefined.
  return 32 - __builtin_clz(x-1);
}

inline unsigned lg_flr(unsigned x) {
  assert(x >= 1);
  return 32 - __builtin_clz(x);
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
