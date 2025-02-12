#include <algorithm>
#include <cstdint>
#include <limits>
#include <numeric>
#include <random>
#include <vector>

#include <benchmark/benchmark.h>
using Key = int64_t;

template <int record_bytes = 8, int pad = 32> class PaddedVector {
  static constexpr int payload_bytes = record_bytes - sizeof(Key);
  using Payload = char[payload_bytes];
  struct Record {
    Key k;
    Payload p;
    Record() {}
    Record(Key k) : k(k) {}
    bool operator<(const Record &r) const { return k < r.k; }
    operator Key() const { return k; }
  };

  std::vector<Record> v;

public:
  PaddedVector(size_t n) : v(n + 2 * pad) {
    std::fill(this->v.begin(), this->v.begin() + pad,
              std::numeric_limits<Key>::min());
    std::fill(this->v.end() - pad, this->v.end(),
              std::numeric_limits<Key>::max());
  }
  PaddedVector(const std::vector<Key> &v) : v(v.size() + 2 * pad) {
    std::copy(v.begin(), v.end(), this->v.begin() + pad);
    std::fill(this->v.begin(), this->v.begin() + pad,
              std::numeric_limits<Key>::min());
    std::fill(this->v.end() - pad, this->v.end(),
              std::numeric_limits<Key>::max());
  }
  Key &operator[](long ix) {
    // allow some inaccuracy to reduce needed precision
    assert(ix >= -pad);
    assert(ix <= size() + pad);
    return v[ix + pad].k;
  }
  const Key &operator[](long ix) const {
    // allow some inaccuracy to reduce needed precision
    assert(ix >= -pad);
    assert(ix <= size() + pad);
    return v[ix + pad].k;
  }
  auto begin() { return v.begin() + pad; }
  auto end() { return v.end() - pad; };
  const Key *cbegin() const { return v.data() + pad; }
  size_t size() const { return v.size() - 2 * pad; }
  Key back() const { return (*this)[size() - 1]; }
  auto get_pad() const { return pad; }
};

template <int roll = 2> class LinearSIMD {
  template <bool reverse = false>
  static int64_t aligned(const int64_t *ptr, int64_t i, const int64_t x) {
    auto vecXX = reverse ? _mm256_set1_epi64x(x) : _mm256_set1_epi64x(x - 1);
    for (;; i = reverse ? (i - 16) : i + 16) {
      auto sign = reverse ? -1 : 1;
      auto av0 = _mm256_load_si256((__m256i *)(ptr + i + sign * 0));
      auto av1 = _mm256_load_si256((__m256i *)(ptr + i + sign * 4));
      auto av2 = _mm256_load_si256((__m256i *)(ptr + i + sign * 8));
      auto av3 = _mm256_load_si256((__m256i *)(ptr + i + sign * 12));
      auto cmp3 = reverse ? _mm256_cmpgt_epi64(vecXX, av3)
                          : _mm256_cmpgt_epi64(av3, vecXX);
      auto msk3 = _mm256_movemask_epi8(cmp3);
      if (!msk3)
        continue;
      auto cmp0 = reverse ? _mm256_cmpgt_epi64(vecXX, av0)
                          : _mm256_cmpgt_epi64(av0, vecXX);
      auto cmp1 = reverse ? _mm256_cmpgt_epi64(vecXX, av1)
                          : _mm256_cmpgt_epi64(av1, vecXX);
      auto cmp2 = reverse ? _mm256_cmpgt_epi64(vecXX, av2)
                          : _mm256_cmpgt_epi64(av2, vecXX);
      auto msk0 = _mm256_movemask_epi8(cmp0);
      auto msk1 = _mm256_movemask_epi8(cmp1);
      auto msk2 = _mm256_movemask_epi8(cmp2);
      if (msk0)
        return reverse ? (i + 4 - _lzcnt_u32(msk0) / 8 - 0 * 4)
                       : i + _tzcnt_u32(msk0) / 8 + 0 * 4;
      if (msk1)
        return reverse ? (i + 4 - _lzcnt_u32(msk1) / 8 - 1 * 4)
                       : i + _tzcnt_u32(msk1) / 8 + 1 * 4;
      if (msk2)
        return reverse ? (i + 4 - _lzcnt_u32(msk2) / 8 - 2 * 4)
                       : i + _tzcnt_u32(msk2) / 8 + 2 * 4;
      if (msk3)
        return reverse ? (i + 4 - _lzcnt_u32(msk3) / 8 - 3 * 4)
                       : i + _tzcnt_u32(msk3) / 8 + 3 * 4;
    }
  }

  template <bool reverse = false>
  static int64_t linSIMD(const int64_t *arr, const int64_t guessIx,
                         const int64_t x) {
    auto ptr = arr;
    auto i = guessIx;
    auto misalignment = ((uintptr_t)(ptr + i) & 31) / sizeof(int64_t);
    for (int j = 0; j < 4 * roll; j++)
      if (reverse ? (arr[i - j] <= x) : arr[i + j] >= x)
        return reverse ? i - j : i + j;
    i = reverse ? (i - 4 * (roll - 1) - misalignment)
                : i + 4 * roll - misalignment;
    return aligned<reverse>(arr, i, x);
  }

public:
  static int64_t forward(const int64_t *a, int len, const int64_t guessIx,
                         const int64_t x) {
    return linSIMD<false>(a, guessIx, x);
  }
  static int64_t reverse(const int64_t *a, int len, const int64_t guessIx,
                         const int64_t x) {
    return linSIMD<true>(a, guessIx, x);
  }
};

template <int n = 8> class LinearUnroll {
  template <bool reverse = false>
  static int64_t linUnroll(const Key *a, int64_t m, Key k) {
    for (;; m = (reverse ? m - n : m + n)) {
      for (int i = 0; i < n; i++) {
        // assert(m+i < 1032); assert((m-i) > -32);
        if (reverse ? (a[m - i] <= k) : (a[m + i] >= k)) {
          return reverse ? (m - i) : (m + i);
        }
      }
    }
  }

public:
  static int64_t forward(const Key *a, int len, const int64_t guessIx,
                         const Key x) {
    return linUnroll<false>(a, guessIx, x);
  }
  static int64_t reverse(const Key *a, int len, const int64_t guessIx,
                         const Key x) {
    return linUnroll<true>(a, guessIx, x);
  }
};

template <int n = 8> struct LinearNoSentinel {
  template <bool reverse = false>
  static int64_t linUnroll(const Key *a, int len, int64_t m, Key k) {
    for (; reverse ? m >= 0 : m < len; m = (reverse ? m - n : m + n)) {
      for (int i = 0; i < n && (reverse ? m - i >= 0 : m + i < len); i++) {
        if (reverse ? (a[m - i] <= k) : (a[m + i] >= k)) {
          return reverse ? (m - i) : (m + i);
        }
      }
    }
    return reverse ? std::numeric_limits<Key>::min()
                   : std::numeric_limits<Key>::max();
  }

public:
  static int64_t forward(const Key *a, int len, const int64_t guessIx,
                         const Key x) {
    return linUnroll<false>(a, len, guessIx, x);
  }
  static int64_t reverse(const Key *a, int len, const int64_t guessIx,
                         const Key x) {
    return linUnroll<true>(a, len, guessIx, x);
  }
};

template <typename Search>
static void BENCHMARK_search(benchmark::State &state) {
  const int search_len = state.range(0);
  PaddedVector v(state.range(1));
  std::iota(v.begin(), v.end(), 0);

  std::array<Key, 1024> keys_to_search;
  std::uniform_int_distribution<int> distribution(0, v.size() - search_len);
  std::default_random_engine rng(42);
  auto gen_num = [&distribution, &rng]() { return distribution(rng); };

  int i = 0;
  for (auto _ : state) {
    state.PauseTiming();
    std::generate(keys_to_search.begin(), keys_to_search.end(), gen_num);
    state.ResumeTiming();
    for (auto k : keys_to_search) {
      auto x = Search::forward(&v[0], v.size(), k - search_len, k);
      benchmark::DoNotOptimize(x);
    }
  }
}

BENCHMARK_TEMPLATE(BENCHMARK_search, LinearUnroll<>)
    ->RangeMultiplier(4)
    ->Ranges({{1, 32}, {4096, 16777216}});
BENCHMARK_TEMPLATE(BENCHMARK_search, LinearUnroll<1>)
    ->RangeMultiplier(4)
    ->Ranges({{1, 32}, {4096, 16777216}});
BENCHMARK_TEMPLATE(BENCHMARK_search, LinearNoSentinel<>)
    ->RangeMultiplier(4)
    ->Ranges({{1, 32}, {4096, 16777216}});
BENCHMARK_TEMPLATE(BENCHMARK_search, LinearNoSentinel<1>)
    ->RangeMultiplier(4)
    ->Ranges({{1, 32}, {4096, 16777216}});
BENCHMARK_TEMPLATE(BENCHMARK_search, LinearSIMD<>)
    ->RangeMultiplier(4)
    ->Ranges({{1, 32}, {4096, 16777216}});
BENCHMARK_TEMPLATE(BENCHMARK_search, LinearSIMD<1>)
    ->RangeMultiplier(4)
    ->Ranges({{1, 32}, {4096, 16777216}});
BENCHMARK_MAIN();
