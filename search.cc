#include "benchmark.h"
#include "util.h"

#include <algorithm>
#include <assert.h>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <numeric>
#include <sstream>
#include <vector>
#include <x86intrin.h>
#include <unordered_map>

std::vector<Key> randNums(int n, long seed) {
  std::mt19937_64 rng(seed);
  std::uniform_int_distribution<Key> dist(1, (1L<<63) - 2);
  std::vector<Key> v(n);
  for (auto& x : v) x = dist(rng);
  return v;
}

struct Input {
  private:
    struct KeyRandomSortedIndex {
      Key key;
      int random_index;
      int sorted_index;

      KeyRandomSortedIndex(Key k) { key = k; }
    };
  public:

  std::vector<Key> keys;
  std::vector<int> indexes;

  Input(std::vector<Key> nums, int n_gets) : indexes(n_gets) {
#ifdef DUMP_INPUT
    keys = nums;
    return;
#endif

    std::vector<KeyRandomSortedIndex> key_indexes(nums.begin(), nums.end());
    for (int i = 0; i < key_indexes.size(); i++)
      key_indexes[i].random_index = i;
    std::sort(key_indexes.begin(), key_indexes.end(), [](auto l, auto r) {
        return l.key < r.key; });
    for (int i = 0; i < key_indexes.size(); i++)
      key_indexes[i].sorted_index = i;

    // save the sorted input while it's here
    keys.resize(nums.size());
    for (auto i = 0; i < nums.size(); i++) keys[i] = key_indexes[i].key;

    std::sort(key_indexes.begin(), key_indexes.end(), [](auto l, auto r) {
        return l.random_index < r.random_index; });
    for (int i = 0; i < n_gets; i++) indexes[i] = key_indexes[i].sorted_index;
  }
};

// ./x datasetSz seed nThds [benchmarks ...]
int main(int argc, char *argv[]) {
  if (3 >= argc) return -1;

  using std::istream_iterator;

  int argi = 1;
  int datasetSz = atoi(argv[argi++]);
  long seed = atol(argv[argi++]);
  int n_threads = atoi(argv[argi++]);
  std::vector<std::string> benchmarks;
  while (argi < argc) benchmarks.push_back(argv[argi++]);

  std::cerr << datasetSz << ' ' << seed << '\n';

  int nGets = SUBSET_SIZE < 1 ? datasetSz : SUBSET_SIZE;

  Input input(randNums(datasetSz, seed), nGets);

#ifdef DUMP_INPUT
  for (auto x : input.keys) std::cout << x << '\n';
  return 0;
#endif

  std::vector<TestStats> tests;
  for (auto& benchmark_name : benchmarks)
    tests.push_back(benchmark(
          input.keys, input.indexes, n_threads, benchmark_name, seed));

  // print report
  for (int i = 0; i < benchmarks.size(); i++)
    std::cout << (0 != i ? "," : "") << benchmarks[i];
  std::cout << '\n';

  constexpr int n_samples = N_SAMPLES < 0 ? N_RUNS : std::min(N_RUNS, N_SAMPLES);
#ifndef NSORT
  for (auto& t : tests)
    std::partial_sort(t.ns.begin(), t.ns.begin() + n_samples, t.ns.end());
#endif
  for (int runIx = 0; runIx < n_samples; runIx++, std::cout << '\n')
    for (int testIx = 0; testIx < tests.size(); testIx++)
      printf("%s%.3f", 0 != testIx ? "," : "", tests[testIx].ns[runIx]);
}
