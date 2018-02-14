#ifndef BENCHMARK_H
#define BENCHMARK_H

#include "oracle.h"
#include "interpolate.h"
#include "bin.h"
#include "lin.h"
//#include "fib.h"

#include "omp.h"
#include "util.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <numeric>
#include <string>
#include <unordered_map>
#include <vector>
#include <fstream>
#include <sstream>
#include <map>

#if IACA == 1
#include <iacaMarks.h>
#else
#define IACA_START 
#define IACA_END 
#endif

struct Run;

struct Input {
  private:
    static std::vector<Key> randNums(long n, long seed) {
      std::mt19937_64 rng(seed);
      std::uniform_int_distribution<Key> dist(1, (1L<<63) - 2);
      std::vector<Key> v(n);
      for (auto& x : v) x = dist(rng);
      return v;
    }

    struct KeyRandomSortedIndex {
      Key key;
      int random_index;
      int sorted_index;

      KeyRandomSortedIndex(Key k) { key = k; }
    };
    Input(size_t n) : keys(n) {
      permuted_keys.resize(n);
    }
    void fillRand(long seed) {
      std::mt19937_64 rng(seed);
      std::uniform_int_distribution<Key> dist(1, (1L<<63) - 2);
      for (size_t i = 0; i < keys.size(); i++)
        keys[i] = permuted_keys[i] = dist(rng);
      std::sort(keys.begin(), keys.end());
      sum = std::accumulate(permuted_keys.begin(), permuted_keys.end(), 0UL);
    }

  public:
    using SeedN = std::tuple<long, long>;
    static std::map<SeedN, Input> load(std::vector<Run> runs);

    std::vector<Key> permuted_keys;
    PaddedVector<> keys;
    unsigned long sum;

    // TODO use n_gets so that we don't have to store whole array permutation
  Input(std::vector<Key> nums, long n_gets) : permuted_keys(nums), keys(nums),
      sum(std::accumulate(nums.begin(), nums.end(), 0UL)) {
#ifdef DUMP_INPUT
    return;
#endif
    std::sort(keys.begin(), keys.end());
  }
};

struct Run {
  static std::vector<Run> load(std::string file_name) {
    std::ifstream f(file_name);
    std::vector<std::string> header;
    std::string line;
    std::getline(f, line);
    assert(f.good());
    for (auto [ss, word] = std::tuple(std::stringstream(line), std::string());
        ss.good(); header.push_back(word)) ss >> word;
    std::vector<Run> runs;
    for (; f.good(); ) {
      std::getline(f, line);
      std::stringstream ss(line);
      Run r;
      std::vector<int> set(4, 0);
      for (auto field : header) {
        if (!ss.good()) break;
        if (field == "n") {
          ss >> r.n;
          set[0] = 1;
        } else if (field == "seed") {
          ss >> r.seed;
          set[1] = 1;
        } else if (field == "thread") {
          ss >> r.n_thds;
          set[2] = 1;
        } else if (field == "algorithm") {
          ss >> r.name;
          set[3] = 1;
        }
      }
      if (std::all_of(set.begin(), set.end(), [](int x){ return 1==x; }))
        runs.push_back(r);
    }
    return runs;
  }

  std::string name;
  long seed, n;
  int n_thds;
  bool ok = true;

  template <class S>
  static std::vector<double> measure(Run& run, const Input& inputC) {
    // I didn't get great results trying to find a particular value for the
    // sample size, but this semed to be not terrible
    constexpr int sample_size = 1000;
    const int n_samples = inputC.keys.size() / sample_size;
    auto& queries = inputC.permuted_keys;

    S search(inputC.keys);

    //std::vector<double> ns( run.n_thds);
    std::vector<double> ns(n_samples * run.n_thds);

    // TODO break apart by threads
#pragma omp parallel default(none) num_threads(run.n_thds) shared(queries, run, search, ns)
    {
      const int tid = omp_get_thread_num();
      const auto& thread_ns = &ns[tid * n_samples];
      thread_ns[0] = 0.0;
      auto valSum = 0UL;
      for (int sample_index = 0, query_index = 0; sample_index < n_samples;
          sample_index++, query_index += sample_size) {
        if (query_index + sample_size > queries.size()) query_index = 0;

        auto t0 = std::chrono::steady_clock::now();
        for (int i = query_index; i < query_index + sample_size; i++) {
          auto val = search(queries[i]);
          valSum += val;
          assert(val == queries[i]);
        }

        auto t1 = std::chrono::steady_clock::now();
        double ns_elapsed = std::chrono::nanoseconds(t1-t0).count();
        //thread_ns[0] += ns_elapsed;
        thread_ns[sample_index] = ns_elapsed / sample_size;
        // should the sample index be included?
      }
#pragma omp critical
      run.ok = run.ok && valSum == inputC.sum;
    }
    return ns;
  }

  using Benchmark = std::vector<double>(Run&, const Input&);
  auto operator()(const Input& input) {
    static std::unordered_map<std::string, Benchmark*> fns{
      {"i", measure<InterpolationNaive>},
        {"i-precompute", measure<InterpolationPrecompute>},
        {"i-lut", measure<InterpolationRecurseLut>},
        {"i-seq-fp", measure<InterpolationLinearFp> },
        {"i-seq", measure<InterpolationLinear>},
        {"i-guard", measure<InterpolationRecurseGuard>},
        {"i-slope", measure<i_slope>},
        {"i-seq-simd", measure<i_simd>},

        {"b", measure<Binary<>>},
        {"b-cond", measure<BinaryCond>},
        {"b-noeq", measure<BinaryNoEq>},
        {"b-for", measure<BinaryFor>},
        {"b-noeq-for", measure<BinaryNoEqFor>},
        {"b-pow", measure<BinaryPow>},
        {"b-lin", measure<BinaryLin>},
    };
    auto ns = fns[this->name](*this, input);
    if (!this->ok)
      std::cerr << "mess up " << this->seed << ' ' << this->name << '\n';
    return ns;
  }
};

std::map<Input::SeedN, Input> Input::load(std::vector<Run> runs) {
  std::vector<SeedN> input_set;
  for (auto r : runs) input_set.emplace_back(r.seed, r.n);
  std::sort(input_set.begin(), input_set.end());
  input_set.resize(std::distance(input_set.begin(),
        std::unique(input_set.begin(), input_set.end())));
  std::map<SeedN, Input> inputs;
  for (auto seed_n : input_set)
    inputs.emplace(seed_n, Input(std::get<1>(seed_n)));
#pragma omp parallel for
  for (size_t i = 0; i < input_set.size(); i++)
    inputs.at(input_set[i]).fillRand(std::get<0>(input_set[i]));

  return inputs;
}


#endif
