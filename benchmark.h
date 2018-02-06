#ifndef BENCHMARK_H
#define BENCHMARK_H

#include "oracle.h"
#include "interpolate.h"
#include "bin.h"
#include "lin.h"
//#include "fib.h"

#include "omp.h"
#include <string>
#include <vector>
#include "util.h"
#include <numeric>
#include <iostream>
#include <unordered_map>
#include <cmath>
#include <vector>
#include <algorithm>
#include <chrono>

#if IACA == 1
#include <iacaMarks.h>
#else
#define IACA_START 
#define IACA_END 
#endif

struct TestStats {
  std::string name;
  std::vector<double> ns;
  bool ok;
};

using Benchmark = TestStats(const std::vector<Key>&, const std::vector<int>&, const int);

double mode(std::vector<double> X) {
  std::unordered_map<int, int> cnt;
  //for (auto x : X) cnt[round(x*100.0)]++;
  for (auto x : X) {
    int ix = round(x*100.0);
    if (0 == cnt.count(ix)) cnt[ix] = 0;
    cnt[ix]++;
  }
  auto it = std::max_element(cnt.begin(), cnt.end(),
      [](std::pair<int, int> l, std::pair<int, int> r) {
      return l.second < r.second;
      });
  std::cerr << it->second << '\n';
  return it->first / 100.0;
}

template <class S>
TestStats benchmark(
    const std::vector<Key>& input,
    const std::vector<int>& indexes,
    const int nThreads = 1) {
  S search(input, indexes);
  constexpr int nRuns = N_RUNS;
  TestStats ts{.ok=true, .ns = std::vector<double>(nRuns*nThreads)};
  std::vector<int> t_ok(nThreads);

  const std::vector<Key> vals = std::move([&]() {
    std::vector<Key> vals;
    vals.reserve(indexes.size());
    for (int i : indexes) vals.push_back(input[i]);
    return vals; }());
  
  // get verification info
  const auto expSum = [&vals](){
    auto expSum = 0UL;
    for (auto j=0;j<nRuns;j++) for (auto v : vals) expSum += v;
    return expSum; }();
  auto runSum = 0UL;
  for (auto v : vals) runSum += v;
#pragma omp parallel num_threads(nThreads) firstprivate(vals) shared(ts)
  {
    const int tid = omp_get_thread_num();
    auto valSum = 0UL;
    for (int runIx = 0; runIx < nRuns; runIx++) {
      auto t0 = std::chrono::steady_clock::now();

      for (int i = 0; i < vals.size(); i++) {
#ifndef NDEBUG
        std::cout << "answer " << indexes[i] << '\n';
#endif
        auto val = search(vals[i]);
#ifndef NDEBUG
        std::cout << '\n';
#endif
        valSum += val;
        assert(val == vals[i]);
      }

      auto t1 = std::chrono::steady_clock::now();
      ts.ns[(tid*nRuns)+runIx] = (std::chrono::nanoseconds(t1-t0)).count() / (double)vals.size();
    }
    t_ok[tid] = valSum == expSum;
  }
  for (auto ok : t_ok) ts.ok = ts.ok && !!(ok);
  //std::cerr << "Mode " << mode(ts.ns) << '\n';
  return ts;
}

TestStats benchmark(
    const std::vector<Key>& keys,
    const std::vector<int>& indexes,
    const int nThreads,
    const std::string& name,
    long seed) {
  static std::unordered_map<std::string, Benchmark*> fns{
    {"b0", benchmark<B0>},
    {"b1", benchmark<B1>},
    {"i", benchmark<InterpolationNaive>},
    {"i-precompute", benchmark<InterpolationPrecompute>},
    {"i-lut", benchmark<InterpolationRecurseLut>},
    {"i-seq-fp", benchmark<InterpolationLinearFp> },
    {"i-seq-fp-intercept", benchmark<i_seq_fp_intercept> },
    {"i-seq-fp-pick", benchmark<i_seq_fp_pick> },
    {"i-seq", benchmark<InterpolationLinear>},
    {"i-seq-intercept", benchmark<i_seq_intercept> },
    {"i-guard", benchmark<InterpolationRecurseGuard>},
    {"i-slope", benchmark<i_slope>},
    {"i-slope-lut", benchmark<i_slope_lut>},
    {"i-seq-simd", benchmark<i_simd>},

    {"b-lr", benchmark<BinaryLR<>>},
    {"b-lr-cond", benchmark<BinaryLRCond>},
    {"b-lr-noeq", benchmark<BinaryLRNoEq>},
    {"b-lr-for", benchmark<BinaryLRFor>},
    {"b-lr-noeq-for", benchmark<BinaryLRNoEqFor>},
    {"b-lr-over", benchmark<BinaryLROver>},
    {"b-lr-lin", benchmark<BinaryLRLin>},

    {"b-sz", benchmark<BinarySz<>>},
    {"b-sz-cond", benchmark<BinarySzCond>},
    {"b-sz-noeq", benchmark<BinarySzNoEq>},
    {"b-sz-for", benchmark<BinarySzFor>},
    {"b-sz-noeq-for", benchmark<BinarySzNoEqFor>},
    {"b-sz-pow", benchmark<BinarySzPow>},
    {"b-sz-lin", benchmark<BinarySzLin>},
    {"b2", benchmark<B2>},
    //{"fib", benchmark<Fib>},
    //{"set", benchmark<InterpolationSet>},

    {"oracle", benchmark<Oracle> }
  };
#ifndef NDEBUG
  std::cout << "BENCHMARK " << name << '\n';
#endif
  auto ts = fns[name](keys, indexes, nThreads);
  if (!ts.ok)
    std::cerr << "mess up " << seed << ' ' << name << '\n';
  ts.name = name;
  return ts;
}

#endif
