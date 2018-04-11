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
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <set>
#include <unordered_map>
#include <vector>

#if IACA == 1
#include <iacaMarks.h>
#else
#define IACA_START 
#define IACA_END 
#endif

template<typename T>
T parse(std::string s) {
  T t;
  std::stringstream(s) >> t;
  return t;
}


struct Run;

struct InputParam {
  std::string distribution, param;
  long n;
  int record_bytes;

  using Tuple = std::tuple<std::string, std::string, long, int>;

  operator Tuple() { return Tuple{distribution, param, n, record_bytes}; }
};

struct InputBase {
  using InputMap = std::map<InputParam::Tuple, std::unique_ptr<InputBase>>;
  static InputMap load(std::vector<Run> runs);
};

template <int record_bytes>
struct Input : public InputBase {
  private:
    auto uniform(long seed) {
      std::mt19937_64 rng(seed);
      std::uniform_int_distribution<Key> dist(1, (1L<<63) - 2);
      std::vector<Key> v(keys.size());
      for (auto& y : v) y = dist(rng);
      std::sort(v.begin(), v.end());
      return v;
    }

    auto gap(long seed, double sparsity) {
      assert(sparsity <= 1.0);

      std::vector<Key> v(keys.size());
      long range = v.size() / sparsity;
      assert(range >= v.size());

      std::mt19937_64 rng(seed);
      // TODO make sure that this is correct
      std::uniform_int_distribution<Key> dist(1, range+1);
      std::set<Key> skips;
      while (skips.size() + v.size() < range) skips.insert(dist(rng));
      for (Key k = 1, i = 0; k < range+1; k++) if (skips.count(k) == 0) v[i++] = k;
      return v;
    }

    auto fal(double shape) {
      std::vector<Key> v(keys.size());
      auto n = v.size();
      double scale = 1.0 / (pow(n-1, -shape) - pow(n, -shape)); // (n-1)**(-s) - n**(-s)
      // [int((n-i)**(-s) / C) for i in range(n,0,-1)]
      for (auto i = 0; i < v.size(); i++)
        v[i] = pow((double)(n-i), -shape) * scale; 
      return v;
    }
    
    auto cfal(double shape) {
      auto v = fal(shape);
      std::partial_sum(v.begin(), v.end(), v.begin());
      return v;
    }

    void fill(const std::vector<Key>&& v, long seed = 42) {
      keys = std::move(v);
      std::copy(keys.begin(), keys.end(), permuted_keys.begin());
      std::shuffle(permuted_keys.begin(), permuted_keys.end(), std::mt19937(seed));
    }

  public:
    using Id = std::tuple<long, std::string, std::string>;

    std::vector<Key> permuted_keys;
    PaddedVector<record_bytes> keys;
    unsigned long sum;
    
  Input(const long n, const std::string& distribution, const std::vector<std::string>& params) : permuted_keys(n), keys(n) {
    auto param = params.begin();
    // uniform - seed
    // gaps    - seed, range
    // FB      - file
    // fal     - shape
    // cfal    - shape
    if (distribution == "uniform") {
      auto seed = parse<long>(param[0]);
      fill(uniform(seed));
    } else if (distribution == "gap") {
      auto [seed, sparsity] = std::tuple{parse<long>(param[0]),
        parse<double>(param[1])};
      fill(gap(seed, sparsity));
    } else if (distribution == "fal") {
      auto shape = parse<double>(param[0]);
      fill(fal(shape));
    } else if (distribution == "cfal") {
      auto shape = parse<double>(param[0]);
      fill(cfal(shape));
    } else {
      assert(!"No distribution found.");
    }

    sum = std::accumulate(keys.begin(), keys.end(), 0UL);
  }

};

struct Run {
  static auto reverse_index(std::vector<std::string> v) {
    std::map<std::string, long> m;
    int i = 0;
    for (auto s : v) m[s] = i++;
    return m;
  }

  static auto load(std::ifstream&& file) {
    auto header = reverse_index(split(read_line(file), '\t'));
    assert(header.size() == 6);
    std::vector<Run> runs;
    for (; file.good(); ) {
      auto fields = split(read_line(file), '\t');
      if (!file.good()) break;
      assert(fields.size() == header.size());
      InputParam input_param{
          .distribution = fields[header["distribution"]],
          .param = fields[header["param"]],
          .n = parse<long>(fields[header["n"]]),
          .record_bytes = parse<int>(fields[header["payload"]])
      };
      runs.emplace_back(input_param, fields[header["algorithm"]],
          parse<int>(fields[header["thread"]]));
    }
    return runs;
  }

  InputParam input_param;
  std::string name;
  int n_thds;
  // TODO consider removing
  bool ok;

  Run(InputParam input_param, std::string name, int n_thds) : input_param(input_param), name(name), n_thds(n_thds), ok(true) {}

  template <typename Search, int record_bytes>
  static std::vector<double> measure2(Run& run, const InputBase& input) {
    const auto& inputC = static_cast<const Input<record_bytes>&>(input);
#ifdef INFINITE_REPEAT
    constexpr bool infinite_repeat = true;
#else
    constexpr bool infinite_repeat = false;
#endif
    // I didn't get great results trying to find a particular value for the
    // sample size, but this seemed to be not terrible
    constexpr int sample_size = 1000;
    const int n_samples = inputC.keys.size() / sample_size;
    auto& queries = inputC.permuted_keys;

    // TODO this can't be a template of a template have to specialize earlier
    // have to specialize in the class itself. Maybe template macros?
    Search search(inputC.keys);

    //std::vector<double> ns( run.n_thds);
    std::vector<double> ns(n_samples * run.n_thds);
    // TODO allow this to run indefinitely if an appropriate flag is set. Ensure memory is O(1)

    // TODO break apart by threads
#pragma omp parallel default(none) num_threads(run.n_thds) shared(queries, run, search, ns)
    {
      const int tid = omp_get_thread_num();
      const auto& thread_ns = &ns[tid * n_samples];
      thread_ns[0] = 0.0;
      auto valSum = 0UL;
      for (int sample_index = 0, query_index = 0;;
          sample_index++, query_index += sample_size) {
        if (query_index + sample_size > queries.size()) {
          if (sample_index == n_samples) {
            if (!infinite_repeat || valSum != inputC.sum) break;
            valSum = sample_index = 0;
          }
          query_index = 0;
        }

        auto t0 = std::chrono::steady_clock::now();
        for (int i = query_index; i < query_index + sample_size; i++) {
          auto val = search(queries[i]);
          valSum += val;
          assert(val == queries[i]);
        }

        auto t1 = std::chrono::steady_clock::now();
        double ns_elapsed = std::chrono::nanoseconds(t1-t0).count();
        //thread_ns[0] += ns_elapsed;
        if (!infinite_repeat) thread_ns[sample_index] = ns_elapsed / sample_size;
        // should the sample index be included?
      }
#pragma omp critical
      run.ok = run.ok && valSum == inputC.sum;
    }
    //std::cerr << search.err << '\n';
    return ns;
  }

  template <int record_bytes>
    static auto measure(Run& run, const InputBase& input) {
      using Fn = std::vector<double>(Run&, const InputBase&);
      static std::unordered_map<std::string, Fn*> fns{
        {"i-naive", measure2<i_naive(record_bytes), record_bytes>},
          {"i-opt", measure2<i_opt(record_bytes), record_bytes>},
          {"i-seq", measure2<i_seq(record_bytes), record_bytes>},
          {"i-recompute", measure2<i_recompute(record_bytes), record_bytes>},
          {"i-no-guard", measure2<i_no_guard(record_bytes), record_bytes>},
          {"i-fp", measure2<i_fp(record_bytes), record_bytes>},
          {"i-idiv", measure2<i_idiv(record_bytes), record_bytes>},
          {"i-hyp", measure2<i_hyp(record_bytes), record_bytes>},
          {"b-lin", measure2<b_lin(record_bytes), record_bytes>},
      };
      auto ns = fns[run.name](run, input);
      return ns;
    }

//  auto measure(const InputBase& input) {
//    return fns[this->name](input);
//  }
//
  auto operator()(const InputBase& input) {
    auto [n, distribution, param, record_bytes] = input_param;
    std::cerr << "run " << n << ' ' << distribution << ' ' << param << ' ' << record_bytes << ' ' << name << '\n';
    std::vector<double> ns;
    switch (input_param.record_bytes) {
#define CASE(N) \
      case N: ns = measure<N>(*this, input); break
      CASE(8);
      CASE(32);
      CASE(128);
#undef CASE

      default: assert(!"payload not supported");
    }
    if (!this->ok)
      std::cerr << "mess up " << param << ' ' << this->name << '\n';
    return ns;
  }
};

InputBase::InputMap InputBase::load(std::vector<Run> runs) {
  InputMap inputs;
  for (auto r : runs) {
    auto input_param = r.input_param;
    std::cerr << "load " << input_param.n << ' ' << input_param.distribution << ' ' << input_param.param << '\n';
    if (inputs.count(input_param) == 0) {
      auto distribution_param = split(input_param.param);
      inputs.emplace((InputParam::Tuple)input_param, [=]() {
          switch (input_param.record_bytes) {
#define CASE(N) \
          case N: return static_cast<std::unique_ptr<InputBase>>( \
                      std::make_unique<Input<N>>(input_param.n, \
                        input_param.distribution, distribution_param))
          CASE(8);
          CASE(32);
          CASE(128);
#undef CASE
          default: assert(!"payload size not supported");
          };
          return std::make_unique<InputBase>();}());
    }
  }
  return inputs;
}


#endif
