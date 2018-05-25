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
#include <unistd.h>
#include <unordered_map>
#include <vector>
#include <x86intrin.h>
#include <map>
#include <set>

// ./x run_tsv
int main(int argc, char *argv[]) {
  using RunTuple = std::tuple<InputParam::Tuple, std::string, int>;

  std::vector<Run> runs = Run::load(std::ifstream(argv[1]));
  auto inputs = InputBase::load(runs);
  std::cout << "run\tn\tdistribution\tparam\tn_thds\talgorithm\trecord\tns\n";
  int run_ix = 0;

  RunTuple old_param;
  auto t0 = std::chrono::steady_clock::now();
  for (auto &run : runs) {
    auto[ distribution, param, n, record_bytes ] = run.input_param;
    RunTuple new_param{ run.input_param, run.name, run.n_thds };
    auto t1 = std::chrono::steady_clock::now();
    if (new_param != old_param) {
      std::cerr << '\n' << n << ' ' << distribution << ' ' << param << ' '
                << record_bytes << ' ' << run.name << ' ' << run.n_thds;
      old_param = new_param;
      t0 = t1;
    } else if (std::chrono::duration<double, std::milli>(t1 - t0).count() >
               1000.0) {
      std::cerr << '.';
      t0 = t1;
    }

    for (auto ns : run(*inputs.at(run.input_param))) {
      if (!run.ok) {
        break;
      }
      printf("%d\t%ld\t%s\t%s\t%d\t%s\t%d\t%.3f\n", run_ix, n,
             distribution.c_str(), param.c_str(), run.n_thds, run.name.c_str(),
             record_bytes, ns);
    }
    run_ix++;
  }
  std::cerr << '\n';
}
