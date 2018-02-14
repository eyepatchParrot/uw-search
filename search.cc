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
  std::vector<Run> runs = Run::load(argv[1]);
  auto inputs = Input::load(runs);
  std::cout << "run,n,seed,n_thds,name,ns\n";
  int run_ix = 0;
  for (auto& run : runs) {
    //double sum = 0.;
    //for (auto ns : run(inputs.at(std::tuple(run.seed, run.n))))
    //  sum += ns;
    //std::cout << sum << '\n';
    for (auto ns : run(inputs.at(std::tuple(run.seed, run.n))))
      printf("%d,%ld,%ld,%d,%s,%.3f\n", run_ix, run.n, run.seed, run.n_thds, run.name.c_str(), ns);
    run_ix++;
  }
}
