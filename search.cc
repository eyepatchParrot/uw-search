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
  std::cout << "run\tn\tseed\tn_thds\tname\tns\n";
  int run_ix = 0;
  for (auto& run : runs) {
    //double sum = 0.;
    //for (auto ns : run(inputs.at(std::tuple(run.seed, run.n))))
    //  sum += ns;
    //std::cout << sum << '\n';
    for (auto ns : run(inputs.at(std::tuple(run.n, run.param))))
      printf("%d\t%ld\t%s\t%d\t%s\t%.3f\n", run_ix, run.n, run.param.c_str(), run.n_thds, run.name.c_str(), ns);
    run_ix++;
  }
}
