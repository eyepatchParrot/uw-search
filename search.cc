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
  std::cout << "n,seed,n_thds,name,ns\n";
  for (auto& run : runs)
    for (auto ns : run(inputs.at(std::tuple(run.seed, run.n))))
      printf("%ld,%ld,%d,%s,%.3f\n", run.n, run.seed, run.n_thds, run.name.c_str(), ns);
}
