#include "benchmark.h"
#include "interpolate.h"

#include <vector>

template <int record_bytes = 8>
void dump_keys(InputBase &in, const std::string &distribution,
               const std::string &param, long n) {
  auto &keys = static_cast<Input<8> &>(in).keys;
  IBase<8>::Float<> interpolate(keys);
  int index = 0;
  for (auto key : keys) {
    std::cout << n << '\t' << distribution << '\t' << param << '\t'
              << record_bytes << '\t' << key << '\t' << index++ << '\t'
              << interpolate(key) << '\n';
    // record_bytes << '\t' << key << '\t' << index++ << '\t' << 0 << '\n';
  }
}

int main(int argc, char *argv[]) {
  std::vector<Run> runs = Run::load(std::ifstream(argv[1]));
  std::cout << "n\tdistribution\tparam\trecord\tkey\tindex\tinterpolate\n";
  for (auto &id_input : InputBase::load(runs)) {
    auto[ distribution, param, n, record_bytes ] = id_input.first;
    dump_keys<8>(*id_input.second, distribution, param, n);
    // for (auto key : static_cast<Input<8>&>(*id_input.second).keys) {
    //	std::cout << n << '\t' << distribution << '\t' << param << '\t' <<
    //record_bytes << '\t' << key << '\n';
    //}
  }
  return 0;
}
