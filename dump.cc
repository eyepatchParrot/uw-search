#include "benchmark.h"

#include <vector>

int main(int argc, char *argv[]) {
  std::vector<Run> runs = Run::load(argv[1]);
	std::cout << "n\tdistribution\tparam\tkey\n";
	for (auto [id, input] : Input::load(runs)) {
		auto [n, distribution, param] = id;
		for (auto key : input.keys) {
			std::cout << n << '\t' << distribution << '\t' << param << '\t' <<  key << '\n';
		}
	}
	return 0;
}
