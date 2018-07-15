CXX=$(HOME)/clang5/bin/clang++
CXXFLAGS=-std=c++17 -O3 -march=native
LDFLAGS=-lbenchmark -lpthread

ssh: sequential-bench.cc
	rsync -r . cl:~/uw-search/bench
	ssh cl "cd uw-search/bench && make sequential-bench"

sequential-bench: sequential-bench.cc
	$(CXX) $(CXXFLAGS) $< $(LDFLAGS) -o $@
