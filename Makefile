CHOME=$(HOME)/clang5
CXX=$(CHOME)/bin/clang++
CXXFLAGS=-std=c++17 -stdlib=libc++
LDFLAGS=-fopenmp

HEADERS=benchmark.h bin.h padded_vector.h interpolate.h util.h div.h lin.h
SOURCES=search.cc

.PHONY: run gdb clean perf

##### Run Targets ######
run : release
run :
	./release i-seq.tsv

gdb : debug
gdb :
	gdb --args ./debug i-seq.tsv 

perf : CXXFLAGS += -O3 -DNDEBUG -DINFINITE_REPEAT
perf :
	$(CXX) $(CXXFLAGS) $(SOURCES) -o$@ $(LDFLAGS)
	perf record -F99 -g ./perf i-seq.tsv

clean:
	rm -f ./release ./debug

####### Build Targets #########

release : CXXFLAGS += -O3 -DNDEBUG
release:
	$(CXX) $(CXXFLAGS) $(SOURCES) -o$@ $(LDFLAGS)

debug : CXXFLAGS += -O0
debug:
	$(CXX) $(CXXFLAGS) $(SOURCES) -o$@ $(LDFLAGS)
