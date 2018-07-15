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

ssh : 
	ssh csl "cd hw/681/ && make release && ./release i-seq.tsv | ./run.py"

ssh_d : 
	ssh csl "cd hw/681/ && make debug"

clean:
	rm -f ./release ./debug ./dump

%.trace : CXXFLAGS += -DIACA -I$(HOME)/iaca/include
%.trace : release
	iaca -arch HSW -trace-cycle-count 50 -trace $@ $<

####### Build Targets #########

release : CXXFLAGS += -O3 -DNDEBUG
release: $(SOURCES) $(HEADERS)
	$(CXX) $(CXXFLAGS) $(SOURCES) -o$@ $(LDFLAGS)

debug : CXXFLAGS += -O0
debug: $(SOURCES) $(HEADERS)
	$(CXX) $(CXXFLAGS) $(SOURCES) -o$@ $(LDFLAGS)

dump : dump.cc benchmark.h
	$(CXX) $(CXXFLAGS) dump.cc -o $@ $(LDFLAGS)

