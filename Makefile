#N_RUNS=10
ifdef NSORT
	DEFINES += -DNSORT
endif
ifdef N_RUNS
	DEFINES += -DN_RUNS=$(N_RUNS)
endif
ifdef N_SAMPLES
	DEFINES += -DN_SAMPLES=$(N_SAMPLES)
endif
ifdef SUBSET_SIZE
	DEFINES += -DSUBSET_SIZE=$(SUBSET_SIZE)
endif

CXX=clang++
CXXFLAGS=-fopenmp -ffast-math -Wall -std=c++1z -fno-omit-frame-pointer -ggdb -march=native $(DEFINES) -I$(HOME)/iaca/include
LIB=-I$(HOME)/include -L$(HOME)/lib  -I$(HOME)/iaca/include
HEADERS=oracle.h interpolate.h benchmark.h bin.h lin.h util.h div.h
OBJ=
IACA=0

N_INTS=1000
SEED=42

N_THREADS=1
#BENCHMARKS=bsEq bs bsLin_32 isRecurse isLin_1 isLin_2 oracle isSub
#BENCHMARKS=isRecurse isFp isFp_slow isLin_1 isLin_1_slow bs
#BENCHMARKS=isFp isFp_slow isIDiv
#BENCHMARKS=binary-naive binary-size binary-linear interpolation-naive interpolation-recurse interpolation-linear-fp interpolation-linear oracle
BENCHMARKS=binary-naive binary-size binary-linear interpolation-naive interpolation-recurse interpolation-linear-fp interpolation-linear
BENCHMARKS=binary-linear interpolation-linear
BENCHMARKS=i-guard i-seq b-sz-lin i-simd
BENCHMARKS=fib
BENCHMARKS=i i-precompute i-seq-fp i-seq-fp-intercept i-seq i-seq-intercept
BENCHMARKS=i-seq-fp i-seq-fp-intercept i-seq-fp-pick i-seq
#BENCHMARKS=i-guard i-guard i-guard i-slope
#BENCHMARKS=b-lr b-lr-cond b-lr-over b-lr-noeq b-lr-for b-lr-noeq-for b-lr-lin
#BENCHMARKS=b-sz b-sz-cond b-sz-noeq b-sz-for b-sz-noeq-for b-sz-pow b2 b-sz-lin
RUN=./search $(N_INTS) $(SEED) $(N_THREADS) $(BENCHMARKS)

.PHONY: run search debug d_lin lin splines
run: release
	$(RUN)

iaca : IACA=1
iaca: release
	~/iaca/bin/iaca -mark 0 -arch HSW search

release : CXXFLAGS += -O3 -DNDEBUG
release : search

#release : LDFLAGS += -flto -L$(HOME)/llvm/lib

debug : CXXFLAGS += -O0
debug : search
	gdb --args $(RUN)

dump : CXXFLAGS += -DDUMP_INPUT
dump : search

# add release identifier for object files
%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

search: $(OBJ) $(HEADERS)
	$(CXX) $(CXXFLAGS) search.cc $(OBJ) -o $@ $(LDFLAGS) -DIACA=$(IACA)

div: div.cc div.h
	$(CXX) $(CXXFLAGS) $< -o $@

asm: search.cc
	$(PROFILE) -g2 -S $< -o $@.s
	$(PROFILE) $@.s -o $@

lin: lin.cc
	#$(PROFILE) -o $@ $(LIB) -fprofile-instr-generate $< -lbenchmark -lpthread 
	#./lin
	#llvm-profdata merge -output=default.profdata default.profraw
	#$(PROFILE) -o $@ $(LIB) -fprofile-instr-use $< -lbenchmark -lpthread
	$(PROFILE) -o $@ $(LIB) $< -lbenchmark -lpthread
	./lin --benchmark_out_format=csv --benchmark_out=lin.csv

d_lin: lin.cc
	$(DEBUG) -o $@ $(LIB) $< -lbenchmark -lpthread && gdb ./$@

splines: splines.cc
	$(PROFILE) $< -o $@
	#$(DEBUG) $< -o $@
	#gdb --args ./splines -f input/uniform.1000.0 -n 2

clean:
	rm -f ./search ./debug ./splines ./lin
