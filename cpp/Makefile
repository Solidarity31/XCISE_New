# Build the C++ port of excise1.5.pl.
# Requires htslib (e.g. apt install libhts-dev, or conda install -c bioconda htslib)
# and a C++17 compiler. Override HTSLIB_PREFIX for a non-system htslib:
#   make HTSLIB_PREFIX=$CONDA_PREFIX

CXX      ?= g++
CXXFLAGS ?= -O3 -march=native -std=c++17 -Wall -Wextra -Wno-unused-parameter
HTSLIB_PREFIX ?=

ifneq ($(HTSLIB_PREFIX),)
  CPPFLAGS += -I$(HTSLIB_PREFIX)/include
  LDFLAGS  += -L$(HTSLIB_PREFIX)/lib -Wl,-rpath,$(HTSLIB_PREFIX)/lib
endif
LDLIBS += -lhts -lz -lpthread

all: excise15

excise15: excise15.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $< $(LDFLAGS) $(LDLIBS)

# test build: adds an INJECT=<file> hook that replaces the optimiser result,
# used only to compare post-optimiser stages against the Perl script
excise15_test: excise15.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -DEXCISE_TEST_INJECT -o $@ $< $(LDFLAGS) $(LDLIBS)

clean:
	rm -f excise15 excise15_test

.PHONY: all clean
