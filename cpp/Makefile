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

SRC = main.cpp excise_core.cpp

excise15: $(SRC) excise_core.h
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $(SRC) $(LDFLAGS) $(LDLIBS)

# test build: adds an INJECT=<file> hook that replaces the optimiser result,
# used only to compare post-optimiser stages against the Perl script
excise15_test: $(SRC) excise_core.h
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -DEXCISE_TEST_INJECT -o $@ $(SRC) $(LDFLAGS) $(LDLIBS)

clean:
	rm -f excise15 excise15_test

.PHONY: all clean

# The R/Bioconductor package (bioc/EXCISE) builds the same core from its own
# src/ copy, because an R source package must be self-contained.
R_SRC = ../bioc/EXCISE/src
sync-r:
	cp excise_core.h excise_core.cpp $(R_SRC)/

check-sync:
	cmp excise_core.h $(R_SRC)/excise_core.h && cmp excise_core.cpp $(R_SRC)/excise_core.cpp && echo "R package core in sync"

.PHONY: sync-r check-sync
