CXX      := g++
LD       := g++
BOOSTDIR := /cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/boost/1.63.0-gnimlf
PYDIR := /cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/python/2.7.14-omkpbe4
CXXFLAGS := -Wall -Wextra -O2 -std=c++1y -g -I../../ $(shell root-config --cflags) -I$(BOOSTDIR)/include -I$(PYDIR)/include/python2.7
LDFLAGS  := -shared -Wall -Wextra $(shell root-config --libs) -L $(BOOSTDIR)/lib -lboost_python -L$(PYDIR)/lib -lpython2.7
SRC 	 := $(wildcard src/*.cc)
BINSRC 	 := $(wildcard test/*.cpp)
OBJS     := $(SRC:src/%.cc=obj/lib/%.o)
BINOBJS  := $(BINSRC:test/%.cpp=obj/bin/%.o)
BINS 	 := $(BINSRC:test/%.cpp=bin/%)
LIB 	 := lib/libAcornNTupler.so

-include $(OBJS:.o=.d)
-include $(BINOBJS:.o=.d)

obj/lib/%.o: src/%.cc
	$(CXX) $(CXXFLAGS) -fPIC -c $< -o $@
	$(CXX) $(CXXFLAGS) -MM -MP -MT $@ $< -o $(@:.o=.d)

obj/bin/%.o: test/%.cpp
	$(CXX) $(CXXFLAGS) -fPIC -c $< -o $@
	$(CXX) $(CXXFLAGS) -MM -MP -MT $@ $< -o $(@:.o=.d)


lib/AcornNTupler_xr.cc:
	genreflex src/classes.h -s src/classes_def.xml -o lib/AcornNTupler_xr.cc --deep --fail_on_warnings \
	--rootmap=lib/AcornNTupler_xr.rootmap --rootmap-lib=lib/libAcornNTupler.so --capabilities=lib/AcornNTupler_xi.cc \
	-DCMS_DICT_IMPL -D_REENTRANT -DGNUSOURCE -D__STRICT_ANSI__ -DGNU_GCC -D_GNU_SOURCE \
	-I../../ -DCMSSW_REFLEX_DICT

lib/AcornNTupler_xr.cc.o: lib/AcornNTupler_xr.cc
	c++ -MMD -MF lib/AcornNTupler_xr.cc.d -c -I../../ -I./ -DGNU_GCC -D_GNU_SOURCE $(shell root-config --cflags)  \
	-DCMSSW_REFLEX_DICT -pthread -Os -Wno-unused-variable -g -fPIC lib/AcornNTupler_xr.cc -o lib/AcornNTupler_xr.cc.o


$(LIB): $(OBJS) lib/AcornNTupler_xr.cc.o
	$(LD) $(LDFLAGS) -o $@ $^

bin/%: obj/bin/%.o $(LIB)
	$(CXX) -o $@ $^ $(shell root-config --libs)

all: $(LIB) $(BINS)

.SECONDARY: $(BINOBJS) obj/lib/AcornNTupler_xr.cc
.PHONY: all clean
.DEFAULT_GOAL := all

clean:
	rm obj/lib/* obj/bin/* lib/*
