SOURCES=runtest.cc
DEPS=thermalEvent.hh pythiaEvent.hh csSubtractor.hh softDropGroomer.hh treeWriter.hh jetMatcher.hh
EXECUTABLE=runtest
OBJECTS=$(SOURCES:.cpp=.o)
FASTJET=/afs/cern.ch/user/m/mverweij/work/soft/toy/fastjet-install
PYTHIA8LOCATION=/afs/cern.ch/sw/lcg/external/MCGenerators_lcgcmt67c/pythia8/226/x86_64-slc6-gcc48-opt
LIBDIRARCH=lib

CC=g++
CFLAGS=-c -I$(ROOTSYS)/include #-I$(PYTHIA8LOCATION)/include 

LDFLAGS=-L$(ROOTSYS)/lib -L$(PYTHIA8LOCATION)/lib -lpythia8 $(shell root-config --glibs) $(shell root-config --cflags) -I$(PYTHIA8LOCATION)/include `$(FASTJET)/bin/fastjet-config --cxxflags --libs --plugins` -lfastjetcontribfragile -L$(PYTHIA8LOCATION)/lib -lpythia8

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) $(DEPS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o: $(PYTHIA8LOCATION)/libpythia8.a
	$(CC) $(CFLAGS) $< -o $@

.PHONY: clean

clean:
	rm -rf *.o *~ $(EXECUTABLE)

