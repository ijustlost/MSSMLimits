ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)

CXX    = g++
CXXFLAGS      = -g -Wall $(ROOTCFLAGS) -std=c++11 

EXTRALIBS      = -lMathMore
LIBS    = -lm $(ROOTLIBS) ${ROOTGLIBS} $(EXTRALIBS)

default: all

SOURCES=limitCalc.cxx LimitUtilities.cxx limits_ZZ.cxx BlobUtils.cxx
OBJECTS=$(SOURCES:.cxx=.o)
EXECUTABLE=limits_ZZ

all: $(OBJECTS) $(EXECUTABLE)
	    
$(EXECUTABLE): $(OBJECTS) 
	$(CXX) $(LDFLAGS) $(OBJECTS) -o $@ $(ROOTLIBS)

%.o: %.cxx
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(ROOTCFLAGS)

clean:
	rm -rf *.o limits_ZZ
