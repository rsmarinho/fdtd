# Rafael Marinho
# fdtd Makefile

VERSION := fdtd1d_v2.cpp
VOBJECT	:= $(subst generator,generated_source,$(VERSION:.cpp=.o))

CXX	:= g++
CXXFLAGS := -Wall -O2
INCLUDES := -lstdc++ -lm -lhdf5 -llapack -lopenblas -larmadillo -ldl
APPFLAGS := -DARMA_DONT_USE_WRAPPER -DARMA_USE_BLAS -DARMA_USE_LAPACK -DARMA_USE_HDF5

SOURCES := $(wildcard *.cpp)
OBJECTS := $(subst generator,generated_source,$(SOURCES:.cpp=.o))
APP := fdtd

all:	version

version:	$(OBJECTS)
	$(CXX) $(INCLUDES) main.o $(VOBJECT) -o $(APP)

fdtd:	$(OBJECTS)
	$(CXX) $(INCLUDES) $^ -o $(APP)

%.o:	%.cpp
	$(CXX) $(CXXFLAGS) $(APPFLAGS) $(INCLUDES) -c $^ -o $@

fd3d: fd3d.cxx
	$(CXX) $(CXXFLAGS) $(INCLUDES) $^ -o $@

clean:
	rm -rf *.o 
	rm -rf $(APP)
	rm -rf fd3d

