.SUFFIXES:
.SUFFIXES: .o .cpp
#============================================================

#Your path to Eigen
#EIGEN=/share/pkg.7/eigen/3.3.5/src/eigen-eigen-b3f3d4950030
EIGEN=/Users/hattrick/Eigen

SPECTRUM_INCLUDES = graph.h util.h cg.h cg_multishift.h eigen.h graph.h

CXX = g++ -std=c++11
CXXFLAGS = -O2 -g -Wall -std=c++11 -I${EIGEN} -I. -I/share/pkg/gsl/2.3/install/include -Wall -Wno-sign-compare

LIBDIRS = -lgsl

#============================================================
all: spectrum

spectrum: ads_graph.o
	$(CXX) $(CXXFLAGS) -o spectrum ads_graph.o $(LIBDIRS) 

ads_graph.o: ads_graph.cpp Makefile ${SPECTRUM_INCLUDES}
	$(CXX) $(CXXFLAGS) -c ads_graph.cpp
#============================================================

#============================================================
all: deltalin

deltalin: delta_lin_fit.o
	$(CXX) $(CXXFLAGS) -o deltalin delta_lin_fit.o $(LIBDIRS)

delta_lin_fit.o: delta_lin_fit.cpp Makefile 
	$(CXX) $(CXXFLAGS) -c delta_lin_fit.cpp
#============================================================

ALL_SOURCES = Makefile ads_graph.cpp delta_lin_fit.cpp $(SPECTRUM_INCLUDES)

clean:
	rm -f spectrum deltalin *.o core*
