.SUFFIXES:
.SUFFIXES: .o .cpp
#============================================================

#Your path to Eigen
#EIGEN=/share/pkg.7/eigen/3.3.5/src/eigen-eigen-b3f3d4950030
EIGEN=/Users/hattrick/Desktop/Projects/Lattice_AdS/C-program/Eigen

TARGET	    = spectrum
C_SOURCES   = ads_graph.cpp 
C_OBJS      = ads_graph.o  
MY_INCLUDES = graph.h util.h cg.h cg_multishift.h eigen.h graph.h 

CCX = g++ -std=c++11
CXXFLAGS = -O2 -g -Wall -std=c++11 -I${EIGEN} -I. -I/share/pkg/gsl/2.3/install/include -Wall -Wno-sign-compare

#LIBDIRS = /share/pkg/gsl/2.3/install/lib/libgsl.a
LIBDIRS = -lgsl

#============================================================
all: $(TARGET)

ads_graph.o: graph.h util.h cg.h cg_multishift.h eigen.h graph.h 

ads_graph.o: $(MY_INCLUDES)

.o:.cpp	$(MY_INCLUDES)
	$(CCX) -c $(CXXFLAGS) $<

$(TARGET) : $(C_OBJS)
	$(CCX) $(CXXFLAGS)  $^ $(LIBDIRS)  -o $@

# Implicit rules: $@ = target name, $< = first prerequisite name, $^ = name of all prerequisites
#============================================================

ALL_SOURCES = Makefile $(C_SOURCES) $(MY_INCLUDES) 

NOTES =
%= otherstuff.np 

clean:
	rm -f $(TARGET) $(C_OBJS) core *.~*~

tar: $(ALL_SOURCES) $(NOTES)
	tar cvf $(TARGET).tar $(ALL_SOURCES) $(NOTES)

$(TARGET).ps: $(ALL SOURCES)
	enscript -pcode.ps $(ALL_SOURCES)
