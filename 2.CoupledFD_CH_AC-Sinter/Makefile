CXX=/usr/bin/g++

CXXFLAGS=-c -Wall -O2 -std=c++11
FFTWLIB=-lfftw3

all: ceq.x

ceq.x: main.o PFUtilities.o CoupledEqFDSinter.o Parameters.o
	$(CXX) $(FFTWLIB) -O2 main.o PFUtilities.o CoupledEqFDSinter.o Parameters.o -o ceq.x

main.o: main.cpp
	$(CXX) $(CXXFLAGS) main.cpp

PFUtilities.o: PFUtilities.cpp PFUtilities.h
	$(CXX) $(CXXFLAGS) PFUtilities.cpp

CoupledEqSinter.o: CoupledEqFDSinter.cpp CoupledEqFDSinter.h
	$(CXX) $(CXXFLAGS) CoupledEqFDSinter.cpp

Parameters.o: Parameters.cpp Parameters.h
	$(CXX) $(CXXFLAGS) Parameters.cpp


run:
	./ceq.x > log

clean:
	rm -f *o *.x *.vtk log* *.txt
	clear

clear:
	rm -f *o *.x *.vtk log* *.txt
	clear
