CXX=/usr/bin/g++

CXXFLAGS=-c -Wall -O2 -std=c++11
FFTWLIB=-lfftw3

all: sp.x

sp.x: main.o PFUtilities.o Spinodal.o Parameters.o
	$(CXX) $(FFTWLIB) main.o PFUtilities.o Spinodal.o Parameters.o -o sp.x

main.o: main.cpp
	$(CXX) $(CXXFLAGS) main.cpp

PFUtilities.o: PFUtilities.cpp PFUtilities.h
	$(CXX) $(CXXFLAGS) PFUtilities.cpp

Spinodal.o: Spinodal.cpp Spinodal.h
	$(CXX) $(CXXFLAGS) Spinodal.cpp

Parameters.o: Parameters.cpp Parameters.h
	$(CXX) $(CXXFLAGS) Parameters.cpp


run:
	./sp.x > log

clean:
	rm -f *o *.x *.vtk log* *.txt
	clear

clear:
	rm -f *o *.x *.vtk log* *.txt
	clear
