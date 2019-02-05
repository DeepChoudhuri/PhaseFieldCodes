Set C++ programs that solves the Cahn-Hilliard (CH) equation, i.e. time evolution of the concentration field

Specific C++ files present in the repository are:

1.Parameters.h and Parameters.cpp: stores computational parameters e.g. time step, etc.

2. Spinodal.h and Spinodal.cpp: CH equations and the potential (I have used a double well potential)

3. PFUtilities.h and PFUtilities.cpp: Utility programs that contaning generic functions meant to solve phase field problems

4. main.cpp: The driver program.

5. Makefile: Tested on Linux (Ubuntu 16.04) and OSX

5. In main.cpp we have the option of using two types of solvers

A) using finite difference, i.e. evolveSpinodalMicrostructureFD(Cxy,cp,mp) and,

B) using FFTW3 library, i.e. evolveSpinodalMicrostructureFFT(Cxy,cp,mp). (You will need to install FFTW3 in your computer)

Comment out either of the two functions to use finite difference or FFTW. In the uploaded version finite differnce is comment out.

NOTE:

1. This is work in progress, and I will not be responsible if you plan to use this code in your research. (Have fun indulging yourself with all the bugs!!)

2. Current implementation is only for 2D systems and MPI has not been implemented.

3. The code was developed using C11 standard
