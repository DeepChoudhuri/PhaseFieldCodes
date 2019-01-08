Purpose: Set of C++ programs that solves the coupled Cahn-Hilliard (CH) and Allen-Cahn equations.

Specific C++/python/other files present in the repository are:

1. Parameters.h and Parameters.cpp: stores computational parameters e.g. time step, etc.

2. CoupledEqFDSinter.h and CoupledEqFDSinter.cpp: Collections of functions for creating a microstructure and solving the coupled CH and AC eqs

3. PFUtilities.h and PFUtilities.cpp: Utility programs that containing generic functions meant to solve phase field problems

4. main.cpp: The driver program.

5. Makefile: Tested on OSX

6. plot.py: For plotting Free energies as a function of time and time steps

7. Result folder: Contains results from two test cases

NOTE:

1. This is work in progress, and I will not be responsible if you plan to use this code in your research. (Have fun indulging yourself with all the bugs!!)

2. Current implementation is only for 2D systems and MPI has not been implemented.

3. The code was developed using the C11 standard

4. View the otput *.vtk files in the "Results" folder using paraview

5. Contact me for more details : deep.choudhuri@gmail.com
