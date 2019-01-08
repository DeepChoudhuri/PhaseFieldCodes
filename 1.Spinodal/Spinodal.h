#ifndef GUARD_Spinodal_h
#define GUARD_Spinodal_h

#include "Parameters.h"

//Funtions related to microstrucutre creation
void prepare2DRandomConcMicrostructure(double **, \
                                       const computeParam *, \
                                       const materialParam *);
                                       
//Bulk Free energy for spinodal decomposition for a 2D system
//MODIFY this for a differnt system
double FBulkSpinodal(double **, const computeParam *, const materialParam *);

//Local gradient of free energy
//MODIFY this for a differnt system
double DFlocSpinodal(double , const materialParam *);

//Evolve the microstrucutre using Finite differnce method
void evolveSpinodalMicrostructureFD(double **, const computeParam *, const materialParam *);

//Evolve the microstrucutre using FFT
void evolveSpinodalMicrostructureFFT(double **, const computeParam *, const materialParam *);


#endif
