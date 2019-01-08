#include<iostream>
#include "Spinodal.h" //CHANGE THIS header file per your NEED
#include "PFUtilities.h"



int main(int argc, char const *argv[])
{
  computeParam *cp = new computeParam(); //Define strucutre for computational parameters
  materialParam *mp = new materialParam(); // Define strucutre for material parameters
  defineParameters(cp,mp);

  // Define 2D concetration field
  double **Cxy = create2DField(cp);

  //Assign an initial microstrucutre, i.e here we create a random concentration fluctuation
  prepare2DRandomConcMicrostructure(Cxy,cp,mp);

  //Evolve the microstrucure using implicit Euler Finite differnece method
  evolveSpinodalMicrostructureFD(Cxy,cp,mp);

  //Evolve the microstrucure using implicit Euler FFT
  //evolveSpinodalMicrostructureFFT(Cxy,cp,mp);

  //Garbage cleaning - delete dangling pointers
  delete cp;
  delete mp;
  delete Cxy;
  return 0;
} //END function
