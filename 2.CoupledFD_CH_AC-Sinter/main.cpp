#include<iostream>
#include "CoupledEqFDSinter.h" //CHANGE THIS header file per your NEED*/
#include "PFUtilities.h"



int main(int argc, char const *argv[])
{
  /*
  Initializing microstrucural parameters
  */
  computeParam *cp = new computeParam(); //Define strucutre for computational parameters
  materialParam *mp = new materialParam(); // Define strucutre for material parameters
  defineParameters(cp,mp);

  // Define 2D concetration field
  double **Cxy = create2DField(cp);

  // Vector containing 2D slices of all eta xy fields
  std::vector<double **> Evec;
  //Create 2D slices inside the vector
  for(int i = 0; i < mp->etaNum; ++i){
    Evec.push_back(create2DField(cp));//Fill with ZERO valued matrix
  }
  std::cout << "Number of Etas: " << Evec.size() << std::endl;
  prepareInitialMicrostructure(Cxy, Evec, cp, mp);
  //writeAllDataToVTKFile(Cxy, Evec, cp, mp, 0);
  /*
  END intial microstrucutre preparation
  */

  evolveMicrostrucutre(Cxy, Evec, cp, mp);

  //Garbage cleaning - delete dangling pointers
  delete cp;
  delete mp;
  //delete2DArray(Cxy,cp->Nx);
  //Evec.clear();
  return 0;
} //END function
