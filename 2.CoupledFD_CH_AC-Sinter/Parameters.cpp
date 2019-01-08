#include "Parameters.h"

// Initializing the compute parameters - MODIFY these variables!!
void defineParameters(computeParam *cp, materialParam *mp)
{
  //Computational parameters
  cp->Nx    = 1000;
  cp->Ny    = 1000;
  cp->Nz    = 1;
  cp->nstep = 20000;
  cp->dx    = 0.4; // any thing below 0.4 is junk!!
  cp->dy    = 0.4;
  cp->dz    = 0.0;
  cp->dtime = 0.0001;

  /*
  //Parameters for spinodal decomposition
  mp->c0 = 0.4; //Alloy composition
  mp->M = 1.0;
  mp->K = 0.5;
  mp->noise = 0.02;
  mp->A = 1.0;
  */

  /* Parameters for coupled equation
  */
  mp->etaNum = 3;
  mp->Kc   = 5.0;
  //mp->M =  1.0; //For defining "M" in CH equation
  mp->Ke   = 2.0;
  mp->L    = 5.0;

  //Parameters for sintering. dvol, dval, dsur and dgb
  //are used for defining "M" in CH equation
  mp->dvol  = 0.04;
  mp->dvap  = 0.002;
  mp->dsur  = 16.0;
  mp->dgb   = 1.6;

  // Parameters in the Free energy expression
  mp->A = 16;
  mp->B = 1; //Works for 1,2,5

}//END function
