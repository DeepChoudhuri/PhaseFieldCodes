#include "Parameters.h"

// Initializng the compute parameters - MODIFY these variables!!
void defineParameters(computeParam *cp, materialParam *mp)
{
  cp->Nx = 64;
  cp->Ny = 64;
  cp->Nz = 1;
  cp->nstep = 10000;
  cp->dx = 1.0;
  cp->dy = 1.0;
  cp->dz = 0.0;
  cp->dtime = 0.02;

  mp->c0 = 0.4; //Alloy composition
  mp->M = 1.0;
  mp->K = 0.5;
  mp->noise = 0.02;
  mp->A = 1.0;
}//END function
