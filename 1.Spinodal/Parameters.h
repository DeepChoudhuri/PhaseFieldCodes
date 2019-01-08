#ifndef GUARD_Parameters_h
#define GUARD_Parameters_h

/* Parameters for numerical solutions*/
struct computeParam
{
  long Nx, Ny, Nz; //defines the 2D grid
  long nstep;
  double dtime;
  double dx, dy, dz; //step sizes
};

/*Material parameters*/
struct materialParam
{
  double c0; //Alloy composition
  double M; //Mobility
  double K; // Coeeficient of gradient
  double noise; // noise term
  double A; //Constant in free energy formulation
};

// Define parameters
void defineParameters(computeParam *, materialParam *);

#endif
