#include "Spinodal.h"
#include "PFUtilities.h"
#include<cmath>




/*Prepare the intial microstrucutre by populating the concentration field
  with random numbers between 0 and 1*/
void prepare2DRandomConcMicrostructure(double **Axy, \
                                       const computeParam *cp, \
                                       const materialParam *mp)
{
  srand((unsigned)time( NULL));

  for(int i = 0; i < cp->Nx; ++i){
    for(int j = 0; j < cp->Ny; ++j){
      Axy[i][j] = mp->c0 + mp->noise*(0.5 - (float)rand()/RAND_MAX);
    }
  }
}//END function


//Bulk Free and gradient energies for spinodal decomposition - for a 2D system
double FBulkSpinodal(double **Cxy, const computeParam *cp, const materialParam *mp)
{
  double Fbulk = 0.0;
  int ip, jp;

  for(int i = 0; i < cp->Nx-1; ++i){
      //ip = i + 1;
    for(int j = 0; j < cp->Ny-1; ++j){
      ip = i + 1;
      jp = j + 1;
      Fbulk = Fbulk + mp->A*(std::pow(Cxy[i][j],2.0)*std::pow((1.0 - Cxy[i][j]),2.0)) + \
              0.5*mp->K*(std::pow((Cxy[ip][j] - Cxy[i][j]),2.0) + \
                          std::pow((Cxy[i][jp] - Cxy[i][j]),2.0));
    }
  }
  return Fbulk;
}//END function



//Local gradient of free energy
double DFlocSpinodal(double Cxy, const materialParam *mp)
{
  double dfCon = 2.0*mp->A*(Cxy*std::pow(1.0-Cxy,2.0) - (1.0-Cxy)*std::pow(Cxy,2.0));
  return dfCon;
}//END function



//Evolve the microstrucutre using Finite differnce method
void evolveSpinodalMicrostructureFD(double **Cxy, const computeParam *cp, const materialParam *mp)
{
  //prepare2DRandomConcMicrostrucutre(Cxy,cp,mp);
  //write2DVTKfile(Cxy,cp,0);
  double *FBulk = new double[cp->nstep];
  double *time = new double[cp->nstep];

  double **lapConc;
  double **lapIdFdc;

  // Interior variational derivative of bulk Free Energy wrt concentration, i.e. dF/dc
  double **idFdc = create2DField(cp);

  double tt = 0;
  //START main time loop
  for(int t = 0; t < cp->nstep; t++){
    time[t] = tt;

    // Laplacian of the concentration field
    lapConc = get2DLaplacian(Cxy,cp);

    //START loop for evaluating the interior varational dF/dc derivative
    for(int i = 0; i < cp->Nx; i++){
      for(int j = 0; j < cp->Ny; j++){
        idFdc[i][j] = DFlocSpinodal(Cxy[i][j],mp) - mp->K*lapConc[i][j];
      }
     }//END loop for evaluating the interior varational dF/dc derivative


    /*Get Laplacian of the inner variational derivative
    in the Cahn-Hillard formulation. To be used
    in the final evaluation of concentration field in the following*/
    lapIdFdc = get2DLaplacian(idFdc,cp);

    //START loop for evaluating gradient Grad(dF/dc)
    for(int i = 0; i < cp->Nx; ++i){
      for(int j = 0; j < cp->Ny; ++j){

        Cxy[i][j] = Cxy[i][j] + mp->M*cp->dtime*lapIdFdc[i][j];

        if(Cxy[i][j] >= 0.9999){ //Set MAX limits for small variations
          Cxy[i][j] = 0.9999;
        }

        if(Cxy[i][j] < 0.00001){ //Set MIN limits for small variations
          Cxy[i][j] = 0.00001;
        }
      }
     }//END loop for evaluating gradient Grad(dF/dc)*/

    if((t%100 == 0) && (t<1001)){
      write2DVTKfile(Cxy,cp,t);
    }
    FBulk[t] = FBulkSpinodal(Cxy, cp, mp);
    tt = tt + cp->dtime;
  }//END main time loop

  write1DTXTfile(time, FBulk, cp->nstep, "energy-time.txt");
  write2DVTKfile(Cxy,cp,cp->nstep);

  //Garbage cleaning - delete dangling pointers
  delete []FBulk; //1D array
  delete []time;  //1D array
  delete2DArray(lapConc, cp->Nx);
  delete2DArray(lapIdFdc, cp->Nx);
  delete2DArray(idFdc, cp->Nx);

} //END function



//Evolve the microstrucutre using FFT method
void evolveSpinodalMicrostructureFFT(double **Cxy, const computeParam *cp, const materialParam *mp)
{
  //displayArray(Cxy,cp->Nx, cp->Ny);
  double **dFdc = create2DField(cp);

  double *K2 = new double[cp->Nx * cp->Nx];
  double *K4 = new double[cp->Nx * cp->Nx];
  prepare2DFFTSquaredBasis(K2, K4, cp);

  fftw_complex *cxyK;
  fftw_complex *dfdcK;

  double *FBulk = new double[cp->nstep];
  double *time = new double[cp->nstep];

  double tt = 0;
  for(int t = 0; t < cp->nstep; t++){
    time[t] = tt;

    /*Fill up 2D array for the gradiaent of Free energy wrt to Concentration
    i.e dF/dc*/
    for(int i = 0; i < cp->Nx; ++i){
      for(int j = 0; j < cp->Ny; ++j){
        dFdc[i][j] = DFlocSpinodal(Cxy[i][j],mp);
      }
    }
    //displayArray(dFdc,cp->Nx, cp->Ny);


    //Get FFT of dF/dc as 1D arrays
    dfdcK = get2DFFT(dFdc, cp->Nx, cp->Ny);
    //display2DFFTdata(dfdcK, cp->Nx, cp->Ny);
    //Get FFT of the concentration field as 1D arrays
    cxyK = get2DFFT(Cxy, cp->Nx, cp->Ny);
    //display2DFFTdata(cxyK, cp->Nx, cp->Ny);

    /*Do time integration*/
    int k = 0;

    /*Numerators*/
    double numeREAL = 0.0;
    double numeCMPLX = 0.0; //Just the values of the complex PORTION

    double denom = 0.0; // Denomenator

    double factor = mp->M*cp->dtime; //Multiplying factor

    for(int i = 0; i < cp->Nx; ++i){
      for(int j = 0; j < cp->Ny; ++j){
        denom = 1 + factor*mp->K*K4[k];
        numeREAL = cxyK[k][REAL] - factor*dfdcK[k][REAL]*K2[k];
        numeCMPLX = cxyK[k][IMAG] - factor*dfdcK[k][IMAG]*K2[k];
        cxyK[k][REAL] = numeREAL/denom;
        cxyK[k][IMAG] = numeCMPLX/denom;
        k++;
      }
    }

    //display2DFFTdata(cxyK, cp->Nx, cp->Ny);
    //Get Real values after IFFT
    Cxy = get2DREAL(get2DiFFT(cxyK,cp->Nx,cp->Ny), cp->Nx, cp->Ny);
    //displayArray(Cxy,cp->Nx, cp->Ny);

    //Set limits
    for(int i = 0; i < cp->Nx; ++i){
      for(int j = 0; j < cp->Ny; ++j){

        if(Cxy[i][j] >= 0.9999){ //Set MAX limits for small variations
          Cxy[i][j] = 0.9999;
        }

        if(Cxy[i][j] < 0.00001){ //Set MIN limits for small variations
          Cxy[i][j] = 0.00001;
        }
      }
     }//END loop for evaluating gradient Grad(dF/dc)*/

    if((t%100 == 0) && (t<1001)){
       write2DVTKfile(Cxy,cp,t);
    }
    FBulk[t] = FBulkSpinodal(Cxy, cp, mp); //Compute bulk energy
    tt = tt + cp->dtime;
  }

  write1DTXTfile(time, FBulk, cp->nstep, "energy-time-FFT.txt");
  write2DVTKfile(Cxy,cp,cp->nstep);

  //Garbage collection

  delete2DArray(dFdc, cp->Nx);

  delete []K2;
  delete []K4;
  delete []FBulk;
  delete []time;

  fftw_free(cxyK);
  fftw_free(dfdcK);

}//END function

















































//
