#include "PFUtilities.h"

/* Get Lapnacian using the 5-point stencil method.
Also check for periodic boundary conditions*/
double **get2DLaplacian(double **Axy, const computeParam *cp)
{
  int im, ip, jm, jp;

  double **dummy = create2DField(cp);
  /*dummy = new double*[cp->Nx];
  for(int i = 0; i < cp->Nx; i++){
    dummy[i] = new double[cp->Ny];
  }*/

  for(int i = 0; i< cp->Nx; i++){
    for(int j = 0; j< cp->Ny; j++){
      jp = j + 1;
      jm = j - 1;
      ip = i + 1;
      im = i - 1;

      // Check for periodic boundary conditions
      if(im < 0){
        im = cp->Nx - 1;
      }
      if(ip == cp->Nx){
        ip = 0;
      }

      if(jm < 0){
        jm = cp->Ny - 1;
      }
      if(jp == cp->Ny){
        jp = 0;
      }
      //Creating five point stencil points
      double aE = Axy[ip][j];
      double aW = Axy[im][j];
      double aN = Axy[i][jp];
      double aS = Axy[i][jm];
      double aC = Axy[i][j];

      //evaluating values usinf all the five stencil points
      dummy[i][j] = (aN + aS + aE + aW - 4.0*aC)/(cp->dx*cp->dy);
    }
  }
  return dummy;
}//END function



/*Prepare the X and Y -axis basis set for FFT computation
NOTE that the double arrays Kx and Ky are modifued within the function*/
void prepare2DFFTBasis(double *Kx, double*Ky, const computeParam *cp)
{
  long nx21 = cp->Nx/2 + 1;
  long nx2 = cp->Nx + 2;

  long ny21 = cp->Ny/2 + 1;
  long ny2 = cp->Ny + 2;

  double delKx = (2.0*PI)/(cp->Nx*cp->dx);
  double delKy = (2.0*PI)/(cp->Ny*cp->dy);

  double tempx = 0.0;
  for(int i = 0; i < nx21; i++){
    tempx = i*delKx;
    Kx[i] = tempx;
    Kx[nx2-i-2] = (-1.0)*tempx;
  }

  double tempy = 0.0;
  for(int i = 0; i < ny21; i++){
    tempy = i*delKy;
    Ky[i] = tempy;
    Ky[ny2-i-2] = (-1.0)*tempy;
  }
} //END Function




/*Prepare the squared version of the X and Y -axis basis set for FFT computation
i.e. K2 = Kx^2 + Ky^2 and K4 = K2^2
NOTE: K2 and K4 are 1D arrays*/
void prepare2DFFTSquaredBasis(double *K2, double *K4, const computeParam *cp)
{
  double *Kx = new double[cp->Nx];
  double *Ky = new double[cp->Ny];

  prepare2DFFTBasis(Kx, Ky, cp); //Prepares the X and Y basis

  int k = 0;
  for(int i = 0; i< cp->Nx; i++){
    for(int j = 0; j< cp->Ny; j++){
      K2[k] = std::pow(Kx[i],2.0) + std::pow(Ky[j],2.0);
      K4[k] = std::pow(K2[k],2.0);
      k++;
    }
  }

  delete []Kx;
  delete []Ky;

}//END function




/* Function returns FFT of a 2D feild contanind real values of the FFT.
Use the get2DREAL function to get a 2D array if real values*/
fftw_complex *get2DFFT(double **Axy, long Nx, long Ny)
{
  //Create the 2D array for FFT
  double **fftData = new double*[3];

  //Allocate memory to the array
  for (int i = 0; i < 3; i++) {
    fftData[i] = new double[3];
  }
  //Create input and output FFTs
  fftw_complex *in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
  fftw_complex *out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nx*Ny);

  //Create plan
  fftw_plan plan = fftw_plan_dft_2d(Nx, Ny, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    // Fill input array
    int k = 0;
    for(int i = 0; i < Nx; i++){
      for(int j = 0; j < Ny; j++){
        in[k][REAL] = Axy[i][j];
        in[k++][IMAG] = 0.0;
      }
    }

  //Execute FFT plan
  fftw_execute(plan);
  //Garbage cleaning
  fftw_destroy_plan(plan);
  fftw_free(in);
  //fftw_free(out);

  return out;
} //End function



/* Function returns the complete inverse FFT of a 2D feild.
Use the get2DREAL function to get a 2D array if real values*/
fftw_complex *get2DiFFT(fftw_complex *inFFT, long Nx, long Ny)
{
  double normalization = Nx*Ny; //Normalization factor

  //Create output FFTs
  fftw_complex *out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
  //Create plan
  fftw_plan plan = fftw_plan_dft_2d(Nx, Ny, inFFT, out, FFTW_BACKWARD, FFTW_ESTIMATE);
  //Execute FFT plan
  fftw_execute(plan);

  //scale the output to the obtain the exact result
  for(int i = 0; i < Nx*Ny; i++){
    out[i][REAL] /= normalization;
    out[i][IMAG] /= normalization;
  }
  //Garbage cleaning
  fftw_destroy_plan(plan);
  return out;
} //End Function



/* Function returns INVERSE FFT of a 2D feild by taking two arrays as inputs
First ayyaryu contains real values while the second contains complex values*/
fftw_complex *get2DiFFTAfterREALManipulation(double **REALxy, double **CMPLXxxy, long Nx, long Ny)
{
  double normalization = Nx*Ny; //Normalization factor
  //Create input and output FFTs
  fftw_complex *in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
  fftw_complex *out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nx*Ny);

  //Create plan
  fftw_plan plan = fftw_plan_dft_2d(Nx, Ny, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

  // Fill input array
  int k = 0;
  for(int i = 0; i < Nx; i++){
    for(int j = 0; j < Ny; j++){
      in[k][REAL]   = REALxy[i][j];
      in[k++][IMAG] = CMPLXxxy[i][j];
    }
  }
  //Execute FFT plan
  fftw_execute(plan);
  //scale the output to the obtain the exact result
  for(int i = 0; i < Nx*Ny; i++){
    out[i][REAL] /= normalization;
    out[i][IMAG] /= normalization;
  }
  //Garbage cleaning
  fftw_destroy_plan(plan);
  fftw_free(in);
  //fftw_free(out);

  return out;
} //End function



/* Function sets up a 2D array of double
NOTE: Be sure to create/allocate memory to the 2D array before passing to this function*/
void set2DiFFTArray(fftw_complex *inFFT, double **ifftArray, long Nx, long Ny)
{
  double normalization = Nx*Ny; //Normalization factor

  //Create output FFTs
  fftw_complex *out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
  //Create plan
  fftw_plan plan = fftw_plan_dft_2d(Nx, Ny, inFFT, out, FFTW_BACKWARD, FFTW_ESTIMATE);
  //Execute FFT plan
  fftw_execute(plan);
  //scale the output to the obtain the exact result
  for(int i = 0; i < Nx*Ny; i++){
    out[i][REAL] /= normalization;
    out[i][IMAG] /= normalization;
  }
  // Copy real portion of the FFT to the new array
  int k = 0;
  for(int i = 0; i < Nx; i++){
    for(int j = 0; j < Ny; j++){
        ifftArray[i][j] =  out[k++][REAL];
    }
  }
  //Garbage cleaning
  fftw_destroy_plan(plan);
  fftw_free(out);
} //End Function



/*This function extracts the REAL values from the FFT, fills up those values in
a 2D array and returns to the calling point*/
double **get2DREAL(fftw_complex *inFFT, long Nx, long Ny)
{
  //Allocate memory to the 2D array
  double **realArray = new double*[Nx];
  for (int i = 0; i < Nx; i++) {
    realArray[i] = new double[Ny];
  }
  // Copy real portion of the FFT to the new array
  int k = 0;
  for(int i = 0; i < Nx; i++){
    for(int j = 0; j < Ny; j++){
        realArray[i][j] =  inFFT[k++][REAL];
    }
  }
  return realArray;
} //END function



/*This function extracts the  COMPLEX values from the FFT, fills up those values in
a 2D array and returns to the calling point*/
double **get2DCMPLX(fftw_complex *inFFT, long Nx, long Ny)
{
  double **cmplxArray = new double*[Nx];
  for (int i = 0; i < Nx; i++) {
    cmplxArray[i] = new double[Ny];
  }
  // Copy real portion of the FFT to the new array
  int k = 0;
  for(int i = 0; i < Nx; i++){
    for(int j = 0; j < Ny; j++){
        cmplxArray[i][j] =  inFFT[k++][IMAG];
    }
  }
  return cmplxArray;
}//END function




void display2DFFTdata(fftw_complex *fftData, long Nx, long Ny)
{
  std::cout << std::endl;
  int k = 0;
  for(int i = 0; i < Nx; i++){
    for(int j = 0; j < Ny; j++){
        std::cout << fftData[k][REAL] << " + i*" << fftData[k][IMAG] <<"   ";
        k++;
    }
    std::cout << std::endl;
  }
} //END Function



/* Set a double 2D array to zero*/
void setZero2DField(double **Fxy, const computeParam *cp){
    for(int i = 0; i < cp->Nx; ++i){
      for(int j = 0; j < cp->Nx; ++j){
        Fxy[i][j] = 0.0;
      }
    }
}//END function



/* Creates a 2D array and sets the value of individual elements to ZERO*/
double **create2DField(const computeParam * cp){
  //Allocate memory for the concentration field
  double **Fxy = new double*[cp->Nx];
  for(int i = 0; i < cp->Nx; i++){
    Fxy[i] = new double[cp->Ny];
  }
  setZero2DField(Fxy, cp);
  return Fxy;
}//END function



/* Write data to VTK file */
void write2DVTKfile(double **Cxy, const computeParam *cp, long tstep)
{
  std::ofstream out;
  std::string fname = std::string("time_") + std::to_string(tstep) + ".vtk";
  out.open(fname);

  long nx, ny, nz;
  nx = cp->Nx;
  ny = cp->Ny;
  nz = 1;

  long nTot = nx*ny*nz; //total number of data or grid points

  //Start writting to file

  //Header information
  out << "# vtk DataFile Version 2.0\n";
  out << "time_10.vtk\n";
  out << "ASCII\n";
  out << "DATASET STRUCTURED_GRID\n";

  //Coordinates of grid points
  out << "DIMENSIONS " << nx << " " << ny << " " << nz <<"\n";
  out << "POINTS " << nTot << " float"<<"\n";

  for(int i = 0; i < nx; i++){
    for(int j = 0; j < ny; j++){
      out << (i-1)*cp->dx << " " << (j-1)*cp->dy << " " << 0.0 <<"\n";
    }
  }//END function

  //Actual grid point values
  out << "POINT_DATA " << nTot <<"\n";
  out << "SCALARS CON float 1\n";
  out << "LOOKUP_TABLE default\n";

  for(int i = 0; i < nx; i++){
    for(int j = 0; j < ny; j++){
      out << Cxy[i][j] << "\n";
    }
  }

  out.close();
}//END function




/* Print file for two arrays contnaing double - basaically x y for plotting
*/
void write1DTXTfile(double *x, double *y, long Nx, std::string filename)
{
  std::ofstream out;
  out.open(filename);
  //std::cout << "Number of items: " << Nx << std::endl;
  for(int i = 0; i < Nx; i++){
    out << x[i] << "  " << y[i] << "\n";
  }
  out.close();
}//END function




/*Funciton for deleting 2D array of typedouble*/
void delete2DArray(double** Arr, long arraySize)
{
  for(int i = 0; i < arraySize; ++i){ //For 2D array of pointers
    delete[] Arr[i];
  }
}//END function

//Display the 2D array oin the screen. MUST be a square matrix
void displayArray(double **Arr, long Nx, long Ny)
{
  std::cout << std::endl;
  for (int i = 0; i < Nx; i++) {
    for (int j = 0; j < Ny; j++) {
      std::cout << Arr[i][j] << " ";
    }
    std::cout << std::endl;
  }
}//END function
















//
