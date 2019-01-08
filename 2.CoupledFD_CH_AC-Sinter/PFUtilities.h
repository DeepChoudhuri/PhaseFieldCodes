#ifndef GUARD_PFUtilities_h
#define GUARD_PFUtilities_h

#include<iostream>
#include<fstream>
#include<fftw3.h>
#include<cmath>
#include "stdlib.h"
#include "time.h"
#include "Parameters.h"

#define REAL 0
#define IMAG 1

const double PI = 3.141592653589793238463;

//Functions related to mathematical manipulation and processing in 2D
double **get2DLaplacian(double **, const computeParam *);

void prepare2DFFTBasis(double *, double*, const computeParam *);
void prepare2DFFTSquaredBasis(double *, double *, const computeParam *cp);

fftw_complex *get2DFFT(double **, long, long);
fftw_complex *get2DiFFT(fftw_complex *, long, long);
fftw_complex *get2DiFFTAfterREALManipulation(double **, double **, long, long);
void set2DiFFTArray(fftw_complex *, double **, long, long);
void display2DFFTdata(fftw_complex *, long, long);

double **get2DREAL(fftw_complex *, long, long);
double **get2DCMPLX(fftw_complex *, long, long);

// Fuctions related to Array memory allocation
double **create2DField(const computeParam *);
void setZero2DField(double **, const computeParam *);
void delete2DArray(double**, long);
void displayArray(double **, const computeParam *);

// Functions related to file I/O
void write1DTXTfile(double *, double *, long, std::string);
void write2DVTKfile(double **, const computeParam *, long);
void write2DVTKfile(double **, double **, const computeParam *, long);

#endif
