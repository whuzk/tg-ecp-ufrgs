/*=========================================================================
 * c_ecg_utils.h
 *=======================================================================*/
#ifndef C_ECG_UTILS
#define C_ECG_UTILS

#include "mex.h"

/*=========================================================================
 * FLOATING-POINT FUNCTIONS
 *=======================================================================*/
double estimate(double x, double r, double v);
double limit(double x, double a, double b);
mwSize findmax(const double *x, mwSize nx);
double maxdiff(const double *x, mwSize nx);
double multacc(const double *x, const double *h, mwSize start, mwSize end);
double multmult(const double *x, mwSize start, mwSize end);
void convolve(const double *x, mwSize nx, const double *h, mwSize nh, double *y);
void fir(const double *x, mwSize nx, const double *h, mwSize nh, double *y);
void nlconvolve(const double *x, mwSize nx, int m, double *y);
void nlfir( const double *x, mwSize nx, int m, double *y);
void absolute(double *x, mwSize nx);

/*=========================================================================
 * INTEGER FUNCTIONS
 *=======================================================================*/
mwSize ifindmax(const int *x, mwSize nx);
int imaxdiff(const int *x, mwSize nx);

#endif