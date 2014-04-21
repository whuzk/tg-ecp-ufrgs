/*=========================================================================
 * c_ecg_utils.c
 *=======================================================================*/
#include "c_ecg_utils.h"
#include <math.h>

/*=========================================================================
 * FLOATING-POINT FUNCTIONS
 *=======================================================================*/
/* Iterative estimator */
double estimate(double x, double r, double v)
{
    return (1-r)*x + r*v;
}

/* Limit a value by lower and upper bounds */
double limit(double x, double a, double b)
{
    return fmin(fmax(x,a),b);
}

/* Find the position of the maximum value */
mwSize findmax(const double *x, mwSize nx)
{
    mwSize pos = 0;
    double y = x[0];
    
    for (mwSize i = 1; i < nx; i++) {
        if (x[i] > y) {
            y = x[i];
            pos = i;
        }
    }
    return pos;
}

/* Get the maximum slope of the signal */
double maxdiff(const double *x, mwSize nx)
{
    double d = 0;
    for (mwSize i = 1; i < nx; i++) {
        if (x[i]-x[i-1] > d) {
            d = x[i]-x[i-1];
        }
    }
    return d;
}

/* Multiply and accumulate */
double multacc(const double *x, const double *h, mwSize start, mwSize end)
{
    double acc = 0;
    for (mwSize k = start; k < end; k++) {
        acc += h[k] * x[-k];
    }
    return acc;
}

/* Multiply together */
double multmult(const double *x, mwSize start, mwSize end)
{
    double prod = 1.0;
    for (mwSize k = start; k < end; k++) {
        prod *= x[-k];
    }
    return prod;
}

/* Normal convolution without delay */
void convolve( const double *x, mwSize nx,
               const double *h, mwSize nh, double *y)
{
    mwSize offset = nh/2;
    
    for (mwSize i = offset; i < nh; i++) {
        y[i-offset] = multacc(x+i, h, 0, i+1);
    }
    for (mwSize i = nh; i < nx; i++) {
        y[i-offset] = multacc(x+i, h, 0, nh);
    }
    for (mwSize i = nx; i < nx+offset; i++) {
        y[i-offset] = multacc(x+i, h, i-nx+1, nh);
    }
}

/* FIR filtering */
void fir( const double *x, mwSize nx,
          const double *h, mwSize nh, double *y)
{
    for (mwSize i = 0; i < nh; i++) {
        y[i] = multacc(x+i, h, 0, i+1);
    }
    for (mwSize i = nh; i < nx; i++) {
        y[i] = multacc(x+i, h, 0, nh);
    }
}

/* Non-linear convolution without delay */
void nlconvolve(const double *x, mwSize nx, int m, double *y)
{
    mwSize offset = m/2;
    
    for (mwSize i = offset; i < m; i++) {
        y[i-offset] = multmult(x+i, 0, i+1);
    }
    for (mwSize i = m; i < nx; i++) {
        y[i-offset] = multmult(x+i, 0, m);
    }
    for (mwSize i = nx; i < nx+offset; i++) {
        y[i-offset] = multmult(x+i, i-nx+1, m);
    }
}

/* Non-linear FIR filtering */
void nlfir( const double *x, mwSize nx, int m, double *y)
{
    for (mwSize i = 0; i < m; i++) {
        y[i] = multmult(x+i, 0, i+1);
    }
    for (mwSize i = m; i < nx; i++) {
        y[i] = multmult(x+i, 0, m);
    }
}

/* Absolute value */
void absolute(double *x, mwSize nx)
{
    for (mwSize i = 0; i < nx; i++) {
        x[i] = fabs(x[i]);
    }
}

/*=========================================================================
 * INTEGER FUNCTIONS
 *=======================================================================*/
/* Find the position of the maximum value */
mwSize ifindmax(const int *x, mwSize nx)
{
    mwSize pos = 0;
    int y = x[0];
    
    for (mwSize i = 1; i < nx; i++) {
        if (x[i] > y) {
            y = x[i];
            pos = i;
        }
    }
    return pos;
}

/* Get the maximum slope of the signal */
int imaxdiff(const int *x, mwSize nx)
{
    int d = 0;
    for (mwSize i = 1; i < nx; i++) {
        if (x[i]-x[i-1] > d) {
            d = x[i]-x[i-1];
        }
    }
    return d;
}