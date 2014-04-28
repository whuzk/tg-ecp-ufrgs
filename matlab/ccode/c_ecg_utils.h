/*=========================================================================
 * c_ecg_utils.h
 *=======================================================================*/
#ifndef C_ECG_UTILS
#define C_ECG_UTILS

#include <math.h>
#include "mex.h"

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
void nlfir(const double *x, mwSize nx, int m, double *y)
{
    for (mwSize i = 0; i < m; i++) {
        y[i] = multmult(x+i, 0, i+1);
    }
    for (mwSize i = m; i < nx; i++) {
        y[i] = multmult(x+i, 0, m);
    }
}

/* Max-filter */
void maxfilter(const double *x, mwSize nx, mwSize m, double *y, double ai)
{
    mwSize *pos = (mwSize*)mxMalloc(m*sizeof(mwSize));
    double *val = (double*)mxMalloc(m*sizeof(double));
    mwSize first = 0;
    mwSize len = 1;
    mwSize j,idx;
    double a;
    
    // initial max value (this can be ignored by providing -HUGE_VAL)
    val[0] = ai;
    pos[0] = -1;
    
    for (mwSize i = 0; i < nx; i++) {
        // get the current sample
        a = x[i];
        // search for the first element greater than the current sample
        j = len;
        while (j > 0 && a >= val[(first + j - 1) % m]) {
            j = j - 1;
        }
        // put the sample next to element found and adjust the length
        idx = (first + j) % m;
        val[idx] = a;
        pos[idx] = i;
        len = j + 1;
        // check if the first in line has gone out of the windows length
        if (len > m || pos[first] == i - m) {
            len = len - 1;
            first = (first + 1) % m;
        }
        // store the result
        y[i] = val[first];
    }
    // deallocate memory
    mxFree(pos);
    mxFree(val);
}

/* Adjust the size of an mxArray */
void adjustSize( double *vector, mxArray *mxarray,
                 mwSize nrows, mwSize ncols, mwSize newlen)
{
    if (vector != NULL && mxarray != NULL) {
        if (ncols == 1) {
            mxSetM(mxarray, newlen);
        }
        else {
            mxSetN(mxarray, newlen);
        }
        mxSetPr(mxarray, mxRealloc(vector, newlen * sizeof(double)));
    }
}

/*=========================================================================
 * INTEGER FUNCTIONS
 *=======================================================================*/
/* Iterative estimator */
int iestimate(int x, int log2r, int v)
{
    return (((1 << log2r) - 1) * x + v) >> log2r;
}

/* Limit a value by lower and upper bounds */
int ilimit(int x, int a, int b)
{
    return min(max(x, a), b);
}

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

/* Adjust the size of an mxArray */
void iadjustSize( int *vector, mxArray *mxarray,
                 mwSize nrows, mwSize ncols, mwSize newlen)
{
    if (vector != NULL && mxarray != NULL) {
        if (ncols == 1) {
            mxSetM(mxarray, newlen);
        }
        else {
            mxSetN(mxarray, newlen);
        }
        mxSetPr(mxarray, mxRealloc(vector, newlen * sizeof(int)));
    }
}

#endif