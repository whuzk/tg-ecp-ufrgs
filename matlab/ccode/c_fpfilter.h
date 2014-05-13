/*=========================================================================
 * c_fpfilter.h
 * 
 *  Title: utilities for the use of filters with floating-point
 *         coefficients.
 *  Author:     Diego Sogari
 *  Modified:   10/May/2014
 *
 *=======================================================================*/
#ifndef C_FPFILTER
#define C_FPFILTER

#include <string.h>
#include <math.h>
#include "mex.h"
#include "c_mathutils.h"

/*=========================================================================
 * Type definitions
 *=======================================================================*/
typedef struct {
    mwSize size;        // size of the coefficient array
    double *val;        // array of coefficient values
    mwSize *idx;        // indices of the corresponding coefficients
                        //      in the filter's recurrence equation
} fpfcoeffs;

typedef struct {
    mwSize size;        // size of the buffer array
    double *val;        // array of buffer values
} fpfbuffer;

typedef struct {
    fpfcoeffs b;        // numerator coefficients
    fpfcoeffs a;        // denominator coefficients
    fpfbuffer x;        // X buffer of the filter
    fpfbuffer y;        // Y buffer of the filter
} fpfobject;

/*=========================================================================
 * Macros
 *=======================================================================*/
#define xval(f,k)    (f).x.val[(k)&((f).x.size-1)]
#define yval(f,k)    (f).y.val[(k)&((f).y.size-1)]

#define initfpfilter(f) \
    (f).b.size = 0;     \
    (f).b.val = NULL;   \
    (f).a.size = 0;     \
    (f).a.val = NULL;   \
    (f).x.size = 0;     \
    (f).x.val = NULL;   \
    (f).y.size = 0;     \
    (f).y.val = NULL;
    
#define endfpfilter(f)  \
    mxFree((f).b.val);  \
    mxFree((f).a.val);  \
    mxFree((f).x.val);  \
    mxFree((f).y.val);

/*=========================================================================
 * Update filter memory with an incoming sample, and return the output.
 *=======================================================================*/
double fpfnewx(fpfobject *filter, unsigned int ci, double sample)
{
    double result = 0;
    
    // update filter X memory
    xval(*filter, ci) = sample;
    
    // compute the first part of the recurrence equation
    for (mwSize i = 0; i < filter->b.size; i++) {
        result += filter->b.val[i] * xval(*filter, ci - filter->b.idx[i]);
    }
    
    // compute the second part of the recurrence equation
    for (mwSize i = 1; i < filter->a.size; i++) {
        result -= filter->a.val[i] * yval(*filter, ci - filter->a.idx[i]);
    }
    
    // update filter Y memory
    yval(*filter, ci) = result;
    
    return result;
}

/*=========================================================================
 * Fill the coefficient object with the non-zero portion of an array. If
 * the number of non-zero elements in the array is different from the size
 * of the coefficient object, this size is adjusted accordingly. The index
 * of the last non-zero coefficient is returned, or zero if no such
 * coefficient exists.
 *=======================================================================*/
mwSize construct_coeff(fpfcoeffs *coeffObj, const double *a, mwSize na)
{
    mwSize count = na - double_count_zeros(a, na);
    
    // adjust the size of the coefficient vectors
    if (count != coeffObj->size) {
        coeffObj->val = (double *)mxRealloc(coeffObj->val, count * sizeof(double));
        coeffObj->idx = (mwSize *)mxRealloc(coeffObj->idx, count * sizeof(mwSize));
        coeffObj->size = count;
    }
    
    // fill in the non-zero portion of the array
    while (count-- > 0) {
        while (fpclassify(a[--na]) == FP_ZERO);
        coeffObj->val[count] = a[na];
        coeffObj->idx[count] = na;
    }
    
    return (coeffObj->size > 0) ? coeffObj->idx[coeffObj->size - 1] : 0;
}

/*=========================================================================
 * Initialize a filter object and allocate memory for its properties
 *=======================================================================*/
void create_filter(fpfobject *filter, const double *b, mwSize nb,
        const double *a, mwSize na)
{
    mwSize xsize, ysize;
    
    // initialize the coefficient objects
    xsize = construct_coeff(&filter->b, b, nb) + 1;
    ysize = construct_coeff(&filter->a, a, na) + 1;
    
    // initialize the X buffer
    if (xsize > 2) {
        xsize = 1 << (1 + ilogb(xsize - 1));
    }
    filter->x.val = (double *)mxRealloc(filter->x.val, xsize * sizeof(double));
    memset(filter->x.val, 0, xsize * sizeof(double));
    filter->x.size = xsize;
    
    // initialize the Y buffer
    if (ysize > 2) {
        ysize = 1 << (1 + ilogb(ysize - 1));
    }
    filter->y.val = (double *)mxRealloc(filter->y.val, ysize * sizeof(double));
    memset(filter->y.val, 0, ysize * sizeof(double));
    filter->y.size = ysize;
}

/*=========================================================================
 * Create an all-pass filter with the delay given in samples (gain = 1)
 *=======================================================================*/
bool create_allpass(fpfobject *filter, int delay)
{
    double *b;
    
    // check argument correctness
    if (delay < 0) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_fpfilter:design_allpass:noSolution",
            "could not design an all-pass for this specfication");
        return false;
    }
    
    // create the numerator of the transfer function
    b = (double *)mxCalloc(delay + 1, sizeof(double));
    b[delay] = 1.0;
    
    // initialize the filter object
    create_filter(filter, b, delay + 1, NULL, 0);
    
    // deallocate memory
    mxFree(b);
    
    return true;
}

#endif