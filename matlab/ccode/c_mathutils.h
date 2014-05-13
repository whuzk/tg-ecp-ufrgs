/*=========================================================================
 * c_mathutils.h
 * 
 *  Title: computations using numerical arrays
 *  Author:     Diego Sogari
 *  Modified:   11/May/2014
 *
 *=======================================================================*/
#ifndef C_MATHUTILS
#define C_MATHUTILS

#include <string.h>
#include <math.h>
#include "mex.h"

/*=========================================================================
 * Count zeros of an integer array
 *=======================================================================*/
mwSize count_zeros(const int *a, mwSize na)
{
    mwSize count = 0;
    for (mwSize i = 0; i < na; i++) {
        if (a[i] == 0) {
            count++;
        }
    }
    return count;
}

/*=========================================================================
 * Count zeros of a double array
 *=======================================================================*/
mwSize double_count_zeros(const double *a, mwSize na)
{
    mwSize count = 0;
    for (mwSize i = 0; i < na; i++) {
        if (fpclassify(a[i]) == FP_ZERO) {
            count++;
        }
    }
    return count;
}

/*=========================================================================
 * Perform the cumulative sum of the elements in a.
 *=======================================================================*/
long long sum(const int *a, mwSize na)
{
    long long sum = 0;
    for (mwSize i = 0; i < na; i++) {
        sum += a[i];
    }
    return sum;
}

/*=========================================================================
 * Perform the cumulative sum of the elements in a.
 *=======================================================================*/
double double_sum(const double *a, mwSize na)
{
    double sum = 0.0;
    for (mwSize i = 0; i < na; i++) {
        sum += a[i];
    }
    return sum;
}

/*=========================================================================
 * Perform the element-wise addition of a and b. All arrays must have
 * the same length. The resulting array does not necessarily need to be a
 * different array (it can be either a or b).
 *=======================================================================*/
void add(const int *a, const int *b, int *result, mwSize na)
{
    for (mwSize i = 0; i < na; i++) {
        result[i] = a[i] + b[i];
    }
}

/*=========================================================================
 * Perform the element-wise addition of a and b. All arrays must have
 * the same length. The resulting array does not necessarily need to be a
 * different array (it can be either a or b).
 *=======================================================================*/
void double_add(const double *a, const double *b, double *result,
        mwSize na)
{
    for (mwSize i = 0; i < na; i++) {
        result[i] = a[i] + b[i];
    }
}

/*=========================================================================
 * Perform the element-wise subtraction of a and b. All arrays must have
 * the same length. The resulting array does not necessarily need to be a
 * different array (it can be either a or b).
 *=======================================================================*/
void subtract(const int *a, const int *b, int *result, mwSize na)
{
    for (mwSize i = 0; i < na; i++) {
        result[i] = a[i] - b[i];
    }
}

/*=========================================================================
 * Perform the element-wise subtraction of a and b. All arrays must have
 * the same length. The resulting array does not necessarily need to be a
 * different array (it can be either a or b).
 *=======================================================================*/
void double_subtract(const double *a, const double *b, double *result,
        mwSize na)
{
    for (mwSize i = 0; i < na; i++) {
        result[i] = a[i] - b[i];
    }
}

/*=========================================================================
 * Perform the multiplication of a with scalar b. All arrays must have
 * the same length. The resulting array does not necessarily need to be a
 * different array (it can be a).
 *=======================================================================*/
void scalar_multiply(const int *a, int b, int *result, mwSize na)
{
    for (mwSize i = 0; i < na; i++) {
        result[i] = a[i] * b;
    }
}

/*=========================================================================
 * Perform the multiplication of a with scalar b. All arrays must have
 * the same length. The resulting array does not necessarily need to be a
 * different array (it can be a).
 *=======================================================================*/
void double_scalar_multiply(const double *a, double b, double *result,
        mwSize na)
{
    for (mwSize i = 0; i < na; i++) {
        result[i] = a[i] * b;
    }
}

/*=========================================================================
 * Perform the element-wise multiplication of a and b. All arrays must have
 * the same length. The resulting array does not necessarily need to be a
 * different array (it can be either a or b).
 *=======================================================================*/
void array_multiply(const int *a, const int *b, int *result, mwSize na)
{
    for (mwSize i = 0; i < na; i++) {
        result[i] = a[i] * b[i];
    }
}

/*=========================================================================
 * Perform the element-wise multiplication of a and b. All arrays must have
 * the same length. The resulting array does not necessarily need to be a
 * different array (it can be either a or b).
 *=======================================================================*/
void double_array_multiply(const double *a, const double *b,
        double *result, mwSize na)
{
    for (mwSize i = 0; i < na; i++) {
        result[i] = a[i] * b[i];
    }
}

/*=========================================================================
 * Calculate the root mean squared difference of a and b. All arrays must
 * have the same length.
 *=======================================================================*/
double calcRmsd(const double *a, const double *b, double *aux, mwSize na)
{
    if (na > 0) {
        // a - b
        double_subtract(a, b, aux, na);
        // (a - b) ^ 2
        double_array_multiply(aux, aux, aux, na);
        // sqrt(sum((a - b) ^ 2) / N)
        return sqrt(double_sum(aux, na) / na);
    }
    else return 0;
}

/*=========================================================================
 * Perform the polynomial multiplication of a and b. The resulting array
 * must have a length of (na + nb - 1).
 *=======================================================================*/
void conv(const int *a, mwSize na, const int *b, mwSize nb, int *result)
{
    memset(result, 0, (na + nb - 1) * sizeof(int));
    
    for (mwSize i = 0; i < na; i++) {
        for (mwSize j = 0; j < nb; j++) {
            result[i + j] += a[i] * b[j];
        }
    }
}

/*=========================================================================
 * Perform the polynomial multiplication of a and b. The resulting array
 * must have a length of (na + nb - 1).
 *=======================================================================*/
void double_conv(const double *a, mwSize na, const double *b, mwSize nb,
        double *result)
{
    memset(result, 0, (na + nb - 1) * sizeof(double));
    
    for (mwSize i = 0; i < na; i++) {
        for (mwSize j = 0; j < nb; j++) {
            result[i + j] += a[i] * b[j];
        }
    }
}

/*=========================================================================
 * Perform n-times the polynomial multiplication of a with itself. The
 * resulting array must have a length of 1 + n * (na - 1).
 *=======================================================================*/
void nconv(const int *a, mwSize na, int n, int *result)
{
    result[0] = 1;
    
    if (n > 0) {
        int *temp = (int *)mxMalloc((1 + n * (na - 1)) * sizeof(int));
        mwSize size = 1;
        for (mwSize i = 0; i < n; i++) {
            conv(result, size, a, na, temp);
            size += na - 1;
            memcpy(result, temp, size * sizeof(int));
        }
        mxFree(temp);
    }
}

/*=========================================================================
 * First-order polynomial fitting
 *=======================================================================*/
void first_order_fit(const double *x, const double *y, double *p)
{
    p[0] = (y[1] - y[0]) / (x[1] - x[0]);
    p[1] = y[0] - x[0]*p[1];
}

#endif