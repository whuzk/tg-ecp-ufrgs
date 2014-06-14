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

#define PI          3.1415926535897932384626433832795
#define ILOG2(x)    (int)(log10((double)x) / log10(2.0))
#define SIGN(x)     (((double)x) >= 0.0 ? 0 : 1)

/*=========================================================================
 * Count zeros of an integer array
 *=======================================================================*/
mwSize count_zeros(const int *a, mwSize na)
{
    mwSize i, count = 0;
    for (i = 0; i < na; i++) {
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
    mwSize i, count = 0;
    for (i = 0; i < na; i++) {
        if (a[i] == 0.0) {
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
    mwSize i;
    
    for (i = 0; i < na; i++) {
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
    mwSize i;
    
    for (i = 0; i < na; i++) {
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
    mwSize i;
    
    for (i = 0; i < na; i++) {
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
    mwSize i;
    
    for (i = 0; i < na; i++) {
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
    mwSize i;
    
    for (i = 0; i < na; i++) {
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
    mwSize i;
    
    for (i = 0; i < na; i++) {
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
    mwSize i;
    
    for (i = 0; i < na; i++) {
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
    mwSize i;
    
    for (i = 0; i < na; i++) {
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
    mwSize i;
    
    for (i = 0; i < na; i++) {
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
    mwSize i;
    
    for (i = 0; i < na; i++) {
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
    mwSize i, j;
    
    memset(result, 0, (na + nb - 1) * sizeof(int));
    
    for (i = 0; i < na; i++) {
        for (j = 0; j < nb; j++) {
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
    mwSize i, j;
    
    memset(result, 0, (na + nb - 1) * sizeof(double));
    
    for (i = 0; i < na; i++) {
        for (j = 0; j < nb; j++) {
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
        mwSize i, size;
        int *temp = (int *)mxMalloc((1 + n * (na - 1)) * sizeof(int));
        
        size = 1;
        for (i = 0; i < n; i++) {
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
    p[1] = y[0] - x[0] * p[0];
}

/*=========================================================================
 * Rational approximation of x by p/q, with the specified tolerance
 *=======================================================================*/
void rational(double x, int *p, int *q, double tol)
{
    int nh[2] = {1, 0};    // [n(k) n(k-1)]
    int dh[2] = {0, 1};    // [d(k) d(k-1)]
    double savex, diff;
    int k, d, saven, saved;
    
    k = 0;
    savex = x;
    while (true) {
        k++;
        d = (int)x;
        x = x - d;
        
        saven = nh[0];
        saved = dh[0];
        nh[0] = nh[0] * d + nh[1];
        dh[0] = dh[0] * d + dh[1];
        nh[1] = saven;
        dh[1] = saved;
        
        diff = fabs(nh[0] / (double)dh[0] - savex);
        if (k == 100 || x == 0 || diff <= tol) {
            break;
        }
        x = 1 / x;
    }
    if (dh[0] < 0) {
        *p = -nh[0];
    }
    else *p = nh[0];
    *q = abs(dh[0]);
}

#endif