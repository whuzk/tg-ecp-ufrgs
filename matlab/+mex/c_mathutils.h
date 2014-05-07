/*=========================================================================
 * c_mathutils.h
 * 
 *  Title: utilities for manipulation of numerical arrays of integer type
 *  Author:     Diego Sogari
 *  Modified:   05/May/2014+
 *
 *=======================================================================*/
#ifndef C_MATHUTILS
#define C_MATHUTILS

#include <string.h>
#include <math.h>
#include "mex.h"

/*=========================================================================
 * Count zeros of an array
 *=======================================================================*/
int count_zeros(const int *a, mwSize na)
{
    int count = 0;
    for (mwSize i = 0; i < na; i++) {
        if (a[i] == 0) {
            count++;
        }
    }
    return count;
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

#endif