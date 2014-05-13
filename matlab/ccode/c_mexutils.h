/*=========================================================================
 * c_mexutils.h
 * 
 *  Title: utility routines for use with the MATLAB mex library
 *  Author:     Diego Sogari
 *  Modified:   05/May/2014
 *
 *  Notes: supported filter types are: low-pass, high-pass, band-pass,
 *         derivative and moving-average. The design for a given
 *         specification is not always realizable.
 *
 *=======================================================================*/
#ifndef C_MEXUTILS
#define C_MEXUTILS

/*=========================================================================
 * Print an integer vector
 *=======================================================================*/
void mexPrintInt(int *vector, mwSize n)
{
    if (n > 0) {
        for (mwSize i = 0; i < n; i++) {
            mexPrintf("%d ", vector[i]);
        }
        mexPrintf("\n");
    }
}

/*=========================================================================
 * Print a double vector
 *=======================================================================*/
void mexPrintDouble(double *vector, mwSize n)
{
    if (n > 0) {
        for (mwSize i = 0; i < n; i++) {
            mexPrintf("%.4f ", vector[i]);
        }
        mexPrintf("\n");
    }
}

/*=========================================================================
 * Adjust the size of an Int32 mxArray
 *=======================================================================*/
void adjustSizeInt( int *vector, mxArray *mxarray,
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

/*=========================================================================
 * Adjust the size of a Double mxArray
 *=======================================================================*/
void adjustSizeDouble( double *vector, mxArray *mxarray,
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

#endif