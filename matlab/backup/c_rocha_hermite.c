/*=========================================================================
 * c_rocha_hermite.c
 * 
 *  Title: computation of hermite coefficients according to Rocha
 *  Author:     Diego Sogari
 *  Modified:   16/May/2014
 *
 *  Intputs:
 *      1. List of segments 1
 *      2. List of segments 2
 *
 *  Outputs:
 *      1. List of hermite coefficients 1
 *      2. List of hermite coefficients 2
 *
 *  *all inputs and outputs required
 *
 *=======================================================================*/
#include "mex.h"
#include "blas.h"
#include "c_hermcoeff.h"
#include "c_timeutils.h"

/*=========================================================================
 * Abbreviations
 *  LEN     length
 *  BUF     buffer
 *  MIN     minimum
 *  MAX     maximum
 *  SEG     segment
 *  SAMP    sample
 *  HERM    hermite
 *  COEF    coefficient
 *=======================================================================*/

/*=========================================================================
 * Constants
 *=======================================================================*/
#define MIN_INPUTS      2       // minimum number of input arguments
#define MAX_INPUTS      2       // maximum number of input arguments
#define MIN_OUTPUTS     2       // minimum number of output arguments
#define MAX_OUTPUTS     2       // maximum number of output arguments

#define NUM_HERM_COEF   6       // number of hermite coefficients
#define NUM_SEG_SAMP    64      // number of samples in the segment

/*=========================================================================
 * Input and output variables
 *=======================================================================*/
static double *segList1;        // list of segments 1
static double *segList2;        // list of segments 2
static mwSize qrsLen;           // number of segments
static double *outCoeff1;       // output list of hermite coefficients 1
static double *outCoeff2;       // output list of hermite coefficients 2

/*=========================================================================
 * Other variables
 *=======================================================================*/
static mwSize bi;               // current beat index

/*=========================================================================
 * Lookup for the segment lists
 *=======================================================================*/
#define seg1(I,J)   (segList1[(I)*NUM_SEG_SAMP+(J)])
#define seg2(I,J)   (segList2[(I)*NUM_SEG_SAMP+(J)])
#define coef1(I,J)  (outCoeff1[(I)*NUM_HERM_COEF+(J)])
#define coef2(I,J)  (outCoeff2[(I)*NUM_HERM_COEF+(J)])

/*=========================================================================
 * Event triggered for a new beat
 *=======================================================================*/
void onNewBeat()
{
    size_t m = NUM_HERM_COEF;
    size_t p = NUM_SEG_SAMP;
    size_t n = 1;
    char *chn = "N";
    double one = 1.0;
    double zero = 0.0;
    double *A, *B, *C;
    
    A = rocha_S1_hcoef;
    B = &seg1(bi, 0);
    C = &coef1(bi, 0);
    dgemm(chn, chn, &m, &n, &p, &one, A, &m, B, &p, &zero, C, &m);
    
    A = rocha_S2_hcoef;
    B = &seg2(bi, 0);
    C = &coef2(bi, 0);
    dgemm(chn, chn, &m, &n, &p, &one, A, &m, B, &p, &zero, C, &m);
    
    bi++;
}

/*=========================================================================
 * Check correct number of arguments and their types 
 *=======================================================================*/
void checkArgs( int nlhs, mxArray *plhs[],
                int nrhs, const mxArray *prhs[])
{
    mwSize i;
    
    // check for proper number of input arguments
    if (nrhs < MIN_INPUTS || nrhs > MAX_INPUTS) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_rocha_hermite:nrhs",
            "%d input(s) required and %d optional",
            MIN_INPUTS, MAX_INPUTS - MIN_INPUTS);
    }
    // check for proper number of output arguments
    if (nlhs < MIN_OUTPUTS || nlhs > MAX_OUTPUTS) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_rocha_hermite:nlhs",
            "%d output(s) required and %d optional",
            MIN_OUTPUTS, MAX_OUTPUTS - MIN_OUTPUTS);
    }
    // make sure all input arguments are of type double
    for (i = 0; i < nrhs; i++) {
        if (!mxIsDouble(prhs[i]) || mxIsComplex(prhs[i])) {
            mexErrMsgIdAndTxt(
                "EcgToolbox:c_rocha_hermite:notDouble",
                "Input #%d must be of type double.", i + 1);
        }
    }
    // make sure the first input argument has the right dimensions
    if (mxGetM(prhs[0]) != NUM_SEG_SAMP) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_rocha_hermite:badDimensions",
            "Input #1 must have %d lines.", NUM_SEG_SAMP);
    }
    // make sure the second input argument has the right dimensions
    if (mxGetM(prhs[1]) != NUM_SEG_SAMP) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_rocha_hermite:badDimensions",
            "Input #2 must have %d lines.", NUM_SEG_SAMP);
    }
}

/*=========================================================================
 * Hanlde input arguments 
 *=======================================================================*/
void handleInputs( int nrhs, const mxArray *prhs[],
                   mwSize *nrows, mwSize *ncols)
{
    // get pointers to the data in the input vectors
    segList1 = mxGetPr(prhs[0]);
    segList2 = mxGetPr(prhs[1]);
    
    // get the dimensions of the first input
    *nrows = (mwSize)mxGetM(prhs[0]);
    *ncols = (mwSize)mxGetN(prhs[0]);
}

/*=========================================================================
 * Hanlde output arguments 
 *=======================================================================*/
void handleOutputs(int nlhs, mxArray *plhs[])
{
    // get a pointer to the first output
    plhs[0] = mxCreateDoubleMatrix(NUM_HERM_COEF, qrsLen, mxREAL);
    outCoeff1 = mxGetPr(plhs[0]);
    
    // get a pointer to the second output
    plhs[1] = mxCreateDoubleMatrix(NUM_HERM_COEF, qrsLen, mxREAL);
    outCoeff2 = mxGetPr(plhs[1]);
}

/*=========================================================================
 * The main routine 
 *=======================================================================*/
void doTheJob()
{
    double time;
    mwSize i;
    
    // start time counter
    tic();
    
    // initialize beat index
    bi = 0;
    
    // process one input sample at a time
    for (i = 0; i < qrsLen; i++) {
        onNewBeat();
    }
    
    // stop time counter
    time = toc();
    
    // display time statistics
    mexPrintf("Total processing time: %.2f ms\n", 1000 * time);
    mexPrintf("Average time per beat: %.2f ns\n", 1000000000 * time / qrsLen);
}

/*=========================================================================
 * The gateway function 
 *=======================================================================*/
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    mwSize nrows;
    mwSize ncols;
    
    // check argument correctness
    checkArgs(nlhs, plhs, nrhs, prhs);
    
    // handle input arguments
    handleInputs(nrhs, prhs, &nrows, &ncols);
    
    // get the number of beats in the list
    qrsLen = ncols;
    
    // handle output arguments
    handleOutputs(nlhs, plhs);
    
    // do the actual processing
    doTheJob();
}
