/*=========================================================================
 * c_pang_jpoints.c
 * 
 *  Title: detection of J points using the method of Pang
 *  Author:     Diego Sogari
 *  Modified:   16/May/2014
 *
 *  Intputs:
 *      1. List of RR intervals
 *      2. Location of R peak (center of the beat frame)
 *      3. Sampling frequency (in Hertz)
 *
 *  Outputs:
 *      1. List of J points
 *
 *  *all inputs and outputs required
 *
 *=======================================================================*/
#include <math.h>
#include "mex.h"
#include "c_timeutils.h"

/*=========================================================================
 * Abbreviations
 *  LEN     length
 *  BUF     buffer
 *  MIN     minimum
 *  MAX     maximum
 *  DEF     default
 *  FREQ    frequency
 *  IDX     index
 *  HR      heart rate
 *=======================================================================*/

/*=========================================================================
 * Constants
 *=======================================================================*/
#define MIN_INPUTS      3       // minimum number of input arguments
#define MAX_INPUTS      3       // maximum number of input arguments
#define MIN_OUTPUTS     1       // minimum number of output arguments
#define MAX_OUTPUTS     1       // maximum number of output arguments

#define NUM_HR_LIMITS   3       // number of heart rate limits
#define NUM_MS_POINTS   NUM_HR_LIMITS+1 // number of measuring points

/*=========================================================================
 * Input and output variables
 *=======================================================================*/
static double *rrList;          // list of RR intervals
static mwSize qrsLen;           // number of beats
static mwSize rPeak;            // location of R peak
static double sampFreq;         // sampling frequency
static double *outJay;          // output list of J points

/*=========================================================================
 * Filter and buffer variables
 *=======================================================================*/
static mwSize hrLimits[] = {100, 110, 120};     // heart rate limits
static double msPoints[] = {0.12, 0.112, 0.104, 0.1};   // measuring points
static mwSize bi;               // current beat index

/*=========================================================================
 * Event triggered for a new beat
 *=======================================================================*/
void onNewBeat()
{
    mwSize i, heartRate = rrList[bi] * 60 / sampFreq;
    
    for (i = 0; i < NUM_HR_LIMITS && heartRate >= hrLimits[i]; i++);
    outJay[bi++] = rPeak + (int)(msPoints[i] * sampFreq);
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
            "EcgToolbox:c_pang_jpoints:nrhs",
            "%d input(s) required and %d optional",
            MIN_INPUTS, MAX_INPUTS - MIN_INPUTS);
    }
    // check for proper number of output arguments
    if (nlhs < MIN_OUTPUTS || nlhs > MAX_OUTPUTS) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_pang_jpoints:nlhs",
            "%d output(s) required and %d optional",
            MIN_OUTPUTS, MAX_OUTPUTS - MIN_OUTPUTS);
    }
    // make sure all input arguments are of type double
    for (i = 0; i < nrhs; i++) {
        if (!mxIsDouble(prhs[i]) || mxIsComplex(prhs[i])) {
            mexErrMsgIdAndTxt(
                "EcgToolbox:c_pang_jpoints:notDouble",
                "Input #%d must be of type double.", i + 1);
        }
    }
    // make sure the first input argument is a vector
    if (mxGetM(prhs[0]) != 1 && mxGetN(prhs[0]) != 1) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_pang_jpoints:notVector",
            "Input #1 must be a vector.");
    }
    // make sure the remaining input arguments are all scalars
    for (i = 1; i < nrhs; i++) {
        if (mxGetNumberOfElements(prhs[i]) != 1) {
            mexErrMsgIdAndTxt(
                "EcgToolbox:c_pang_jpoints:notScalar",
                "Input #%d must be a scalar.", i + 1);
        }
    }
}

/*=========================================================================
 * Hanlde input arguments 
 *=======================================================================*/
void handleInputs( int nrhs, const mxArray *prhs[],
                   mwSize *nrows, mwSize *ncols)
{
    // get pointers to the data in the input vectors
    rrList = mxGetPr(prhs[0]);
    
    // get the dimensions of the first input
    *nrows = (mwSize)mxGetM(prhs[0]);
    *ncols = (mwSize)mxGetN(prhs[0]);
    
    // get the location of the R peak
    rPeak = (int)mxGetScalar(prhs[1]) - 1;
    
    // get the sampling frequency
    sampFreq = mxGetScalar(prhs[2]);
}

/*=========================================================================
 * Hanlde output arguments 
 *=======================================================================*/
void handleOutputs(int nlhs, mxArray *plhs[])
{
    // get a pointer to the first output
    plhs[0] = mxCreateDoubleMatrix(qrsLen, 1, mxREAL);
    outJay = mxGetPr(plhs[0]);
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
    qrsLen = max(nrows,ncols);
    
    // handle output arguments
    handleOutputs(nlhs, plhs);
    
    // do the actual processing
    doTheJob();
}
