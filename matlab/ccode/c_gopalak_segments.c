/*=========================================================================
 * c_gopalak_segments.c
 * 
 *  Title: extraction of segments according to Gopalakrishnan
 *  Author:     Diego Sogari
 *  Modified:   16/May/2014
 *
 *  Intputs:
 *      1. List of beats (in matrix form)
 *      2. List of RR intervals
 *
 *  Outputs:
 *      1. List of segments
 *
 *  *all inputs and outputs required
 *
 *=======================================================================*/
#include <string.h>
#include <math.h>
#include "mex.h"
#include "c_upfirdn.h"
#include "c_timeutils.h"

/*=========================================================================
 * Abbreviations
 *  LEN     length
 *  BUF     buffer
 *  MIN     minimum
 *  MAX     maximum
 *  DEF     default
 *  FREQ    frequency
 *=======================================================================*/

/*=========================================================================
 * Constants
 *=======================================================================*/
#define MIN_INPUTS      2       // minimum number of input arguments
#define MAX_INPUTS      2       // maximum number of input arguments
#define MIN_OUTPUTS     1       // minimum number of output arguments
#define MAX_OUTPUTS     1       // maximum number of output arguments

#define OUT_SEG_SIZE    250     // size of the output segment (in samples)

/*=========================================================================
 * Input and output variables
 *=======================================================================*/
static double *beatList;        // list of beats
static double *rrList;          // list of RR intervals
static mwSize qrsLen;           // number of beats
static double *outSeg;          // output list of segments

/*=========================================================================
 * Other variables
 *=======================================================================*/
static mwSize bi;               // current beat index
static mwSize frameSize;        // size of one beat
static double filterH[OUT_SEG_SIZE] = {1.0};    // filter impulse response

/*=========================================================================
 * Lookup for vectors
 *=======================================================================*/
#define beat(I,J)   (beatList[(I)*frameSize+(J)])
#define seg(I,J)    (outSeg[(I)*OUT_SEG_SIZE+(J)])

/*=========================================================================
 * Event triggered for a new beat
 *=======================================================================*/
void onNewBeat()
{
    mwSize rr = min(frameSize, (mwSize)rrList[bi]);
    mwSize istart = (frameSize - rr) >> 1;
    mwSize n = OUT_SEG_SIZE;
    double *x = &beat(bi, istart);
    double *y = &seg(bi, 0);
    
    upfirdn(y, n, x, rr, filterH, n, n, rr);
    
    bi++;
}

/*=========================================================================
 * Check correct number of arguments and their types 
 *=======================================================================*/
void checkArgs( int nlhs, mxArray *plhs[],
                int nrhs, const mxArray *prhs[])
{
    // check for proper number of input arguments
    if (nrhs < MIN_INPUTS || nrhs > MAX_INPUTS) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_gopalak_segments:nrhs",
            "%d input(s) required and %d optional",
            MIN_INPUTS, MAX_INPUTS - MIN_INPUTS);
    }
    // check for proper number of output arguments
    if (nlhs < MIN_OUTPUTS || nlhs > MAX_OUTPUTS) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_gopalak_segments:nlhs",
            "%d output(s) required and %d optional",
            MIN_OUTPUTS, MAX_OUTPUTS - MIN_OUTPUTS);
    }
    // make sure all input arguments are of type double
    for (mwSize i = 0; i < nrhs; i++) {
        if (!mxIsDouble(prhs[i]) || mxIsComplex(prhs[i])) {
            mexErrMsgIdAndTxt(
                "EcgToolbox:c_gopalak_segments:notDouble",
                "Input #%d must be of type double.", i + 1);
        }
    }
    // make sure the first input argument has more than one line
    if (mxGetM(prhs[0]) == 1) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_gopalak_segments:badDimensions",
            "Input #1 must have more than one line.");
    }
    // make sure the second input argument has exactly one column
    if (mxGetN(prhs[1]) != 1) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_gopalak_segments:badDimensions",
            "Input #2 must have exactly one column.");
    }
}

/*=========================================================================
 * Hanlde input arguments 
 *=======================================================================*/
void handleInputs( int nrhs, const mxArray *prhs[],
                   mwSize *nrows, mwSize *ncols)
{
    // get pointers to the data in the input vectors
    beatList = mxGetPr(prhs[0]);
    rrList = mxGetPr(prhs[1]);
    
    // get the dimensions of the first input
    *nrows = (mwSize)mxGetM(prhs[0]);
    *ncols = (mwSize)mxGetN(prhs[0]);
    
    // make sure the second input argument has compatible dimensions
    if (mxGetM(prhs[1]) != *ncols) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_gopalak_segments:badDimensions",
            "Inputs #1 and #2 must have compatible dimensions.");
    }
}

/*=========================================================================
 * Hanlde output arguments 
 *=======================================================================*/
void handleOutputs(int nlhs, mxArray *plhs[])
{
    // get a pointer to the first output
    plhs[0] = mxCreateDoubleMatrix(OUT_SEG_SIZE, qrsLen, mxREAL);
    outSeg = mxGetPr(plhs[0]);
}

/*=========================================================================
 * The main routine 
 *=======================================================================*/
void doTheJob()
{
    double time;
    
    // start time counter
    tic();
    
    // initialize beat index
    bi = 0;
    
    // process one input sample at a time
    for (mwSize i = 0; i < qrsLen; i++) {
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
    
    // get the frame size
    frameSize = nrows;
    
    // get the number of beats in the list
    qrsLen = ncols;
    
    // handle output arguments
    handleOutputs(nlhs, plhs);
    
    // do the actual processing
    doTheJob();
}
