/*=========================================================================
 * c_rocha_segments.c
 * 
 *  Title: extraction of segments according to Rocha
 *  Author:     Diego Sogari
 *  Modified:   16/May/2014
 *
 *  Intputs:
 *      1. List of beats (in matrix form)
 *      2. List of starting points
 *      3. List of J points
 *      4. List of ending points
 *
 *  Outputs:
 *      1. List of segments 1
 *      2. List of segments 2
 *
 *  *all inputs and outputs required
 *
 *=======================================================================*/
#include <string.h>
#include <math.h>
#include "mex.h"
#include "c_upfirdn.h"
#include "c_mathutils.h"
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
#define MIN_INPUTS      4       // minimum number of input arguments
#define MAX_INPUTS      4       // maximum number of input arguments
#define MIN_OUTPUTS     2       // minimum number of output arguments
#define MAX_OUTPUTS     2       // maximum number of output arguments

#define OUT_SEG_SIZE    64      // size of the output segment (in samples)

/*=========================================================================
 * Input and output variables
 *=======================================================================*/
static double *beatList;        // list of beats
static double *begList;         // list of starting points
static double *jayList;         // list of J points
static double *endList;         // list of ending points
static mwSize qrsLen;           // number of beats
static double *outSeg1;         // output list of segments
static double *outSeg2;         // output list of segments

/*=========================================================================
 * Other variables
 *=======================================================================*/
static mwSize bi;               // current beat index
static mwSize frameSize;        // size of one beat
static double filterH[OUT_SEG_SIZE];    // filter impulse response

/*=========================================================================
 * Lookup for vectors
 *=======================================================================*/
#define beat(I,J)   (beatList[(I)*frameSize+(J)])
#define seg1(I,J)   (outSeg1[(I)*OUT_SEG_SIZE+(J)])
#define seg2(I,J)   (outSeg2[(I)*OUT_SEG_SIZE+(J)])

/*=========================================================================
 * Event triggered for a new beat
 *=======================================================================*/
void onNewBeat()
{
    mwSize istart, iend, len;
    mwSize n = OUT_SEG_SIZE;
    double *x, *y;
    int p, q;
    
    istart = (mwSize)begList[bi] - 1;
    iend = (mwSize)jayList[bi] - 1;
    if (istart < iend) {
        len = iend - istart + 1;
        x = &beat(bi, istart);
        y = &seg1(bi, 0);
        rational(n / (double)len, &p, &q, 1.0e-5);
        upfirdn(y, n, x, len, filterH, p, p, q);
    }
    
    istart = (mwSize)jayList[bi] - 1;
    iend = (mwSize)endList[bi] - 1;
    if (istart < iend) {
        len = iend - istart + 1;
        x = &beat(bi, istart);
        y = &seg2(bi, 0);
        rational(n / (double)len, &p, &q, 1.0e-5);
        upfirdn(y, n, x, len, filterH, p, p, q);
    }
    
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
            "EcgToolbox:c_rocha_segments:nrhs",
            "%d input(s) required and %d optional",
            MIN_INPUTS, MAX_INPUTS - MIN_INPUTS);
    }
    // check for proper number of output arguments
    if (nlhs < MIN_OUTPUTS || nlhs > MAX_OUTPUTS) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_rocha_segments:nlhs",
            "%d output(s) required and %d optional",
            MIN_OUTPUTS, MAX_OUTPUTS - MIN_OUTPUTS);
    }
    // make sure all input arguments are of type double
    for (mwSize i = 0; i < nrhs; i++) {
        if (!mxIsDouble(prhs[i]) || mxIsComplex(prhs[i])) {
            mexErrMsgIdAndTxt(
                "EcgToolbox:c_rocha_segments:notDouble",
                "Input #%d must be of type double.", i + 1);
        }
    }
    // make sure the first input argument has more than one line
    if (mxGetM(prhs[0]) == 1) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_rocha_segments:badDimensions",
            "Input #1 must have more than one line.");
    }
    // make sure the remaining input arguments have exactly one column
    for (mwSize i = 1; i < nrhs; i++) {
        if (mxGetN(prhs[i]) != 1) {
            mexErrMsgIdAndTxt(
                "EcgToolbox:c_rocha_segments:badDimensions",
                "Input #%d must have exactly one column.", i + 1);
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
    beatList = mxGetPr(prhs[0]);
    begList = mxGetPr(prhs[1]);
    jayList = mxGetPr(prhs[2]);
    endList = mxGetPr(prhs[3]);
    
    // get the dimensions of the first input
    *nrows = (mwSize)mxGetM(prhs[0]);
    *ncols = (mwSize)mxGetN(prhs[0]);
    
    // make sure the input arguments 2-4 have compatible dimensions
    for (mwSize i = 1; i < nrhs; i++) {
        if (mxGetM(prhs[i]) != *ncols) {
            mexErrMsgIdAndTxt(
                "EcgToolbox:c_rocha_segments:badDimensions",
                "Inputs #1 and #%d must have compatible dimensions.", i + 1);
        }
    }
}

/*=========================================================================
 * Hanlde output arguments 
 *=======================================================================*/
void handleOutputs(int nlhs, mxArray *plhs[])
{
    // get a pointer to the first output
    plhs[0] = mxCreateDoubleMatrix(OUT_SEG_SIZE, qrsLen, mxREAL);
    outSeg1 = mxGetPr(plhs[0]);
    
    // get a pointer to the first output
    plhs[1] = mxCreateDoubleMatrix(OUT_SEG_SIZE, qrsLen, mxREAL);
    outSeg2 = mxGetPr(plhs[1]);
}

/*=========================================================================
 * The initialization routine 
 *=======================================================================*/
void init()
{
    // initialize beat index
    bi = 0;
    
    // initialize filter impulse response
    for (mwSize i = 0; i < OUT_SEG_SIZE; i++) {
        filterH[i] = 1.0;
    }
}

/*=========================================================================
 * The main routine 
 *=======================================================================*/
void doTheJob()
{
    double time;
    
    // start time counter
    tic();
    
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
    
    // make some initializations
    init();
    
    // do the actual processing
    doTheJob();
}
