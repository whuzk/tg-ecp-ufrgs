/*=========================================================================
 * c_mohebbi_stdiff.c
 * 
 *  Title: computation of ST segment differences ofr the Mohebbi method
 *  Author:     Diego Sogari
 *  Modified:   16/May/2014
 *
 *  Intputs:
 *      1. Reference segment
 *      2. List of test segments
 *
 *  Outputs:
 *      1. List of ST segment differences
 *
 *  *all inputs and outputs required
 *
 *=======================================================================*/
#include <math.h>
#include "mex.h"
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
 *  IDX     index
 *=======================================================================*/

/*=========================================================================
 * Constants
 *=======================================================================*/
#define MIN_INPUTS      2       // minimum number of input arguments
#define MAX_INPUTS      2       // maximum number of input arguments
#define MIN_OUTPUTS     1       // minimum number of output arguments
#define MAX_OUTPUTS     1       // maximum number of output arguments

/*=========================================================================
 * Input and output variables
 *=======================================================================*/
static double *stsegRef;        // reference ST segment
static double *stsegList;       // list of test ST segments
static mwSize qrsLen;           // number of segments
static double *outSTdiff;       // output list of ST segment differences

/*=========================================================================
 * Other variables
 *=======================================================================*/
static mwSize bi;               // current beat index
static mwSize frameSize;        // size of one segment

/*=========================================================================
 * Lookup for the segment lists
 *=======================================================================*/
#define stseg(I,J)  (stsegList[(I)*frameSize+(J)])
#define stout(I,J)  (outSTdiff[(I)*frameSize+(J)])

/*=========================================================================
 * Event triggered for a new beat
 *=======================================================================*/
void onNewBeat()
{
    double_subtract(&stseg(bi, 0), stsegRef, &stout(bi, 0), frameSize);
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
            "EcgToolbox:c_mohebbi_stdiff:nrhs",
            "%d input(s) required and %d optional",
            MIN_INPUTS, MAX_INPUTS - MIN_INPUTS);
    }
    // check for proper number of output arguments
    if (nlhs < MIN_OUTPUTS || nlhs > MAX_OUTPUTS) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_mohebbi_stdiff:nlhs",
            "%d output(s) required and %d optional",
            MIN_OUTPUTS, MAX_OUTPUTS - MIN_OUTPUTS);
    }
    // make sure all input arguments are of type double
    for (i = 0; i < nrhs; i++) {
        if (!mxIsDouble(prhs[i]) || mxIsComplex(prhs[i])) {
            mexErrMsgIdAndTxt(
                "EcgToolbox:c_mohebbi_stdiff:notDouble",
                "Input #%d must be of type double.", i + 1);
        }
    }
    // make sure the first input argument has exactly one column
    if (mxGetN(prhs[0]) != 1) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_mohebbi_stdiff:badDimensions",
            "Input #1 must have exactly one column.");
    }
    // make sure the second input argument has more than one line
    if (mxGetM(prhs[1]) == 1) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_mohebbi_stdiff:badDimensions",
            "Input #2 must have more than one line.");
    }
}

/*=========================================================================
 * Hanlde input arguments 
 *=======================================================================*/
void handleInputs( int nrhs, const mxArray *prhs[],
                   mwSize *nrows, mwSize *ncols)
{
    // get pointers to the data in the input vectors
    stsegRef = mxGetPr(prhs[0]);
    stsegList = mxGetPr(prhs[1]);
    
    // get the dimensions of the second input
    *nrows = (mwSize)mxGetM(prhs[1]);
    *ncols = (mwSize)mxGetN(prhs[1]);
    
    // make sure the input arguments have compatible dimensions
    if (mxGetM(prhs[0]) != *nrows) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_mohebbi_stdiff:badDimensions",
            "Inputs #1 and #2 must have the same number of lines.");
    }
}

/*=========================================================================
 * Hanlde output arguments 
 *=======================================================================*/
void handleOutputs(int nlhs, mxArray *plhs[])
{
    // get a pointer to the first output
    plhs[0] = mxCreateDoubleMatrix(frameSize, qrsLen, mxREAL);
    outSTdiff = mxGetPr(plhs[0]);
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
    
    // get the frame size
    frameSize = nrows;
    
    // get the number of beats in the list
    qrsLen = ncols;
    
    // handle output arguments
    handleOutputs(nlhs, plhs);
    
    // do the actual processing
    doTheJob();
}
