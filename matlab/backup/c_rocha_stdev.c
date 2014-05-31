/*=========================================================================
 * c_rocha_stdev.c
 * 
 *  Title: computation of ST deviation based on isoelectric and J points
 *  Author:     Diego Sogari
 *  Modified:   16/May/2014
 *
 *  Intputs:
 *      1. List of Beats
 *      2. Location of isoelectric point
 *      3. Location of J point (measure 1)
 *      4. Location of J point (measure 2)
 *
 *  Outputs:
 *      1. measures of ST deviation (2xM matrix)
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
 *=======================================================================*/

/*=========================================================================
 * Constants
 *=======================================================================*/
#define MIN_INPUTS      4       // minimum number of input arguments
#define MAX_INPUTS      4       // maximum number of input arguments
#define MIN_OUTPUTS     1       // minimum number of output arguments
#define MAX_OUTPUTS     1       // maximum number of output arguments

/*=========================================================================
 * Input and output variables
 *=======================================================================*/
static double *beatList;        // list of beats
static mwSize qrsLen;           // number of beats
static double *isoList;         // list of isoelectric points
static double *jayList1;        // list of jay points (measure 1)
static double *jayList2;        // list of jay points (measure 2)
static double *outSTdev;        // output list of ST deviation

/*=========================================================================
 * Other variables
 *=======================================================================*/
static mwSize bi;               // current beat index
static mwSize frameSize;        // size of one beat

/*=========================================================================
 * Lookup for the beat list
 *=======================================================================*/
#define beat(I,J)   (beatList[(I)*frameSize+(J)])
#define stdev(I,J)  (outSTdev[(I)*2+(J)])

/*=========================================================================
 * Event triggered for a new beat
 *=======================================================================*/
void onNewBeat()
{
    mwSize iso = (mwSize)isoList[bi] - 1;
    mwSize jay1 = (mwSize)jayList1[bi] - 1;
    mwSize jay2 = (mwSize)jayList2[bi] - 1;
    
    stdev(bi,0) = beat(bi, jay1) - beat(bi, iso);
    stdev(bi,1) = beat(bi, jay2) - beat(bi, iso);
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
            "EcgToolbox:c_rocha_stdev:nrhs",
            "%d input(s) required and %d optional",
            MIN_INPUTS, MAX_INPUTS - MIN_INPUTS);
    }
    // check for proper number of output arguments
    if (nlhs < MIN_OUTPUTS || nlhs > MAX_OUTPUTS) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_rocha_stdev:nlhs",
            "%d output(s) required and %d optional",
            MIN_OUTPUTS, MAX_OUTPUTS - MIN_OUTPUTS);
    }
    // make sure all input arguments are of type double
    for (i = 0; i < nrhs; i++) {
        if (!mxIsDouble(prhs[i]) || mxIsComplex(prhs[i])) {
            mexErrMsgIdAndTxt(
                "EcgToolbox:c_rocha_stdev:notDouble",
                "Input #%d must be of type double.", i + 1);
        }
    }
    // make sure the first input argument has more than one line
    if (mxGetM(prhs[0]) == 1) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_rocha_stdev:badDimensions",
            "Input #1 must have more than one line.");
    }
    // make sure the input arguments 2-4 have exactly one column
    for (i = 1; i < nrhs; i++) {
        if (mxGetN(prhs[i]) != 1) {
            mexErrMsgIdAndTxt(
                "EcgToolbox:c_rocha_stdev:badDimensions",
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
    mwSize i;
    
    // get pointers to the data in the input vectors
    beatList = mxGetPr(prhs[0]);
    isoList = mxGetPr(prhs[1]);
    jayList1 = mxGetPr(prhs[2]);
    jayList2 = mxGetPr(prhs[3]);
    
    // get the dimensions of the first input
    *nrows = (mwSize)mxGetM(prhs[0]);
    *ncols = (mwSize)mxGetN(prhs[0]);
    
    // make sure the input arguments 2-4 have compatible dimensions
    for (i = 1; i < nrhs; i++) {
        if (mxGetM(prhs[i]) != *ncols) {
            mexErrMsgIdAndTxt(
                "EcgToolbox:c_rocha_stdev:badDimensions",
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
    plhs[0] = mxCreateDoubleMatrix(2, qrsLen, mxREAL);
    outSTdev = mxGetPr(plhs[0]);
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
