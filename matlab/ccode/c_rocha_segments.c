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
 *      5. Sampling frequency (in Hertz)
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
#define MIN_INPUTS      5       // minimum number of input arguments
#define MAX_INPUTS      5       // maximum number of input arguments
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
static double sampFreq;         // sampling frequency
static double *outSeg1;         // output list of segments
static double *outSeg2;         // output list of segments

/*=========================================================================
 * Other variables
 *=======================================================================*/
static mwSize bi;               // current beat index
static mwSize frameSize;        // size of one beat
static double *windowBuf;       // window buffer
static double *filterBuf;       // filter buffer
static mwSize filterLen;        // length of filter buffer
static double *resampBuf;       // resampling buffer
static mwSize resampLen;        // length of resampling buffer

/*=========================================================================
 * Lookup for vectors
 *=======================================================================*/
#define beat(I,J)   (beatList[(I)*frameSize+(J)])
#define seg1(I,J)   (outSeg1[(I)*OUT_SEG_SIZE+(J)])
#define seg2(I,J)   (outSeg2[(I)*OUT_SEG_SIZE+(J)])

/*=========================================================================
 * Resample
 *=======================================================================*/
void resamp(double *y, double *x, mwSize Lx, int p, int q)
{
    double wc, win, c;
    mwSize Ly, delay;
    int pqmax, M2;
    
    // compute parameters
    rational(p / (double)q, &p, &q, 1.0e-5);
    pqmax = max(p, q);
    wc = PI / (double)pqmax;
    M2 = (filterLen - 1) >> 1;
    
    // compute filter impulse response
    for (mwSize n = 0; n < filterLen; n++) {
        if (n != M2) {
            c = n - M2;
            filterBuf[n] = p * (sin(wc * c) / (PI * c)) * windowBuf[n];
        }
        else filterBuf[n] = p * (wc / PI) * windowBuf[n];
    }
    
    // perform the resampling
    Ly = (int)ceil(((Lx - 1) * p + filterLen) / (double)q);
    memset(resampBuf, 0, filterLen * sizeof(double));
    if (Ly > resampLen) {
        mexPrintf("Warning: Ly > resampLen\n");
        Ly = resampLen;
    }
    upfirdn(resampBuf, Ly, x, Lx, filterBuf, filterLen, p, q);
    
    // remove delay and write to output
    delay = (int)floor(M2 / (double)q);
    Ly = (int)ceil(Lx * p / (double)q);
    memcpy(y, &resampBuf[delay], Ly * sizeof(double));
}

/*=========================================================================
 * Event triggered for a new beat
 *=======================================================================*/
void onNewBeat()
{
    mwSize istart, iend, len;
    double *x, *y;
    
    istart = (mwSize)begList[bi] - 1;
    iend = (mwSize)jayList[bi] - 1;
    if (istart < iend) {
        len = iend - istart + 1;
        x = &beat(bi, istart);
        y = &seg1(bi, 0);
        resamp(y, x, len, OUT_SEG_SIZE, len);
    }
    
    istart = (mwSize)jayList[bi] - 1;
    iend = (mwSize)endList[bi] - 1;
    if (istart < iend) {
        len = iend - istart + 1;
        x = &beat(bi, istart);
        y = &seg2(bi, 0);
        resamp(y, x, len, OUT_SEG_SIZE, len);
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
    for (mwSize i = 1; i < 4; i++) {
        if (mxGetN(prhs[i]) != 1) {
            mexErrMsgIdAndTxt(
                "EcgToolbox:c_rocha_segments:badDimensions",
                "Input #%d must have exactly one column.", i + 1);
        }
    }
    // make sure the remaining input arguments are all scalars
    for (mwSize i = 4; i < nrhs; i++) {
        if (mxGetNumberOfElements(prhs[i]) != 1) {
            mexErrMsgIdAndTxt(
                "EcgToolbox:c_mohebbi_segments:notScalar",
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
    beatList = mxGetPr(prhs[0]);
    begList = mxGetPr(prhs[1]);
    jayList = mxGetPr(prhs[2]);
    endList = mxGetPr(prhs[3]);
    
    // get the dimensions of the first input
    *nrows = (mwSize)mxGetM(prhs[0]);
    *ncols = (mwSize)mxGetN(prhs[0]);
    
    // make sure the input arguments 2-4 have compatible dimensions
    for (mwSize i = 1; i < 4; i++) {
        if (mxGetM(prhs[i]) != *ncols) {
            mexErrMsgIdAndTxt(
                "EcgToolbox:c_rocha_segments:badDimensions",
                "Inputs #1 and #%d must have compatible dimensions.", i + 1);
        }
    }
    
    // get the sampling frequency
    sampFreq = mxGetScalar(prhs[4]);
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
    
    // initialize filter and resampling buffers
    filterLen = (int)round(2 * sampFreq) + 1;
    filterBuf = (double *)mxMalloc(filterLen * sizeof(double));
    resampLen = OUT_SEG_SIZE + filterLen - 1;
    resampBuf = (double *)mxMalloc(resampLen * sizeof(double));
    
    // initialize window buffer
    windowBuf = (double *)mxMalloc(filterLen * sizeof(double));
    for (mwSize n = 0; n < filterLen; n++) {
        windowBuf[n] = 0.54 - 0.46 * cos(2*PI*n / (double)(filterLen - 1));
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
 * The finalization routine 
 *=======================================================================*/
void finalize()
{
    // deallocate memory
    mxFree(filterBuf);
    mxFree(resampBuf);
    mxFree(windowBuf);
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
    
    // perform final adjustments
    finalize();
}
