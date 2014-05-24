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
 *      3. Sampling frequency (in Hertz)
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
#define MIN_INPUTS      3       // minimum number of input arguments
#define MAX_INPUTS      3       // maximum number of input arguments
#define MIN_OUTPUTS     1       // minimum number of output arguments
#define MAX_OUTPUTS     1       // maximum number of output arguments

#define OUT_SEG_SIZE    250     // size of the output segment (in samples)

/*=========================================================================
 * Input and output variables
 *=======================================================================*/
static double *beatList;        // list of beats
static double *rrList;          // list of RR intervals
static mwSize qrsLen;           // number of beats
static double sampFreq;         // sampling frequency
static double *outSeg;          // output list of segments

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
#define seg(I,J)    (outSeg[(I)*OUT_SEG_SIZE+(J)])

/*=========================================================================
 * Resample
 *=======================================================================*/
void resamp(double *y, double *x, mwSize Lx, int p, int q)
{
    mwSize n, Ly, delay;
    int pqmax, M2;
    double wc, c;
    
    // compute parameters
    rational(p / (double)q, &p, &q, 1.0e-5);
    pqmax = max(p, q);
    wc = PI / (double)pqmax;
    M2 = (filterLen - 1) >> 1;
    
    // compute filter impulse response
    for (n = 0; n < filterLen; n++) {
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
    mwSize rr, istart;
    double *x, *y;
    
    rr = min(frameSize, (mwSize)rrList[bi]);
    if (rr % 2 != 0) {
        rr--;
    }
    istart = (frameSize - rr) >> 1;
    x = &beat(bi, istart);
    y = &seg(bi, 0);
    resamp(y, x, rr, OUT_SEG_SIZE, rr);
    
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
    for (i = 0; i < nrhs; i++) {
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
    // make sure the remaining input arguments are all scalars
    for (i = 2; i < nrhs; i++) {
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
    
    // get the sampling frequency
    sampFreq = mxGetScalar(prhs[2]);
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
 * The initialization routine 
 *=======================================================================*/
void init()
{
    mwSize n;
    
    // initialize beat index
    bi = 0;
    
    // initialize filter and resampling buffers
    filterLen = (int)(2 * sampFreq) + 1;
    filterBuf = (double *)mxMalloc(filterLen * sizeof(double));
    resampLen = OUT_SEG_SIZE + filterLen - 1;
    resampBuf = (double *)mxMalloc(resampLen * sizeof(double));
    
    // initialize window buffer
    windowBuf = (double *)mxMalloc(filterLen * sizeof(double));
    for (n = 0; n < filterLen; n++) {
        windowBuf[n] = 0.54 - 0.46 * cos(2*PI*n / (double)(filterLen - 1));
    }
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
