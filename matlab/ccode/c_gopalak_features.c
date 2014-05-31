/*=========================================================================
 * c_gopalak_features.c
 * 
 *  Title: features extraction according to Gopalakrishnan et al.
 *  Author:     Diego Sogari
 *  Modified:   16/May/2014
 *
 *  Intputs:
 *      1. List of Beats
 *      2. List of RR intervals
 *      3. Sampling frequency (in Hertz)
 *
 *  Outputs:
 *      1. List of hermite coefficients
 *      2. List of segments
 *
 *  *all inputs and outputs required
 *
 *=======================================================================*/
#include <string.h>
#include <math.h>
#include "mex.h"
#include "c_upfirdn.h"
#include "c_hermcoeff.h"
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
#define MIN_INPUTS      3       // minimum number of input arguments
#define MAX_INPUTS      3       // maximum number of input arguments
#define MIN_OUTPUTS     1       // minimum number of output arguments
#define MAX_OUTPUTS     2       // maximum number of output arguments

#define NUM_HERM_COEF   50      // number of hermite coefficients
#define SEGMENT_SIZE    250     // size of the segments (in samples)

/*=========================================================================
 * Input and output variables
 *=======================================================================*/
static double *beatList;        // list of beats
static double *rrList;          // list of RR intervals
static mwSize qrsLen;           // number of beats
static double sampFreq;         // sampling frequency
static double *outSeg;         // output list of segments
static double *outCoeff;       // output list of hermite coefficients

/*=========================================================================
 * Other variables
 *=======================================================================*/
static mwSize bi;               // current beat index
static mwSize frameSize;        // size of one beat
static double *segBuf;          // buffer for the segment
static double *windowBuf;       // window buffer
static double *filterBuf;       // filter buffer
static mwSize filterLen;        // length of filter buffer
static double *resampBuf;       // resampling buffer
static mwSize resampLen;        // length of resampling buffer

/*=========================================================================
 * Lookup for the beat list
 *=======================================================================*/
#define beat(I,J)   (beatList[(I)*frameSize+(J)])
#define coef(I,J)   (outCoeff[(I)*NUM_HERM_COEF+(J)])
#define seg(I,J)    (outSeg[(I)*SEGMENT_SIZE+(J)])

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
    M2 = (int)(filterLen - 1) >> 1;
    
    // compute filter impulse response
    for (n = 0; n < filterLen; n++) {
        if (n != M2) {
            c = n - (double)M2;
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
 * Extract beat segment
 *=======================================================================*/
void exctractSegment()
{
    mwSize rr, istart;
    
    rr = min(frameSize, (mwSize)rrList[bi]);
    if (rr % 2 != 0) {
        rr--;
    }
    istart = (frameSize - rr) >> 1;
    resamp(segBuf, &beat(bi, istart), rr, SEGMENT_SIZE, rr);
    if (outSeg != NULL) {
        memcpy(&seg(bi, 0), segBuf, SEGMENT_SIZE * sizeof(double));
    }
}

/*=========================================================================
 * Calculate Hermite expansion
 *=======================================================================*/
void hermiteExpansion()
{
    size_t m = NUM_HERM_COEF;
    size_t p = SEGMENT_SIZE;
    size_t n = 1;
    char *chn = "N";
    double one = 1.0;
    double zero = 0.0;
    double *A = gopalak_hcoef;
    double *C = &coef(bi, 0);
    
    dgemm(chn, chn, &m, &n, &p, &one, A, &m, segBuf, &p, &zero, C, &m);
}

/*=========================================================================
 * Event triggered for a new beat
 *=======================================================================*/
void onNewBeat()
{
    // Hermite coefficients
    exctractSegment();
    hermiteExpansion();
    
    // increment beat index
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
            "EcgToolbox:c_gopalak_features:nrhs",
            "%d input(s) required and %d optional",
            MIN_INPUTS, MAX_INPUTS - MIN_INPUTS);
    }
    // check for proper number of output arguments
    if (nlhs < MIN_OUTPUTS || nlhs > MAX_OUTPUTS) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_gopalak_features:nlhs",
            "%d output(s) required and %d optional",
            MIN_OUTPUTS, MAX_OUTPUTS - MIN_OUTPUTS);
    }
    // make sure all input arguments are of type double
    for (i = 0; i < nrhs; i++) {
        if (!mxIsDouble(prhs[i]) || mxIsComplex(prhs[i])) {
            mexErrMsgIdAndTxt(
                "EcgToolbox:c_gopalak_features:notDouble",
                "Input #%d must be of type double.", i + 1);
        }
    }
    // make sure the first input argument has more than one line
    if (mxGetM(prhs[0]) == 1) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_gopalak_features:badDimensions",
            "Input #1 must have more than one line.");
    }
    // make sure the input arguments 2 is a vector
    if (mxGetM(prhs[1]) != 1 && mxGetN(prhs[1]) != 1) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_gopalak_features:notVector",
            "Input #2 must be a vector.");
    }
    // make sure the remaining input arguments are all scalars
    for (i = 2; i < nrhs; i++) {
        if (mxGetNumberOfElements(prhs[i]) != 1) {
            mexErrMsgIdAndTxt(
                "EcgToolbox:c_gopalak_features:notScalar",
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
    
    // make sure the input argument 2 has compatible dimensions
    if (max(mxGetM(prhs[1]), mxGetN(prhs[1])) != *ncols) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_gopalak_features:badDimensions",
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
    plhs[0] = mxCreateDoubleMatrix(NUM_HERM_COEF, qrsLen, mxREAL);
    outCoeff = mxGetPr(plhs[0]);
    
    // optional output
    if (nlhs > 1) {
        plhs[1] = mxCreateDoubleMatrix(SEGMENT_SIZE, qrsLen, mxREAL);
        outSeg = mxGetPr(plhs[1]);
    }
    else outSeg = NULL;
}

/*=========================================================================
 * The initialization routine 
 *=======================================================================*/
void init()
{
    mwSize n;
    
    // initialize beat index
    bi = 0;
    
    // initialize buffers
    segBuf = (double *)mxMalloc(SEGMENT_SIZE * sizeof(double));
    filterLen = (int)(2 * sampFreq) + 1;
    filterBuf = (double *)mxMalloc(filterLen * sizeof(double));
    resampLen = SEGMENT_SIZE + filterLen - 1;
    resampBuf = (double *)mxMalloc(resampLen * sizeof(double));
    
    // initialize resampling window buffer
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
 * The finalization routine 
 *=======================================================================*/
void finalize()
{
    // deallocate memory
    mxFree(segBuf);
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
