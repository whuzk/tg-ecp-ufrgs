/*=========================================================================
 * c_mohebbi_features_inst.c
 * 
 *  Title: features extraction according to Mohebbi and Moghadam
 *  Author:     Diego Sogari
 *  Modified:   16/May/2014
 *
 *  Intputs:
 *      1. List of Beats
 *      2. Reference beat
 *      3. Default value of J points
 *      4. Default value of template J point
 *      5. Sampling frequency (in Hertz)
 *      6. Gain of the input signal
 *
 *  Outputs:
 *      1. List of ST segment differences
 *      2. List of ST segments
 *      3. processing time 1*
 *      4. processing time 2*
 *      5. processing time 3*
 *      6. processing time 4*
 *
 *  *all inputs and outputs required
 *
 *=======================================================================*/
#include <string.h>
#include <math.h>
#include "mex.h"
#include "c_upfirdn.h"
#include "c_intfilter.h"
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
#define MIN_INPUTS      6       // minimum number of input arguments
#define MAX_INPUTS      6       // maximum number of input arguments
#define MIN_OUTPUTS     6       // minimum number of output arguments
#define MAX_OUTPUTS     6       // maximum number of output arguments

#define MAF_WIDTH       0.05    // width of the moving-average (in seconds)
#define DER_ORDERN      1       // derivative filter order N
#define DER_ORDERM      0       // derivative filter order M
#define SEGMENT_SIZE    20      // size of the segments (in samples)
#define HALF_SEG_SIZE   0.08    // half the size of a segment (in seconds)

/*=========================================================================
 * Input and output variables
 *=======================================================================*/
static double *beatList;        // list of beats
static double *refBeat;         // reference beat
static double *defJay;          // list of default J points
static double defTempJay;       // default template J point
static mwSize qrsLen;           // number of beats
static double sampFreq;         // sampling frequency
static double inputGain;        // gain of the input signal
static double *outSTdiff;       // output list of ST segment differences
static double *outSeg;          // output list of ST segments
static double *procTime1;       // processing time 1
static double *procTime2;       // processing time 2
static double *procTime3;       // processing time 3
static double *procTime4;       // processing time 4

/*=========================================================================
 * Other variables
 *=======================================================================*/
static mwSize bi;               // current beat index
static mwSize frameSize;        // size of one beat
static mwSize stsegSize;        // size of a segment
static mwSize rPeak;            // location of R peak
static intfobject maFilter;     // moving-average filter object
static intfobject deFilter;     // derivative filter object
static int *beatBuf;            // buffer for the current beat
static double *segRefBuf;       // buffer for the reference segment
static double *segBuf;          // buffer for the segment
static mwSize jayPoint;         // jay point
static int thresh;              // threshold for detection
static mwSize delay;            // sum of detection filter delays
static mwSize searchL1;         // search limit 1
static mwSize searchL2;         // search limit 2
static double *windowBuf;       // window buffer
static double *filterBuf;       // filter buffer
static mwSize filterLen;        // length of filter buffer
static double *resampBuf;       // resampling buffer
static mwSize resampLen;        // length of resampling buffer

/*=========================================================================
 * Instrumentation variables
 *=======================================================================*/
long long totStartClock;
long long jayStartClock;
long long segStartClock;
long long difStartClock;

/*=========================================================================
 * Lookup for the beat list
 *=======================================================================*/
#define beat(I,J)       (beatList[(I)*frameSize+(J)])
#define stdiff(I,J)     (outSTdiff[(I)*SEGMENT_SIZE+(J)])
#define stseg(I,J)      (outSeg[(I)*SEGMENT_SIZE+(J)])

/*=========================================================================
 * Filter the beat signal
 *=======================================================================*/
void filterBeat(double *x)
{
    int sample;
    mwSize i;
    
    for (i = 0; i < delay; i++) {
        sample = intfnewx(&maFilter, i, (int)x[i]);
        intfnewx(&deFilter, i, sample);
    }
    for (i = delay; i < frameSize; i++) {
        sample = intfnewx(&maFilter, i, (int)x[i]);
        beatBuf[i - delay] = intfnewx(&deFilter, i, sample);
    }
    for (i = frameSize; i < frameSize + delay; i++) {
        sample = intfnewx(&maFilter, i, 0);
        beatBuf[i - delay] = intfnewx(&deFilter, i, sample);
    }
}

/*=========================================================================
 * Search a point where the signal is below the threshold
 *=======================================================================*/
mwSize get_point(mwSize istart, mwSize iend, mwSize len, mwSize def)
{
    mwSize i, m, pos, inc, firstAbove;
    
    if (istart <= iend) {
        inc = 1;
        m = len >> 1;
    }
    else {
        inc = -1;
        m = -(len >> 1);
    }
    
    pos = def;
    firstAbove = 0;
    for (i = istart; i != iend + inc; i += inc) {
        if (abs(beatBuf[i]) > thresh) {
            firstAbove = 1;
        }
        else if (firstAbove == len) {
            pos = i - m;
            break;
        }
        else firstAbove++;
    }
    
    return pos;
}

/*=========================================================================
 * detect J point
 *=======================================================================*/
void detectJpoint(mwSize defJ)
{
    mwSize istart = rPeak + searchL1;
    mwSize iend = rPeak + searchL2;
    jayPoint = get_point(istart, iend, searchL1, defJ);
}

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
 * Extract ST segment
 *=======================================================================*/
void exctractSegment(double *x, double *y)
{
    resamp(segBuf, x, stsegSize, SEGMENT_SIZE, stsegSize);
    memcpy(y, segBuf, SEGMENT_SIZE * sizeof(double));
}

/*=========================================================================
 * Calculate ST segment difference
 *=======================================================================*/
void stsegmentDiff()
{
    double_subtract(segBuf, segRefBuf, &stdiff(bi, 0), SEGMENT_SIZE);
}

/*=========================================================================
 * Event triggered for a new beat
 *=======================================================================*/
void onNewBeat()
{
    // detect J point
    jayStartClock = tic();
    filterBeat(&beat(bi, 0));
    detectJpoint((int)defJay[bi] - 1);
    procTime2[bi] = toc(jayStartClock);
    
    // extract segment
    segStartClock = tic();
    exctractSegment(&beat(bi, jayPoint), &stseg(bi, 0));
    procTime3[bi] = toc(segStartClock);
    
    // calculate difference
    difStartClock = tic();
    stsegmentDiff();
    procTime4[bi] = toc(difStartClock);
    
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
            "EcgToolbox:c_mohebbi_features:nrhs",
            "%d input(s) required and %d optional",
            MIN_INPUTS, MAX_INPUTS - MIN_INPUTS);
    }
    // check for proper number of output arguments
    if (nlhs < MIN_OUTPUTS || nlhs > MAX_OUTPUTS) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_mohebbi_features:nlhs",
            "%d output(s) required and %d optional",
            MIN_OUTPUTS, MAX_OUTPUTS - MIN_OUTPUTS);
    }
    // make sure all input arguments are of type double
    for (i = 0; i < nrhs; i++) {
        if (!mxIsDouble(prhs[i]) || mxIsComplex(prhs[i])) {
            mexErrMsgIdAndTxt(
                "EcgToolbox:c_mohebbi_features:notDouble",
                "Input #%d must be of type double.", i + 1);
        }
    }
    // make sure the first input argument has more than one line
    if (mxGetM(prhs[0]) == 1) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_mohebbi_features:badDimensions",
            "Input #1 must have more than one line.");
    }
    // make sure the input arguments 2-3 are vectors
    for (i = 1; i < 3; i++) {
        if (mxGetM(prhs[i]) != 1 && mxGetN(prhs[i]) != 1) {
            mexErrMsgIdAndTxt(
                "EcgToolbox:c_mohebbi_features:notVector",
                "Input #%d must be a vector.", i + 1);
        }
    }
    // make sure the remaining input arguments are all scalars
    for (i = 3; i < nrhs; i++) {
        if (mxGetNumberOfElements(prhs[i]) != 1) {
            mexErrMsgIdAndTxt(
                "EcgToolbox:c_mohebbi_features:notScalar",
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
    refBeat = mxGetPr(prhs[1]);
    defJay = mxGetPr(prhs[2]);
    
    // get the dimensions of the first input
    *nrows = (mwSize)mxGetM(prhs[0]);
    *ncols = (mwSize)mxGetN(prhs[0]);
    
    // make sure the input argument 2 has compatible dimensions
    if (max(mxGetM(prhs[1]), mxGetN(prhs[1])) != *nrows) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_mohebbi_features:badDimensions",
            "Input #1 and #2 must have compatible dimensions.");
    }
    // make sure the input argument 3 has compatible dimensions
    if (max(mxGetM(prhs[2]), mxGetN(prhs[2])) != *ncols) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_mohebbi_features:badDimensions",
            "Inputs #1 and #3 must have compatible dimensions.");
    }
    
    // get the default template J point
    defTempJay = mxGetScalar(prhs[3]);
    
    // get the sampling frequency
    sampFreq = mxGetScalar(prhs[4]);
    
    // get the gain of the input
    inputGain = mxGetScalar(prhs[5]);
}

/*=========================================================================
 * Hanlde output arguments 
 *=======================================================================*/
void handleOutputs(int nlhs, mxArray *plhs[])
{
    // get a pointer to the first output
    plhs[0] = mxCreateDoubleMatrix(SEGMENT_SIZE, qrsLen, mxREAL);
    outSTdiff = mxGetPr(plhs[0]);
    
    // get a pointer to the second output
    plhs[1] = mxCreateDoubleMatrix(SEGMENT_SIZE, qrsLen, mxREAL);
    outSeg = mxGetPr(plhs[1]);
    
    // get a pointer to the instrumentation outputs
    plhs[2] = mxCreateDoubleMatrix(qrsLen, 1, mxREAL);
    procTime1 = mxGetPr(plhs[2]);
    plhs[3] = mxCreateDoubleMatrix(qrsLen, 1, mxREAL);
    procTime2 = mxGetPr(plhs[3]);
    plhs[4] = mxCreateDoubleMatrix(qrsLen, 1, mxREAL);
    procTime3 = mxGetPr(plhs[4]);
    plhs[5] = mxCreateDoubleMatrix(qrsLen, 1, mxREAL);
    procTime4 = mxGetPr(plhs[5]);
}

/*=========================================================================
 * The initialization routine 
 *=======================================================================*/
void init()
{
    mwSize n;
    
    // initialize beat index
    bi = 0;
    
    // initialize detection filter objeects
    initintfilter(maFilter);
    initintfilter(deFilter);
    
    // design detection filters
    design_maverage(&maFilter, sampFreq, MAF_WIDTH);
    design_derivative(&deFilter, DER_ORDERN, DER_ORDERM);
    
    // initialize numeric data
    stsegSize = (int)(HALF_SEG_SIZE * sampFreq) * 2;
    thresh = (int)(maFilter.gain * inputGain * 1.25 / sampFreq);
    delay = (int)ceil(maFilter.delay + deFilter.delay);
    searchL1 = (int)(0.02 * sampFreq);
    searchL2 = (int)(0.12 * sampFreq);
    
    // initialize buffers
    beatBuf = (int *)mxMalloc(frameSize * sizeof(int));
    segRefBuf = (double *)mxMalloc(SEGMENT_SIZE * sizeof(double));
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
    mwSize i;
    
    // initialize beat index
    bi = 0;
    
    // handle the template
    filterBeat(refBeat);
    detectJpoint((int)defTempJay - 1);
    exctractSegment(&refBeat[jayPoint], segRefBuf);
    
    // process one input sample at a time
    for (i = 0; i < qrsLen; i++) {
        totStartClock = tic();
        onNewBeat();
        procTime1[i] = toc(totStartClock);
    }
}

/*=========================================================================
 * The finalization routine 
 *=======================================================================*/
void finalize()
{
    // deallocate memory of the filter objects
    endintfilter(maFilter);
    endintfilter(deFilter);
    
    // deallocate memory
    mxFree(beatBuf);
    mxFree(segRefBuf);
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
    rPeak = frameSize >> 1;
    
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
