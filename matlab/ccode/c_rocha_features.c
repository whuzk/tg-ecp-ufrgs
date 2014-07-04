/*=========================================================================
 * c_rocha_features.c
 * 
 *  Title: features extraction according to Rocha et al.
 *  Author:     Diego Sogari
 *  Modified:   16/May/2014
 *
 *  Intputs:
 *      1. List of Beats
 *      2. List of RR intervals
 *      3. Default value of isoelectric points
 *      4. Default value of J points
 *      5. List of beat ending points
 *      6. Sampling frequency (in Hertz)
 *      7. Gain of the input signal
 *
 *  Outputs:
 *      1. List of ST deviation measures
 *      2. List of hermite coefficients 1
 *      3. List of hermite coefficients 2
 *      4. List of segments 1
 *      5. List of segments 2
 *
 *  *all inputs and outputs required
 *
 *=======================================================================*/
#include <string.h>
#include <math.h>
#include "mex.h"
#include "blas.h"
#include "c_upfirdn.h"
#include "c_hermcoeff.h"
#include "c_intfilter.h"
#include "c_mathutils.h"

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
#define MIN_INPUTS      7       // minimum number of input arguments
#define MAX_INPUTS      7       // maximum number of input arguments
#define MIN_OUTPUTS     3       // minimum number of output arguments
#define MAX_OUTPUTS     5       // maximum number of output arguments

#define NUM_HR_LIMITS   3       // number of heart rate limits
#define NUM_MS_POINTS   3+1     // number of measuring points
#define MAF_WIDTH       0.05    // width of the moving-average (in seconds)
#define DER_ORDERN      1       // derivative filter order N
#define DER_ORDERM      0       // derivative filter order M
#define NUM_HERM_COEF   6       // number of hermite coefficients
#define SEGMENT_SIZE    64      // size of the segments (in samples)

/*=========================================================================
 * Input and output variables
 *=======================================================================*/
static double *beatList;        // list of beats
static double *rrList;          // list of RR intervals
static double *defIso;          // list of default isoelectric points
static double *defJay;          // list of default J points
static double *endList;         // list of beat ending points
static mwSize qrsLen;           // number of beats
static double sampFreq;         // sampling frequency
static double inputGain;        // gain of the input signal
static double *outSTdev;        // output list of ST deviation
static double *outSeg1;         // output list of segments 1
static double *outSeg2;         // output list of segments 2
static double *outCoeff1;       // output list of hermite coefficients 1
static double *outCoeff2;       // output list of hermite coefficients 2

/*=========================================================================
 * Other variables
 *=======================================================================*/
static mwSize bi;               // current beat index
static mwSize frameSize;        // size of one beat
static mwSize rPeak;            // location of R peak
static double hrLimits[] = {100.0, 110.0, 120.0};       // heart rate limits
static double msPoints[] = {0.12, 0.112, 0.104, 0.1};   // measuring points
static intfobject maFilter;     // moving-average filter object
static intfobject deFilter;     // derivative filter object
static int *beatBuf;            // buffer for the current beat
static double *seg1Buf;         // buffer for the segment 1
static double *seg2Buf;         // buffer for the segment 2
static mwSize isoPoint;         // isoelectric point
static mwSize jayPoint1;        // jay point 1
static mwSize jayPoint2;        // jay point 2
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
 * Lookup for the beat list
 *=======================================================================*/
#define beat(I,J)   (beatList[(I)*frameSize+(J)])
#define stdev(I,J)  (outSTdev[(I)*2+(J)])
#define coef1(I,J)  (outCoeff1[(I)*NUM_HERM_COEF+(J)])
#define coef2(I,J)  (outCoeff2[(I)*NUM_HERM_COEF+(J)])
#define seg1(I,J)   (outSeg1[(I)*SEGMENT_SIZE+(J)])
#define seg2(I,J)   (outSeg2[(I)*SEGMENT_SIZE+(J)])

/*=========================================================================
 * Filter the beat signal
 *=======================================================================*/
void filterBeat()
{
    int sample;
    mwSize i;
    
    for (i = 0; i < delay; i++) {
        sample = intfnewx(&maFilter, i, (int)beat(bi, i));
        intfnewx(&deFilter, i, sample);
    }
    for (i = delay; i < frameSize; i++) {
        sample = intfnewx(&maFilter, i, (int)beat(bi, i));
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
 * detect I and J points
 *=======================================================================*/
void detectIJpoints()
{
    mwSize i, istart, iend;
    double heartRate;
    
    // pang J point
    heartRate = 60.0 * sampFreq / rrList[bi];
    for (i = 0; i < NUM_HR_LIMITS && heartRate >= hrLimits[i]; i++);
    jayPoint1 = rPeak + (int)(msPoints[i] * sampFreq);
    
    // isoelectric point
    istart = rPeak - searchL1;
    iend = rPeak - searchL2;
    isoPoint = get_point(istart, iend, searchL1, (int)defIso[bi] - 1);
    
    // J point
    istart = rPeak + searchL1;
    iend = rPeak + searchL2;
    jayPoint2 = get_point(istart, iend, searchL1, (int)defJay[bi] - 1);
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
    pqmax = MAX(p, q);
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
    if (Ly > resampLen) {
        mexPrintf("Warning: Ly > resampLen\n");
        Ly = resampLen;
    }
    memset(resampBuf, 0, Ly * sizeof(double));
    upfirdn(resampBuf, Ly, x, Lx, filterBuf, filterLen, p, q);
    
    // remove delay and write to output
    delay = (int)floor(M2 / (double)q);
    Ly = (int)ceil(Lx * p / (double)q);
    memcpy(y, &resampBuf[delay], Ly * sizeof(double));
}

/*=========================================================================
 * Extract beat segments
 *=======================================================================*/
void exctractSegments()
{
    mwSize istart, iend, len;
    
    istart = isoPoint;
    iend = jayPoint2;
    if (istart < iend) {
        len = iend - istart + 1;
        resamp(seg1Buf, &beat(bi, istart), len, SEGMENT_SIZE, (int)len);
        if (outSeg1 != NULL) {
            memcpy(&seg1(bi, 0), seg1Buf, SEGMENT_SIZE * sizeof(double));
        }
    }
    
    istart = jayPoint2;
    iend = (int)endList[bi] - 1;
    if (istart < iend) {
        len = iend - istart + 1;
        resamp(seg2Buf, &beat(bi, istart), len, SEGMENT_SIZE, (int)len);
        if (outSeg2 != NULL) {
            memcpy(&seg2(bi, 0), seg2Buf, SEGMENT_SIZE * sizeof(double));
        }
    }
}

/*=========================================================================
 * Calculate Hermite expansion
 *=======================================================================*/
void hermiteExpansion()
{
    ptrdiff_t m = NUM_HERM_COEF;
    ptrdiff_t p = SEGMENT_SIZE;
    ptrdiff_t n = 1;
    char *chn = "N";
    double one = 1.0;
    double zero = 0.0;
    double *A, *C;
    
    A = rocha_S1_hcoef;
    C = &coef1(bi, 0);
    dgemm(chn, chn, &m, &n, &p, &one, A, &m, seg1Buf, &p, &zero, C, &m);
    
    A = rocha_S2_hcoef;
    C = &coef2(bi, 0);
    dgemm(chn, chn, &m, &n, &p, &one, A, &m, seg2Buf, &p, &zero, C, &m);
}

/*=========================================================================
 * Event triggered for a new beat
 *=======================================================================*/
void onNewBeat()
{
    // ST deviation
    filterBeat();
    detectIJpoints();
    stdev(bi,0) = beat(bi, jayPoint1) - beat(bi, isoPoint);
    stdev(bi,1) = beat(bi, jayPoint2) - beat(bi, isoPoint);
    
    // Hermite coefficients
    exctractSegments();
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
            "EcgToolbox:c_rocha_features:nrhs",
            "%d input(s) required and %d optional",
            MIN_INPUTS, MAX_INPUTS - MIN_INPUTS);
    }
    // check for proper number of output arguments
    if (nlhs < MIN_OUTPUTS || nlhs > MAX_OUTPUTS) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_rocha_features:nlhs",
            "%d output(s) required and %d optional",
            MIN_OUTPUTS, MAX_OUTPUTS - MIN_OUTPUTS);
    }
    // make sure all input arguments are of type double
    for (i = 0; i < nrhs; i++) {
        if (!mxIsDouble(prhs[i]) || mxIsComplex(prhs[i])) {
            mexErrMsgIdAndTxt(
                "EcgToolbox:c_rocha_features:notDouble",
                "Input #%d must be of type double.", i + 1);
        }
    }
    // make sure the first input argument has more than one line
    if (mxGetM(prhs[0]) == 1) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_rocha_features:badDimensions",
            "Input #1 must have more than one line.");
    }
    // make sure the input arguments 2-5 are vectors
    for (i = 1; i < 5; i++) {
        if (mxGetM(prhs[i]) != 1 && mxGetN(prhs[i]) != 1) {
            mexErrMsgIdAndTxt(
                "EcgToolbox:c_rocha_features:notVector",
                "Input #%d must be a vector.", i + 1);
        }
    }
    // make sure the remaining input arguments are all scalars
    for (i = 5; i < nrhs; i++) {
        if (mxGetNumberOfElements(prhs[i]) != 1) {
            mexErrMsgIdAndTxt(
                "EcgToolbox:c_rocha_features:notScalar",
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
    mwSize i;
    
    // get pointers to the data in the input vectors
    beatList = mxGetPr(prhs[0]);
    rrList = mxGetPr(prhs[1]);
    defIso = mxGetPr(prhs[2]);
    defJay = mxGetPr(prhs[3]);
    endList = mxGetPr(prhs[4]);
    
    // get the dimensions of the first input
    *nrows = (mwSize)mxGetM(prhs[0]);
    *ncols = (mwSize)mxGetN(prhs[0]);
    
    // make sure the input arguments 2-5 have compatible dimensions
    for (i = 1; i < 5; i++) {
        if (MAX(mxGetM(prhs[i]), mxGetN(prhs[i])) != *ncols) {
            mexErrMsgIdAndTxt(
                "EcgToolbox:c_rocha_features:badDimensions",
                "Inputs #1 and #%d must have compatible dimensions.", i + 1);
        }
    }
    
    // get the sampling frequency
    sampFreq = mxGetScalar(prhs[5]);
    
    // get the gain of the input
    inputGain = mxGetScalar(prhs[6]);
}

/*=========================================================================
 * Hanlde output arguments 
 *=======================================================================*/
void handleOutputs(int nlhs, mxArray *plhs[])
{
    // get a pointer to the first output
    plhs[0] = mxCreateDoubleMatrix(2, qrsLen, mxREAL);
    outSTdev = mxGetPr(plhs[0]);
    
    // get a pointer to the first output
    plhs[1] = mxCreateDoubleMatrix(NUM_HERM_COEF, qrsLen, mxREAL);
    outCoeff1 = mxGetPr(plhs[1]);
    
    // get a pointer to the second output
    plhs[2] = mxCreateDoubleMatrix(NUM_HERM_COEF, qrsLen, mxREAL);
    outCoeff2 = mxGetPr(plhs[2]);
    
    // optional output 1
    if (nlhs > 3) {
        plhs[3] = mxCreateDoubleMatrix(SEGMENT_SIZE, qrsLen, mxREAL);
        outSeg1 = mxGetPr(plhs[3]);
    }
    else outSeg1 = NULL;
    
    // optional output 2
    if (nlhs > 4) {
        plhs[4] = mxCreateDoubleMatrix(SEGMENT_SIZE, qrsLen, mxREAL);
        outSeg2 = mxGetPr(plhs[4]);
    }
    else outSeg2 = NULL;
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
    thresh = (int)(maFilter.gain * inputGain * 1.25 / sampFreq);
    delay = (int)ceil(maFilter.delay + deFilter.delay);
    searchL1 = (int)(0.02 * sampFreq);
    searchL2 = (int)(0.12 * sampFreq);
    
    // initialize buffers
    beatBuf = (int *)mxMalloc(frameSize * sizeof(int));
    seg1Buf = (double *)mxMalloc(SEGMENT_SIZE * sizeof(double));
    seg2Buf = (double *)mxMalloc(SEGMENT_SIZE * sizeof(double));
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
    
    // process one input sample at a time
    for (i = 0; i < qrsLen; i++) {
        onNewBeat();
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
    mxFree(seg1Buf);
    mxFree(seg2Buf);
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
