/*=========================================================================
 * c_mohebbi_ijpoints.c
 * 
 *  Title: detection of isoelectric and J points
 *  Author:     Diego Sogari
 *  Modified:   15/May/2014
 *
 *  Intputs:
 *      1. List of beats (in matrix form)
 *      2. Default value of isoelectric points
 *      3. Default value of J points
 *      4. Sampling frequency (in Hertz)
 *      5. Gain of the input signal
 *
 *  Outputs:
 *      1. List of isoelectric points
 *      2. List of J points
 *
 *  *all inputs and outputs required
 *
 *=======================================================================*/
#include <math.h>
#include "mex.h"
#include "c_mathutils.h"
#include "c_intfilter.h"
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
 *  MAF     moving-average filter
 *  DER     derivative filter
 *=======================================================================*/

/*=========================================================================
 * Constants
 *=======================================================================*/
#define MIN_INPUTS      5       // minimum number of input arguments
#define MAX_INPUTS      5       // maximum number of input arguments
#define MIN_OUTPUTS     2       // minimum number of output arguments
#define MAX_OUTPUTS     2       // maximum number of output arguments

#define MAF_WIDTH       0.05    // width of the moving-average (in seconds)
#define DER_ORDERN      1       // derivative filter order N
#define DER_ORDERM      0       // derivative filter order M

/*=========================================================================
 * Input and output variables
 *=======================================================================*/
static double *beatList;        // list of beats
static double *defIso;          // list of default isoelectric points
static double *defJay;          // list of default J points
static mwSize qrsLen;           // number of beats
static double sampFreq;         // sampling frequency
static double inputGain;        // gain of the input signal
static double *outIso;          // output list of isoelectric points
static double *outJay;          // output list of J points

/*=========================================================================
 * Filter and buffer variables
 *=======================================================================*/
static mwSize bi;               // current beat index
static intfobject maFilter;     // moving-average filter object
static intfobject deFilter;     // derivative filter object
static int *beatBuf;            // buffer for the current beat
static mwSize frameSize;        // size of one beat
static mwSize isoPoint;         // isoelectric point
static mwSize jayPoint;         // jay point
static int thresh;              // thresold for detection
static mwSize delay;            // sum of filter delays (in samples)
static mwSize center;           // center of the beat
static mwSize searchL1;         // search limit 1
static mwSize searchL2;         // search limit 2

/*=========================================================================
 * Lookup for the beat list
 *=======================================================================*/
#define beat(I,J)   (beatList[(I)*frameSize+(J)])

/*=========================================================================
 * Update outputs
 *=======================================================================*/
void updateOutputs()
{
    outIso[bi] = isoPoint + 1;
    outJay[bi] = jayPoint + 1;
    bi++;
}

/*=========================================================================
 * Filter one beat signal
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
mwSize get_point(mwSize istart, mwSize iend, mwSize len, int thr,
        mwSize def)
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
        if (abs(beatBuf[i]) > thr) {
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
 * Event triggered for a new beat
 *=======================================================================*/
void onNewBeat()
{
    mwSize istart, iend;
    
    // filter the current beat
    filterBeat();
    
    // detect the isoelectric point
    istart = center - searchL1;
    iend = center - searchL2;
    isoPoint = get_point(istart, iend, searchL1, thresh, (mwSize)defIso[bi]);
    
    // detect the J point
    istart = center + searchL1;
    iend = center + searchL2;
    jayPoint = get_point(istart, iend, searchL1, thresh, (mwSize)defJay[bi]);
    
    // update output vectors
    updateOutputs();
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
            "EcgToolbox:c_mohebbi_ijpoints:nrhs",
            "%d input(s) required and %d optional",
            MIN_INPUTS, MAX_INPUTS - MIN_INPUTS);
    }
    // check for proper number of output arguments
    if (nlhs < MIN_OUTPUTS || nlhs > MAX_OUTPUTS) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_mohebbi_ijpoints:nlhs",
            "%d output(s) required and %d optional",
            MIN_OUTPUTS, MAX_OUTPUTS - MIN_OUTPUTS);
    }
    // make sure all input arguments are of type double
    for (i = 0; i < nrhs; i++) {
        if (!mxIsDouble(prhs[i]) || mxIsComplex(prhs[i])) {
            mexErrMsgIdAndTxt(
                "EcgToolbox:c_mohebbi_ijpoints:notDouble",
                "Input #%d must be of type double.", i + 1);
        }
    }
    // make sure the first input argument has more than one line
    if (mxGetM(prhs[0]) == 1) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_mohebbi_ijpoints:badDimensions",
            "Input #1 must have more than one line.");
    }
    // make sure the second input argument has exactly one column
    if (mxGetN(prhs[1]) != 1) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_mohebbi_ijpoints:badDimensions",
            "Input #2 must have exactly one column.");
    }
    // make sure the third input argument has exactly one column
    if (mxGetN(prhs[2]) != 1) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_mohebbi_ijpoints:badDimensions",
            "Input #3 must have exactly one column.");
    }
    // make sure the remaining input arguments are all scalars
    for (i = 3; i < nrhs; i++) {
        if (mxGetNumberOfElements(prhs[i]) != 1) {
            mexErrMsgIdAndTxt(
                "EcgToolbox:c_mohebbi_ijpoints:notScalar",
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
    defIso = mxGetPr(prhs[1]);
    defJay = mxGetPr(prhs[2]);
    
    // get the dimensions of the first input
    *nrows = (mwSize)mxGetM(prhs[0]);
    *ncols = (mwSize)mxGetN(prhs[0]);
    
    // make sure the second and third input arguments have compatible dimensions
    if (mxGetM(prhs[1]) != *ncols || mxGetM(prhs[2]) != *ncols) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_mohebbi_ijpoints:badDimensions",
            "Inputs #1, #2 and #3 must have compatible dimensions.");
    }
    
    // get the sampling frequency
    sampFreq = mxGetScalar(prhs[3]);
    
    // get the gain of the input
    inputGain = mxGetScalar(prhs[4]);
}

/*=========================================================================
 * Hanlde output arguments 
 *=======================================================================*/
void handleOutputs(int nlhs, mxArray *plhs[])
{
    // get a pointer to the first output
    plhs[0] = mxCreateDoubleMatrix(qrsLen, 1, mxREAL);
    outIso = mxGetPr(plhs[0]);
    
    // get a pointer to the second output
    plhs[1] = mxCreateDoubleMatrix(qrsLen, 1, mxREAL);
    outJay = mxGetPr(plhs[1]);
}

/*=========================================================================
 * The initializarion routine
 *=======================================================================*/
void init()
{
    // initialize beat index
    bi = 0;
    
    // initialize filter objeects
    initintfilter(maFilter);
    initintfilter(deFilter);
    
    // design moving-average filter
    design_maverage(&maFilter, sampFreq, MAF_WIDTH);
    
    // design derivative filter
    design_derivative(&deFilter, DER_ORDERN, DER_ORDERM);
    
    // initialize numeric data
    thresh = (int)(maFilter.gain * inputGain * 1.25 / sampFreq);
    delay = (int)ceil(maFilter.delay + deFilter.delay);
    center = frameSize >> 1;
    searchL1 = (int)(0.02 * sampFreq);
    searchL2 = (int)(0.12 * sampFreq);
    
    // create the buffers
    beatBuf = (int *)mxMalloc(frameSize * sizeof(int));
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
    // deallocate memory of the filter objects
    endintfilter(maFilter);
    endintfilter(deFilter);
    
    // deallocate memory
    mxFree(beatBuf);
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
