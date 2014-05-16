/*=========================================================================
 * c_extract_beats.c
 * 
 *  Title: extraction of beats
 *  Author:     Diego Sogari
 *  Modified:   10/May/2014
 *
 *  Intputs:
 *      1. input ECG signal
 *      2. list of fiducial points (in matrix form)
 *      3. sampling frequency (in Hz)
 *
 *  Outputs:
 *      1. List of beats (in matrix form)
 *      2. List of corrected of fiducial points (in matrix form)
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
 *  SEC     seconds
 *  FREQ    frequency
 *  SIG     signal
 *  IDX     index
 *  FDP     fiducial points
 *=======================================================================*/

/*=========================================================================
 * Constants
 *=======================================================================*/
#define MIN_INPUTS      3       // minimum number of input arguments
#define MAX_INPUTS      3       // maximum number of input arguments
#define MIN_OUTPUTS     2       // minimum number of output arguments
#define MAX_OUTPUTS     2       // maximum number of output arguments

#define BL_NUM_PTS      5       // number of points for baseline removal
#define FDP_NUM_COL     7       // number of columns in the FDP matrix
#define HALF_FRAME      0.6     // half the size of the frame

/*=========================================================================
 * Input and output variables
 *=======================================================================*/
static double *inputSig;        // input signal
static mwSize inputLen;         // input signal length
static double *fdpList;         // list of fiducial points
static mwSize qrsLen;           // number of beats
static double sampFreq;         // sampling frequency
static double *beatList;        // list of extracted beats
static double *fdpListOut;      // list of corrected fiducial points

/*=========================================================================
 * Buffer variables
 *=======================================================================*/
static unsigned int ci;             // current sample index
static mwSize bi;                   // current beat index
static double *buffer;              // buffer for the signal
static mwSize bufLen;               // length of the buffer
static mwSize frameSize;            // length of beat frame

/*=========================================================================
 * Fast lookup for buffers
 *=======================================================================*/
#define buf(I)      (buffer[(ci+(I))&(bufLen-1)])
#define fdp(I,J)    (fdpList[(I)+qrsLen*(J)])
#define beat(I,J)   (beatList[(I)*frameSize+(J)])
#define fdpout(I,J) (fdpListOut[(I)+qrsLen*(J)])

/*=========================================================================
 * Check if index is valid in detection buffers
 *=======================================================================*/
bool isvalidindex(mwSize idx)
{
    return (idx <= 0 && 0 < bufLen + idx);
}

/*=========================================================================
 * Copy from buffer to output
 *=======================================================================*/
void copyBuffer()
{
    mwSize startIdx, endIdx, rpeakIdx, absIdx1, absIdx2;
    mwSize len1, len2, r, half, beatsize, pad, cut1, cut2;
    
    startIdx = (mwSize)fdp(bi, 0) - 1 - ci;
    rpeakIdx = (mwSize)fdp(bi, 3) - 1 - ci;
    endIdx = (mwSize)fdp(bi, FDP_NUM_COL-1) - 1 - ci;
    
    half = frameSize >> 1;
    beatsize = endIdx - startIdx + 1;
    r = rpeakIdx - startIdx + 1;
    
    pad = max(0, half - r + 1);
    cut1 = max(0, r - 1 - half);
    cut2 = max(0, beatsize - r - half);
    
    startIdx += cut1;
    endIdx -= cut2;
    if (startIdx > endIdx || !isvalidindex(startIdx) ||
            !isvalidindex(endIdx)) {
        return;
    }
    
    absIdx1 = (ci + startIdx) & (bufLen - 1);
    absIdx2 = (ci + endIdx) & (bufLen - 1);
    
    if (absIdx1 <= absIdx2) {
        len1 = (absIdx2 - absIdx1 + 1) * sizeof(double);
        memcpy(&beat(bi, pad), buffer + absIdx1, len1);
    }
    else {
        len1 = (bufLen - absIdx1) * sizeof(double);
        len2 = (absIdx2 + 1) * sizeof(double);
        memcpy(&beat(bi, pad), buffer + absIdx1, len1);
        memcpy(&beat(bi, pad + bufLen - absIdx1), buffer, len2);
    }
}

/*=========================================================================
 * Update outputs
 *=======================================================================*/
void updateOutputs()
{
    mwSize center;
    
    if (bi < qrsLen) {
        copyBuffer();
        
        center = (frameSize >> 1) + 1;
        fdpout(bi,0) = max(1, center - (fdp(bi, 3) - fdp(bi,0)));
        fdpout(bi,1) = max(1, center - (fdp(bi, 3) - fdp(bi,1)));
        fdpout(bi,2) = max(1, center - (fdp(bi, 3) - fdp(bi,2)));
        fdpout(bi,3) = center;
        fdpout(bi,4) = min(frameSize, center + (fdp(bi,4) - fdp(bi, 3)));
        fdpout(bi,5) = min(frameSize, center + (fdp(bi,5) - fdp(bi, 3)));
        fdpout(bi,6) = min(frameSize, center + (fdp(bi,6) - fdp(bi, 3)));
        
        bi++;
    }
}

/*=========================================================================
 * Perform the baseline removal procedure
 *=======================================================================*/
void removeBaseline()
{
    double x[2], y[2], p[2], mean;
    mwSize startIdx, endIdx;
    
    // calculate start and end points
    startIdx = (mwSize)fdp(bi, 0) - 1 - (BL_NUM_PTS >> 1) - ci;
    endIdx = (mwSize)fdp(bi, FDP_NUM_COL - 1) - 1 + (BL_NUM_PTS >> 1) - ci;
    
    if (!isvalidindex(startIdx) || !isvalidindex(endIdx)) {
        return;
    }
    
    // average of first BL_NUM_PTS samples
    mean = 0.0;
    for (mwSize i = startIdx; i < startIdx + BL_NUM_PTS; i++) {
        mean += buf(i);
    }
    x[0] = 0.0;
    y[0] = mean / BL_NUM_PTS;
    
    // average of last BL_NUM_PTS samples
    mean = 0.0;
    for (mwSize i = endIdx - BL_NUM_PTS + 1; i <= endIdx; i++) {
        mean += buf(i);
    }
    x[1] = endIdx - startIdx;
    y[1] = mean / BL_NUM_PTS;
    
    // create polynomial fit
    first_order_fit(x, y, p);
    
    // calculate new amplitudes
    for (mwSize i = startIdx; i <= endIdx; i++) {
        buf(i) -= p[0] * (i - startIdx) + p[1];
    }
}

/*=========================================================================
 * Event triggered for an arriving sample
 *=======================================================================*/
void onNewSample(double sample, bool newBeat)
{
    // update buffer
    buf(0) = sample;
    
    // simulate beat event
    if (newBeat) {
        // perform baseline removal
        removeBaseline();
        // update the outputs
        updateOutputs();
    }
    
    // increment sample index
    ci++;
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
            "EcgToolbox:c_extract_beats:nrhs",
            "%d input(s) required and %d optional",
            MIN_INPUTS, MAX_INPUTS - MIN_INPUTS);
    }
    // check for proper number of output arguments
    if (nlhs < MIN_OUTPUTS || nlhs > MAX_OUTPUTS) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_extract_beats:nlhs",
            "%d output(s) required and %d optional",
            MIN_OUTPUTS, MAX_OUTPUTS - MIN_OUTPUTS);
    }
    // make sure all input arguments are of type double
    for (mwSize i = 0; i < nrhs; i++) {
        if (!mxIsDouble(prhs[i]) || mxIsComplex(prhs[i])) {
            mexErrMsgIdAndTxt(
                "EcgToolbox:c_extract_beats:notDouble",
                "Input #%d must be of type double.", i + 1);
        }
    }
    // make sure the first input argument is a vector
    if (mxGetM(prhs[0]) != 1 && mxGetN(prhs[0]) != 1) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_extract_beats:notVector",
            "Input #1 must be a vector.");
    }
    // make sure the second input argument has a correct number of columns
    if (mxGetM(prhs[1]) == 1 || mxGetN(prhs[1]) != FDP_NUM_COL) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_extract_beats:notMatrix",
            "Input #2 must have %d columns.", FDP_NUM_COL);
    }
    // make sure the remaining input arguments are all scalars
    for (mwSize i = 2; i < nrhs; i++) {
        if (mxGetNumberOfElements(prhs[i]) != 1) {
            mexErrMsgIdAndTxt(
                "EcgToolbox:c_extract_beats:notScalar",
                "Input #%d must be a scalar.", i + 1);
        }
    }
}

/*=========================================================================
 * Hanlde input arguments 
 *=======================================================================*/
void handleInputs( int nrhs, const mxArray *prhs[],
                   mwSize *nrows1, mwSize *ncols1,
                   mwSize *nrows2, mwSize *ncols2)
{
    // get pointers to the data in the inputs
    inputSig = mxGetPr(prhs[0]);
    fdpList = mxGetPr(prhs[1]);
    
    // get the dimensions of the first input
    *nrows1 = (mwSize)mxGetM(prhs[0]);
    *ncols1 = (mwSize)mxGetN(prhs[0]);
    
    // get the dimensions of the second input
    *nrows2 = (mwSize)mxGetM(prhs[1]);
    *ncols2 = (mwSize)mxGetN(prhs[1]);
    
    // get the sampling frequency
    sampFreq = mxGetScalar(prhs[2]);
}

/*=========================================================================
 * Hanlde output arguments 
 *=======================================================================*/
void handleOutputs( int nlhs, mxArray *plhs[],
                    mwSize nrows, mwSize ncols)
{
    // get a pointer to the first output
    plhs[0] = mxCreateDoubleMatrix(frameSize, qrsLen, mxREAL);
    beatList = mxGetPr(plhs[0]);
    
    // get a pointer to the second output
    plhs[1] = mxCreateDoubleMatrix(nrows, ncols, mxREAL);
    fdpListOut = mxGetPr(plhs[1]);
}

/*=========================================================================
 * The initializarion routine
 *=======================================================================*/
void init()
{
    // initialize indices
    ci = 0;
    bi = 0;
    
    // calculate the frame size
    frameSize = 2 * (int)round(HALF_FRAME * sampFreq) + 1;
    
    // create the buffer
    bufLen = 1 << (1 + ilogb(2 * sampFreq - 1));
    buffer = (double *)mxMalloc(bufLen * sizeof(double));
}

/*=========================================================================
 * The main routine 
 *=======================================================================*/
void doTheJob()
{
    double time;
    mwSize newbeatend;
    mwSize beatCount = 0;
    
    // start time counter
    tic();
    
    // initialize qrs variables
    newbeatend = (mwSize)fdp(beatCount, FDP_NUM_COL - 1) - 1 + BL_NUM_PTS;
    
    // process one input sample at a time
    for (mwSize i = 0; i < inputLen; i++) {
        onNewSample(inputSig[i], newbeatend == i);
        if (newbeatend == i && beatCount < qrsLen - 1) {
            newbeatend = (mwSize)fdp(++beatCount, FDP_NUM_COL - 1) - 1 + BL_NUM_PTS;
        }
    }
    
    // stop time counter
    time = toc();
    
    // display time statistics
    mexPrintf("Total processing time: %.2f ms\n", 1000 * time);
    mexPrintf("Average time per sample: %.2f ns\n", 1000000000 * time / inputLen);
}

/*=========================================================================
 * The finalization routine 
 *=======================================================================*/
void finalize()
{
    // deallocate memory
    mxFree(buffer);
}

/*=========================================================================
 * The gateway function 
 *=======================================================================*/
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    mwSize nrows1;
    mwSize ncols1;
    mwSize nrows2;
    mwSize ncols2;
    
    // check argument correctness
    checkArgs(nlhs, plhs, nrhs, prhs);
    
    // handle input arguments
    handleInputs(nrhs, prhs, &nrows1, &ncols1, &nrows2, &ncols2);
    
    // calculate the length of the signal
    inputLen = max(nrows1,ncols1);
    
    // get the number of beats in the list
    qrsLen = nrows2;
    
    // make some initializations
    init();
    
    // handle output arguments
    handleOutputs(nlhs, plhs, nrows2, ncols2);
    
    // do the actual processing
    doTheJob();
    
    // perform final adjustments
    finalize();
}
