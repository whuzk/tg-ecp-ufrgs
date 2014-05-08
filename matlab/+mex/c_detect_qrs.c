/*=========================================================================
 * c_detect_qrs.c
 * 
 *  Title: real-time detection of QRS complex using integer arithmetic
 *  Author:     Diego Sogari
 *  Modified:   07/May/2014
 *
 *  Intputs:
 *      1. ECG signal*
 *      2. sampling frequency (in Hz)*
 *
 *  Outputs:
 *      1. location of detected R peaks (in samples)*
 *      2. location of R peaks detected by search back (in samples)
 *      3. history of signal threshold
 *      4. history of secondary signal threshold
 *      5. history of averaged R-R intervals (in samples)
 *      6. history of secondary averaged R-R intervals (in samples)
 *
 *  *required
 *
 *=======================================================================*/
#include <math.h>
#include "mex.h"
#include "c_mexutils.h"
#include "c_timeutils.h"

/*=========================================================================
 * Abbreviations
 *  LEN     length
 *  BUF     buffer
 *  MIN     minimum
 *  MAX     maximum
 *  DEF     default
 *  SEC     seconds
 *  DEL     delay
 *  WIN     window
 *  FREQ    frequency
 *  SIG     signal
 *  QRS     QRS complex
 *  RR      R-R interval
 *  AMP     amplitude
 *  IDX     index
 *  HIST    history
 *  THR     threshold
 *  EST     estimation
 *  FILT    filtered signal
 *=======================================================================*/

/*=========================================================================
 * Constants
 *=======================================================================*/
#define MIN_INPUTS      2       // minimum number of input arguments
#define MAX_INPUTS      2       // maximum number of input arguments
#define MIN_OUTPUTS     1       // minimum number of output arguments
#define MAX_OUTPUTS     6       // maximum number of output arguments

#define MIN_SAMP_FREQ   50*2    // minimum sampling frequency
#define MAX_SAMP_FREQ   60*40   // maximum sampling frequency

#define RRI_BUFLEN      8       // eight most recent RR

/*=========================================================================
 * Input and output variables
 *=======================================================================*/
static double *inputSig;        // input signal
static mwSize inputLen;         // input signal length
static double sampFreq;         // sampling frequency
static double *qrs1Hist;        // history of QRS location
static double *qrs2Hist;        // history of QRS location (search back)
static double *th1Hist;         // history of signal threshold
static double *th2Hist;         // history of secondary signal threshold
static double *rr1Hist;         // history of averaged RR
static double *rr2Hist;         // history of secondary averaged RR

/*=========================================================================
 * Buffer variables
 *=======================================================================*/
static unsigned int ci;             // current sample index
static mwSize rr1Buf[RRI_BUFLEN];   // primary RR buffer
static mwSize rr2Buf[RRI_BUFLEN];   // secondary RR buffer
static int *detBuf;                 // detection buffer
static mwSize detBufLen;            // detection buffer length

/*=========================================================================
 * Fast lookup for detection buffers
 *=======================================================================*/
#define rr1b(I) (rr1Buf[(I)&(RRI_BUFLEN-1)])
#define rr2b(I) (rr2Buf[(I)&(RRI_BUFLEN-1)])
#define detb(I) (detBuf[(ci+(I))&(detBufLen-1)])

/*=========================================================================
 * Detection variables
 *=======================================================================*/
static int sigThreshold;        // signal threshold
static int signalLevel;         // signal level
static int noiseLevel;          // noise level
static int estRatio;            // ratio of signal/noise level estimation
static int peakAmp;             // currently detected peak amplitude
static mwSize rrIntMean1;       // primary running average of RR
static mwSize rrIntMean2;       // secondary running average of RR
static mwSize rrIntMiss;        // interval for QRS to be assumed as missed
static mwSize trainCountDown;   // countdown of training period
static mwSize qrsHalfLength;    // half the length of a QRS (predefined)
static mwSize twaveTolerance;   // tolerance for T-wave detection (predefined)
static mwSize refractoryPeriod; // length of refractory period (predefined)
static mwSize qrsCount1;        // count of QRS complex in main QRS buffer
static mwSize qrsCount2;        // count of QRS complex in QRS buffer 2
static mwSize lastQrsIdx;       // location of last detected QRS
static mwSize searchBackIdx;    // location of searchback starting point
static mwSize peakIdx;          // currently detected peak location
static bool isSignalRising;     // flag to indicate a rise in the signal

/*=========================================================================
 * Check if index is valid in detection buffer 
 *=======================================================================*/
bool isvalidindex(mwSize idx)
{
    return (idx <= 0 && 0 < detBufLen + idx);
}

/*=========================================================================
 * Iterative estimator 
 *=======================================================================*/
int estimate(int x, int log2r, int v)
{
    return (((1 << log2r) - 1) * x + v) >> log2r;
}

/*=========================================================================
 * Get the maximum slope on the detection signal 
 *=======================================================================*/
int maxdiff(mwSize start, mwSize len)
{
    int d = 0;
    int newd;
    
    for (mwSize i = 1; i < len; i++) {
        newd = detb(start + i) - detb(start + i - 1);
        if (newd > d) {
            d = newd;
        }
    }
    return d;
}

/*=========================================================================
 * Find the position of the maximum value on the detection signal 
 *=======================================================================*/
mwSize findmax(mwSize start, mwSize len)
{
    mwSize pos = 0;
    int y, newy;
    
    y = detb(start);
    for (mwSize i = 1; i < len; i++) {
        newy = detb(start + i);
        if (newy > y) {
            y = newy;
            pos = i;
        }
    }
    return start + pos;
}

/*=========================================================================
 * Check for T wave 
 *=======================================================================*/
bool istwave(mwSize candQrs, mwSize lastQrs)
{
    mwSize rr = candQrs - lastQrs;
    
    if (rr <= refractoryPeriod) {
        return true;
    }
    else if (rr > twaveTolerance) {
        return false;
    }
    else {
        // calculate starting indices
        mwSize start1 = candQrs - qrsHalfLength + 1;
        mwSize start2 = lastQrs - qrsHalfLength + 1;
        
        // check index validity
        if (isvalidindex(start1) && isvalidindex(start2)) {
            // check condition for T wave
            int slope1 = maxdiff(start1, qrsHalfLength);
            int slope2 = maxdiff(start2, qrsHalfLength);
            return (slope1 < (slope2 >> 1));
        }
        else return false;
    }
}

/*=========================================================================
 * Peak detection 
 *=======================================================================*/
bool detectPeak(int newAmp)
{
    static mwSize lastPeakIdx = 0;      // index of the last peak
    static int lastPeakAmp = 0;       // amplitude of the last peak
    
    // check if the new amplitude is greater than that of the last peak
    if (newAmp > lastPeakAmp) {
        // signalize beginning or continuation of positive slope
        isSignalRising = true;
        // update peak info
        lastPeakAmp = newAmp;
        lastPeakIdx = 0;
        return false;
    }
    else if (!isSignalRising) {
        // update current amplitude
        lastPeakAmp = newAmp;
        return false;
    }
    else if (newAmp < (lastPeakAmp >> 1)) {
        // report the new peak amplitude and location
        peakAmp = lastPeakAmp;
        peakIdx = lastPeakIdx - 1;
        // reset state of positive slope
        isSignalRising = false;
        // reset peak info
        lastPeakAmp = newAmp;
        lastPeakIdx = 0;
        // signalize detection
        return true;
    }
    else {
        lastPeakIdx--;
        return false;
    }
}

/*=========================================================================
 * Peak evaluation 
 *=======================================================================*/
bool evaluatePeak()
{
    if (trainCountDown > 0) {
        if (peakAmp >= sigThreshold) {
            signalLevel = max(signalLevel, peakAmp);
            lastQrsIdx = peakIdx;
        }
        else {
            noiseLevel = max(noiseLevel, peakAmp);
        }
        return false;
    }
    else if (peakAmp >= sigThreshold) {
        signalLevel = estimate(signalLevel, estRatio, peakAmp);
        return !istwave(peakIdx, lastQrsIdx);
    }
    else {
        noiseLevel = estimate(noiseLevel, estRatio, peakAmp);
        return false;
    }
}

/*=========================================================================
 * Search back procedure 
 *=======================================================================*/
bool searchBack()
{
    if ((!isSignalRising || detb(0) < (sigThreshold >> 1)) &&
            searchBackIdx != 0 && 0 >= searchBackIdx + rrIntMiss) {
        // calculate starting point for the search
        mwSize begin = 0 - rrIntMean2 + 1;
        
        if (isvalidindex(begin)) {
            // search back and locate the max in this interval
            peakIdx = findmax(begin, rrIntMean2);
            peakAmp = detb(peakIdx);
            
            // check if candidate peak is from qrs
            if (peakAmp < sigThreshold &&
                    peakAmp >= (sigThreshold >> 1) &&
                    !istwave(peakIdx, lastQrsIdx)) {
                // adjust signal level
                signalLevel = estimate(signalLevel, 2, peakAmp);
                // signalize qrs detection
                return true;
            }
            else {
                // reduce levels by half
                signalLevel >>= 1;
                noiseLevel >>= 1;
                // postpone searchback
                searchBackIdx += rrIntMean2;
                return false;
            }
        }
        else return false;
    }
    else return false;
}

/*=========================================================================
 * Update RR interval info 
 *=======================================================================*/
void updateRRInfo(mwSize newRR)
{
    static bool initialized = false;
    static mwSize rrIntLow, rrIntHigh;
    static mwSize rrIdx1 = 0, rrIdx2 = 0;
    static mwSize y1, y2;
    
    if (!initialized) {
        for (mwSize i = 0; i < RRI_BUFLEN; i++) {
            rr1b(i) = newRR;
            rr2b(i) = newRR;
        }
        y1 = y2 = newRR << 3;
        initialized = true;
    }
    else {
        y1 += newRR - rr1b(rrIdx1);
        rr1b(rrIdx1++) = newRR;
        if (rrIntLow <= newRR && newRR <= rrIntHigh) {
            y2 += newRR - rr2b(rrIdx2);
            rr2b(rrIdx2++) = newRR;
        }
    }
    
    // update RR averages
    rrIntMean1 = y1 >> 3;
    rrIntMean2 = y2 >> 3;
    rrIntLow = rrIntMean2 - ((y2 - rrIntMean2) >> 6);   // 0.89*RRmean
    rrIntHigh = rrIntLow + ((y2 + rrIntMean2) >> 5);    // 1.17*RRmean
    rrIntMiss = rrIntHigh + (rrIntMean2 >> 2);          // 1.67*RRmean
}

/*=========================================================================
 * Detect QRS complex 
 *=======================================================================*/
bool detectQrs(bool *detectedBySearchback)
{
    bool wasPeakDetected;   // flag to indicate that a peak was detected
    bool wasQrsDetected;    // flag to indicate that a qrs was detected
    
    // detect peak
    wasPeakDetected = detectPeak(detb(0));
    
    // detect qrs
    wasQrsDetected = wasPeakDetected && evaluatePeak();
    
    // search back
    *detectedBySearchback = (!wasQrsDetected) && searchBack();
    
    // update threshold
    sigThreshold = noiseLevel + (abs(signalLevel - noiseLevel) >> 2);
    
    // update qrs info
    if (wasQrsDetected || *detectedBySearchback) {
        updateRRInfo(peakIdx - lastQrsIdx);
        lastQrsIdx = peakIdx;
        searchBackIdx = peakIdx;
    }
    
    // update training countdown
    if (trainCountDown > 0) {
        trainCountDown--;
        // decrease estimation ratio for times beyond the training period
        if (trainCountDown == 0) {
            estRatio = estRatio + 1;
        }
    }
    
    // update indices
    lastQrsIdx--;
    searchBackIdx--;
    
    return wasQrsDetected || *detectedBySearchback;
}

/*=========================================================================
 * Update output buffers
 *=======================================================================*/
void updateOutputs(bool qrsDetected, bool qrsDetectedBySearchback)
{
    // detect qrs history
    if (qrsDetected && qrs1Hist != NULL) {
        qrs1Hist[qrsCount1++] = ci + peakIdx + 1;
    }
    
    // update qrs2 history
    if (qrsDetectedBySearchback && qrs2Hist != NULL) {
        qrs2Hist[qrsCount2++] = ci + peakIdx + 1;
    }
    
    // update main threshold history
    if (th1Hist != NULL) {
        th1Hist[ci] = sigThreshold;
    }
    
    // update threshold history 2
    if (th2Hist != NULL) {
        th2Hist[ci] = sigThreshold >> 1;
    }
    
    // update RR interval history
    if (rr1Hist != NULL) {
        rr1Hist[ci] = rrIntMean1;
    }
    
    // update RR interval history 2
    if (rr2Hist != NULL) {
        rr2Hist[ci] = rrIntMean2;
    }
}

/*=========================================================================
 * Event triggered for an arriving sample
 *=======================================================================*/
void onNewSample(double sample)
{
    bool qrsDetected, qrsDetectedBySearchback;
    
    // update detection buffer
    detb(0) = (int)sample;
    
    // detect qrs
    qrsDetected = detectQrs(&qrsDetectedBySearchback);
    
    // update outputs
    updateOutputs(qrsDetected, qrsDetectedBySearchback);
    
    // increment global index
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
            "EcgToolbox:c_detect_qrs:nrhs",
            "%d input(s) required and %d optional",
            MIN_INPUTS, MAX_INPUTS - MIN_INPUTS);
    }
    // check for proper number of output arguments
    if (nlhs < MIN_OUTPUTS || nlhs > MAX_OUTPUTS) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_detect_qrs:nlhs",
            "%d output(s) required and %d optional",
            MIN_OUTPUTS, MAX_OUTPUTS - MIN_OUTPUTS);
    }
    // make sure the first input argument is a vector
    if (mxGetM(prhs[0]) != 1 && mxGetN(prhs[0]) != 1) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_detect_qrs:notVector",
            "First input must be a vector.");
    }
    // make sure the first input argument is of type double
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_detect_qrs:notDouble",
            "First input must be of type double.");
    }
    // make sure the remaining arguments are all scalars
    for (mwSize i = 1; i < nrhs; i++) {
        if (mxGetNumberOfElements(prhs[i]) != 1) {
            mexErrMsgIdAndTxt(
                "EcgToolbox:c_detect_qrs:notScalar",
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
    // get a pointer to the data in the input vector
    inputSig = mxGetPr(prhs[0]);
    
    // get the dimensions of the input vector
    *nrows = (mwSize)mxGetM(prhs[0]);
    *ncols = (mwSize)mxGetN(prhs[0]);
    
    // get the sampling frequency
    sampFreq = (int)mxGetScalar(prhs[1]);
    
    // make sure the sampling frequency is within pre-defined limits
    if (sampFreq < MIN_SAMP_FREQ || sampFreq > MAX_SAMP_FREQ) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_detect_qrs:badSampFreq",
            "Sampling frequency must be between %d and %d.",
            MIN_SAMP_FREQ, MAX_SAMP_FREQ);
    }
}

/*=========================================================================
 * Hanlde output arguments 
 *=======================================================================*/
void handleOutputs( int nlhs, mxArray *plhs[],
                    mwSize nrows, mwSize ncols)
{
    double *outVectors[MAX_OUTPUTS] = {NULL};
    
    // create the output vectors
    for (mwSize i = 0; i < nlhs; i++) {
        plhs[i] = mxCreateDoubleMatrix(nrows, ncols, mxREAL);
        outVectors[i] = mxGetPr(plhs[i]);
    }
    
    // get pointers to the output vectors
    qrs1Hist = outVectors[0];
    qrs2Hist = outVectors[1];
    th1Hist = outVectors[2];
    th2Hist = outVectors[3];
    rr1Hist = outVectors[4];
    rr2Hist = outVectors[5];
}

/*=========================================================================
 * The initializarion routine
 *=======================================================================*/
void init()
{
    // initialize sample index
    ci = 0;
    
    // threshold, indices, etc.
    sigThreshold = 0;
    signalLevel = 0;
    noiseLevel = 0;
    estRatio = 2;
    peakIdx = 0;
    peakAmp = 0;
    qrsCount1 = 0;
    qrsCount2 = 0;
    lastQrsIdx = 0;
    searchBackIdx = 0;
    isSignalRising = false;
    
    // predefined limits
    trainCountDown = (mwSize)(2 * sampFreq);
    qrsHalfLength = (mwSize)round(0.10 * sampFreq);
    twaveTolerance = (mwSize)round(0.36 * sampFreq);
    refractoryPeriod = (mwSize)round(0.20 * sampFreq);
    
    // create buffer for the filtered signal
    detBufLen = 1 << (1 + ilogb(2 * sampFreq - 1));
    detBuf = (int *)mxMalloc(detBufLen * sizeof(int));
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
    for (mwSize i = 0; i < inputLen; i++) {
        onNewSample(inputSig[i]);
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
void finalize( int nlhs, mxArray *plhs[],
               mwSize nrows, mwSize ncols)
{
    // adjust the size of qrs history output
    adjustSizeDouble(qrs1Hist, plhs[0], nrows, ncols, qrsCount1);
    
    // adjust the size of qrs2 history output
    if (nlhs > 1) {
        adjustSizeDouble(qrs2Hist, plhs[1], nrows, ncols, qrsCount2);
    }
    
    // deallocate memory
    mxFree(detBuf);
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
    
    // handle output arguments
    handleOutputs(nlhs, plhs, nrows, ncols);
    
    // calculate the length of the signal
    inputLen = max(nrows,ncols);
    
    // make some initializations
    init();
    
    // do the actual processing
    doTheJob();
    
    // perform final adjustments
    finalize(nlhs, plhs, nrows, ncols);
}
