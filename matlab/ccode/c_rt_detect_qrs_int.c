/*=========================================================================
 * c_rt_detect_qrs_int.c
 * 
 *  Title: real-time detection of QRS complex using integer arithmetic
 *  Author: Diego Sogari
 *  Last modification:  27/Apr/2014
 *
 *  Outputs:
 *      1. ECG signal*
 *      2. sampling frequency*
 *      3. mains frequency (is guessed based on sampling frequency)
 *      4. MBD filter order (defaults to 3)
 *
 *  Outputs:
 *      1. filtered signal*
 *      2. location of detected R peak (in samples)*
 *      3. location of R peaks detected by search back (in samples)
 *      4. history of main signal threshold
 *      5. history of secondary signal threshold
 *      6. history of estimated R-R intervals (in samples)
 *      7. overall preprocessing delay (in samples)
 *
 *  *required
 *
 *=======================================================================*/
#include <math.h>
#include "mex.h"
#include "c_time_utils.h"

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
 *  LPF     low-pass filter (conv. with second-order backward difference)
 *  WMF     windowed-maximum filter
 *  AGC     automatic gain control
 *  MBD     multiplication of absolute backward differences
 *  MAF     moving-average filter
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
#define MIN_INPUTS      2
#define MAX_INPUTS      4
#define MIN_OUTPUTS     1
#define MAX_OUTPUTS     7

#define MIN_MAINS_FREQ  50
#define MAX_MAINS_FREQ  60
#define DEF_MAINS_FREQ  60

#define MIN_MBD_ORDER   1
#define MAX_MBD_ORDER   8
#define DEF_MBD_ORDER   3

#define ACG_LEN_SEC     0.05
#define MAF_LEN_SEC     0.10

#define LPF_BUFLEN      256     // supports sampling frequencies of up to (256*mainsFreq) Hz
#define WMF_BUFLEN      1024    // 2^nextpow2(MAX_WMF_LIST_SIZE) [determined experimentaly]
#define ACG_BUFLEN      1024    // 2^nextpow2(MAX_MAINS_FREQ * LPF_BUFLEN * ACG_LEN_SEC)
#define MBD_BUFLEN      8       // 2^nextpow2(MAX_MBD_ORDER)
#define MAF_BUFLEN      2048    // 2^nextpow2(MAX_MAINS_FREQ * LPF_BUFLEN * MAF_LEN_SEC)

/*=========================================================================
 * Input variables
 *=======================================================================*/
static short *inputSig;     // input signal
static int sampFreq;        // sampling frequency
static int mainsFreq;       // mains frequency
static int mbdOrder;        // MBD filter order

/*=========================================================================
 * Output variables
 *=======================================================================*/
static int *filtSig;        // filtered signal
static int *qrsHist;        // indices of detected QRS complex
static int *qrs2Hist;       // indices of QRS detected by searchback
static int *thr1Hist;       // history of signal threshold
static int *thr2Hist;       // history of secondary signal threshold
static int *rrintHist;      // history of average RR interval
static double delay;        // overall preprocessing delay

/*=========================================================================
 * Preprocessing variables
 *=======================================================================*/
static mwSize lpfWin;                       // LPF window length
static mwSize wmfWin;                       // WMF window length
static mwSize agcWin;                       // AGC window length
static short agcMin;                        // AGC lower limit
static short agcMax;                        // AGC upper limit
static mwSize mbdWin;                       // MBD window length
static mwSize mafWin;                       // MAF window length
static mwSize mafWinLog2;                   // log2 of MAF window length
static short lpfBuf[LPF_BUFLEN] = {0};      // LPF buffer
static short wmfBuf[WMF_BUFLEN] = {0};      // WMF buffer
static unsigned int wmfAux[WMF_BUFLEN];     // WMF aux buffer
static short agcBuf[ACG_BUFLEN] = {0};      // AGC buffer
static int mbdBuf[MBD_BUFLEN] = {1};        // MBD buffer
static int mafBuf[MAF_BUFLEN] = {0};        // MAF buffer
static int *detBuf;                         // detection buffer
static mwSize detBufLen;                    // detection buffer length

/*=========================================================================
 * Fast lookup for filter buffers
 *=======================================================================*/
#define lpfb(I) (lpfBuf[(I)&(LPF_BUFLEN-1)])
#define wmfb(I) (wmfBuf[(I)&(WMF_BUFLEN-1)])
#define wmfa(I) (wmfAux[(I)&(WMF_BUFLEN-1)])
#define agcb(I) (agcBuf[(I)&(ACG_BUFLEN-1)])
#define mbdb(I) (mbdBuf[(I)&(MBD_BUFLEN-1)])
#define mafb(I) (mafBuf[(I)&(MAF_BUFLEN-1)])
#define detb(I) (detBuf[(I)&(detBufLen-1)])

/*=========================================================================
 * QRS Detection variables
 *=======================================================================*/
static int sigThreshold;        // signal threshold
static int signalLevel;         // signal level
static int noiseLevel;          // noise level
static int peakAmp;             // currently detected peak amplitude
static mwSize rrIntMean;        // running average of RR intervals
static mwSize rrIntMiss;        // interval for qrs to be assumed as missed
static mwSize rrIntLow;         // lower bound for acceptance of new RR
static mwSize rrIntHigh;        // upper bound for acceptance of new RR
static mwSize rrIntMin;         // lowest bound applied to estimated RR
static mwSize rrIntMax;         // highest bound applied to estimated RR
static mwSize estRatio;         // ratio of signal/noise level estimation
static mwSize trainCountDown;   // training period countdown
static mwSize qrsHalfLength;    // half the length of a QRS complex
static mwSize twaveTolerance;   // tolerance for T-wave detection
static mwSize refractoryPeriod; // length of refractory period for QRS detection
static mwSize qrsCount;         // count of QRS complex in main QRS buffer
static mwSize qrsCount2;        // current position in second QRS buffer
static mwSize lastQrsIdx;       // index of last detected QRS complex
static mwSize searchBackIdx;    // index of searchback starting point
static mwSize peakIdx;          // currently detected peak index
static mwSize ci;               // index of the current filtered sample
static bool isSignalRising;     // flag to indicate a rise in the signal

/*=========================================================================
 * Low-pass filter and second-order backward difference
 *=======================================================================*/
short lpf(short sample)
{
    static unsigned short n = 0;
    int y0;
    
    // update filter output: y[n] = x[n] - 2*x[n-M/2] + x[n-M+1]
    y0 = sample - (lpfb(n - (lpfWin >> 1)) << 1) + lpfb(n - lpfWin + 1);
    
    // update filter memory
    lpfb(n++) = sample;
    
    // compute result: y[n]/4
    return y0 >> 2;
}

/*=========================================================================
 * Windowed-maximum
 *=======================================================================*/
short wmf(short sample)
{
    static unsigned short first = 0;
    static unsigned short count = 0;
    static unsigned int i = 0;
    unsigned short j,k;
    
    // search for the first element greater than the current sample
    j = count;
    k = first + count - 1;
    while (j > 0 && sample >= wmfb(k)) {
        j--;
        k--;
    }
    // put the sample next to the element found and adjust the length
    wmfb(k + 1) = sample;
    wmfa(k + 1) = i;
    count = j + 1;
    // check if the first in line has gone out of the window length
    if (count > wmfWin || wmfa(first) == i - wmfWin) {
        first++;
        count--;
    }
    // increment global index
    i++;
    // return the max in the window
    return wmfb(first);
}

/*=========================================================================
 * Automatic gain control
 *=======================================================================*/
short agc(short sample, short gain)
{
    static unsigned short n = 0;
    int y0;
    
    // update filter output: y[n] = x[n-M] * maxG / max(G, minG)
    y0 = (agcb(n - agcWin) * (int)agcMax) / (int)max(gain, agcMin);
    
    // update filter memory
    agcb(n++) = sample;
    
    // compute result
    return (short)min(y0, agcMax);
}

/*=========================================================================
 * Multiplication of absolute backward differences
 *=======================================================================*/
int mdb(short sample)
{
    static unsigned short n = 0;
    int y0 = sample;
    
    // update filter output: y[n] = x[n] * ... * x[n-M+1]
    for (mwSize k = 1; k < mbdWin; k++) {
        y0 *= mbdb(n-k);
    }
    
    // update filter memory
    mbdb(n++) = sample;
    
    // compute result
    return y0;
}

/*=========================================================================
 * Moving-average 
 *=======================================================================*/
int maf(int sample)
{
    static unsigned short n = 0;
    static long long y0 = 0;
    
    // update filter output: y[n] = y[n-1] + x[n] - x[n-M]
    y0 += sample - mafb(n - mafWin);
    
    // update filter memory
    mafb(n++) = sample;
    
    // compute result: y[n]/M
    return (int)min(y0 >> mafWinLog2, INT_MAX);
}

/*=========================================================================
 * Preprocess one input sample 
 *=======================================================================*/
int preprocessSample(short sample)
{
    short gain;
    int sample2;
    
    sample = lpf(sample);
    sample = abs(sample);
    gain = wmf(sample);
    sample = agc(sample, gain);
    sample2 = mdb(sample);
    sample2 = maf(sample2);
    
    return sample2;
}

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
 * Limit a value by lower and upper bounds 
 *=======================================================================*/
int limit(int x, int a, int b)
{
    return min(max(x, a), b);
}

/*=========================================================================
 * Get the maximum slope on the detection signal 
 *=======================================================================*/
int maxdiff(mwSize start, mwSize len)
{
    int d = 0;
    int newd;
    
    for (mwSize i = 1; i < len; i++) {
        newd = detb(start+i)-detb(start+i-1);
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
    int y = detb(start);
    int newy;
    
    for (mwSize i = 1; i < len; i++) {
        newy = detb(start+i);
        if (newy > y) {
            y = newy;
            pos = i;
        }
    }
    return pos;
}

/*=========================================================================
 * Check for T wave 
 *=======================================================================*/
bool istwave(mwSize candQrs, mwSize lastQrs)
{
    // check if the canditate QRS occurs near the previous one
    if (candQrs - lastQrs < twaveTolerance) {
        // calculate starting indices
        mwSize start1 = candQrs - qrsHalfLength + 1;
        mwSize start2 = lastQrs - qrsHalfLength + 1;
        
        // check index validity
        if (isvalidindex(start1) && isvalidindex(start2)) {
            // max slope of waveforms
            int slope1 = maxdiff(ci + start1, qrsHalfLength);
            int slope2 = maxdiff(ci + start2, qrsHalfLength);

            // check condition for T wave
            return (slope1 < (slope2 >> 1));
        }
        else return false;
    }
    else return false;
}

/*=========================================================================
 * Peak detection 
 *=======================================================================*/
bool detectPeak()
{
    static mwSize lastPeakIdx = 0;      // index of the last peak
    static int lastPeakAmp = 0.0;       // amplitude of the last peak
    int newAmp = detb(ci);              // current signal amplitude
    
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
        // update peak info
        lastPeakAmp = newAmp;
        lastPeakIdx = 0;
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
    if (lastQrsIdx == 0 || peakIdx - lastQrsIdx > refractoryPeriod) {
        // adjust levels
        if (peakAmp >= sigThreshold) {
            signalLevel = estimate(signalLevel, estRatio, peakAmp);
        }
        else {
            noiseLevel = estimate(noiseLevel, estRatio, peakAmp);
        }
        
        // assert qrs
        return (trainCountDown == 0 && peakAmp >= sigThreshold &&
                (lastQrsIdx == 0 || !istwave(peakIdx, lastQrsIdx)));
    }
    else return false;
}

/*=========================================================================
 * Search back procedure 
 *=======================================================================*/
bool searchBack()
{
    if ((!isSignalRising) && searchBackIdx != 0 &&
            0 >= searchBackIdx + rrIntMiss) {
        // search back and locate the max in this interval
        mwSize len = min(detBufLen,rrIntMean);
        peakIdx = findmax(ci - len + 1, len) - len + 1;
        peakAmp = detb(ci + peakIdx);
        
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
            //signalLevel = 0.5*signalLevel;
            //noiseLevel = 0.5*noiseLevel;
            // postpone searchback
            searchBackIdx += len;
            return false;
        }
    }
    else return false;
}

/*=========================================================================
 * Update QRS and RR interval info 
 *=======================================================================*/
void updateQrsInfo()
{
    // update RR interval
    if (lastQrsIdx != 0) {
        mwSize newRR = peakIdx - lastQrsIdx;
        if (rrIntLow < newRR && newRR < rrIntHigh) {
            rrIntMean = estimate(rrIntMean, 3, newRR);
            rrIntMean = limit(rrIntMean, rrIntMin, rrIntMax);
        }
    }
    
    // calculate RR missed limit
    rrIntLow = rrIntMean >> 1;
    rrIntHigh = rrIntMean + rrIntLow;
    rrIntMiss = rrIntHigh + (rrIntMean >> 2);
    
    // update indices
    lastQrsIdx = peakIdx;
    searchBackIdx = peakIdx;
}  

/*=========================================================================
 * Detect QRS complex 
 *=======================================================================*/
bool detectQrs(bool *detectedBySearchback)
{
    bool wasPeakDetected;   // flag to indicate that a peak was detected
    bool wasQrsDetected;    // flag to indicate that a qrs was detected
    
    // detect peak
    wasPeakDetected = detectPeak();
    
    // detect qrs
    wasQrsDetected = wasPeakDetected && evaluatePeak();
    
    // search back
    *detectedBySearchback = (!wasQrsDetected) && searchBack();
    
    // update threshold
    sigThreshold = noiseLevel + (abs(signalLevel - noiseLevel) >> 2);
    
    // update qrs info
    if (wasQrsDetected || *detectedBySearchback) {
        updateQrsInfo();
    }
    
    // update training countdown
    if (trainCountDown > 0) {
        trainCountDown--;
        // decrease estimation ratio for times beyond the training period
        if (trainCountDown == 0) {
            estRatio = estRatio - 1;
        }
    }
    
    // update indices
    if (lastQrsIdx != 0) {
        lastQrsIdx--;
    }
    if (searchBackIdx != 0) {
        searchBackIdx--;
    }
    
    return wasQrsDetected || *detectedBySearchback;
}

/*=========================================================================
 * Event triggered for an arriving sample
 *=======================================================================*/
void onNewSample(short sample)
{
    bool detectedBySearchback;
    
    // preprocess sample
    filtSig[ci] = preprocessSample(sample);
    
    // update detection buffer
    detb(ci) = filtSig[ci];
    
    // detect qrs and
    if (detectQrs(&detectedBySearchback) && qrsHist != NULL) {
        qrsHist[qrsCount++] = ci + peakIdx + 1;
    }
    
    // update qrs2 history
    if (detectedBySearchback && qrs2Hist != NULL) {
        qrs2Hist[qrsCount2++] = ci + peakIdx + 1;
    }
    
    // update main threshold history
    if (thr1Hist != NULL) {
        thr1Hist[ci] = sigThreshold;
    }
    
    // update main threshold history
    if (thr2Hist != NULL) {
        thr2Hist[ci] = sigThreshold >> 1;
    }
    
    // update secondary threshold history
    if (rrintHist != NULL) {
        rrintHist[ci] = rrIntMean;
    }
    
    // increment global index
    ci++;
}

/*=========================================================================
 * Design the preprocessing filters 
 *=======================================================================*/
void designPreprocessingFilters()
{
    int maxbits;
    
    // low-pass and second-order backward difference
    lpfWin = ((int)round(sampFreq / mainsFreq) << 1) + 1;
    
    // windowed-maximum
    wmfWin = sampFreq << 1;
    mafWinLog2 = (int)ceil(log2(wmfWin));
    
    // automatic gain control
    agcWin = (int)(sampFreq * (double)ACG_LEN_SEC);
    maxbits = (int)((sizeof(int) << 3) - 1) / mbdOrder;
    maxbits = (int)min((sizeof(int) << 2) - 1, maxbits);
    agcMax = (1 << maxbits) - 1;
    agcMin = 8;
    
    // multiplication of absolute backward differences
    mbdWin = mbdOrder;
    
    // moving-average
    mafWin = (int)(sampFreq * (double)MAF_LEN_SEC);
    
    // calculate overall preprocessing delay
    delay = (lpfWin - 1) / 2.0 +    // LPF
            0 +                     // WMF
            agcWin +                // AGC
            (mbdWin - 1) / 2.0 +    // MBD
            (mafWin - 1) / 2.0;     // MAF
}

/*=========================================================================
 * Initialize detection variables 
 *=======================================================================*/
void initDetectionVariables()
{
    ci = 0;
    
    // levels and threshold
    sigThreshold = 0;
    signalLevel = 0;
    noiseLevel = 0;
    estRatio = 3;
    
    // peak
    peakIdx = 0;
    peakAmp = 0;
    
    // QRS
    qrsCount = 0;
    qrsCount2 = 0;
    lastQrsIdx = 0;
    searchBackIdx = 0;
    
    // RR interval
    rrIntMean = sampFreq;
    rrIntLow = rrIntMean >> 1;
    rrIntHigh = rrIntMean + rrIntLow;
    rrIntMiss = rrIntHigh + (rrIntMean >> 2);
    rrIntMin = sampFreq / 5;
    rrIntMax = sampFreq << 1;
    
    // flags
    isSignalRising = false;
    
    // periods
    trainCountDown = sampFreq << 1;
    qrsHalfLength = sampFreq / 20;
    twaveTolerance = (sampFreq * 9) / 25;
    refractoryPeriod = sampFreq / 5;
    
    // create buffer for the filtered signal
    detBufLen = 1 << (int)ceil(log2(2 * sampFreq));
    detBuf = (int *)mxMalloc(detBufLen * sizeof(int));
}

/*=========================================================================
 * Check correct number of arguments and their types 
 *=======================================================================*/
void checkArgs( int nlhs, mxArray *plhs[],
                int nrhs, const mxArray *prhs[])
{
    /* check for proper number of input arguments */
    if (nrhs < MIN_INPUTS || nrhs > MAX_INPUTS) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_detect_qrs_int:nrhs",
                "%d input(s) required.", MIN_INPUTS);
    }
    /* check for proper number of output arguments */
    if (nlhs < MIN_OUTPUTS || nlhs > MAX_OUTPUTS) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_detect_qrs_int:nlhs",
                "%d output(s) required.", MIN_OUTPUTS);
    }
    /* make sure the first input argument is a vector */
    if (mxGetM(prhs[0]) != 1 && mxGetN(prhs[0]) != 1) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_detect_qrs_int:notVector",
                "First input must be a Vector.");
    }
    /* make sure the first input argument is of type int16 */
    if (!mxIsInt16(prhs[0]) || mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_detect_qrs_int:notInt16",
                "First input must be of type Int16.");
    }
    /* make sure the remaining arguments are all scalars */
    for (mwSize i = 1; i < nrhs; i++) {
        if (mxGetNumberOfElements(prhs[i]) != 1) {
            mexErrMsgIdAndTxt(
                    "EcgToolbox:c_detect_qrs_int:notScalar",
                    "Input #%d must be a scalar.", i+1);
        }
    }
}

/*=========================================================================
 * Hanlde input arguments 
 *=======================================================================*/
void handleInputs( int nrhs, const mxArray *prhs[],
                   mwSize *nrows, mwSize *ncols)
{
    /* get a pointer to the data in the input vector  */
    inputSig = (short *)mxGetData(prhs[0]);
    
    /* get the dimensions of the input vector */
    *nrows = (mwSize)mxGetM(prhs[0]);
    *ncols = (mwSize)mxGetN(prhs[0]);
    
    /* get the sampling frequency  */
    sampFreq = (int)mxGetScalar(prhs[1]);
    
    /* get the mains frequency  */
    if (nrhs > 2) {
        mainsFreq = (int)mxGetScalar(prhs[2]);
    }
    else if (sampFreq % 50 == 0) {
        mainsFreq = 50; // guess based on the multiplicity of Fs by 50 Hz
    }
    else if (sampFreq % 60 == 0) {
        mainsFreq = 60; // guess based on the multiplicity of Fs by 60 Hz
    }
    else {
        mainsFreq = DEF_MAINS_FREQ;
    }
    
    /* get the MBD filter order  */
    if (nrhs > 3) {
        mbdOrder = (int)mxGetScalar(prhs[3]);
    }
    else {
        mbdOrder = DEF_MBD_ORDER;
    }
    
    /* make sure the mains frequency is within pre-defined limits */
    if (mainsFreq < MIN_MAINS_FREQ || mainsFreq > MAX_MAINS_FREQ) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_detect_qrs_int:badMainsFreq",
                "Mains frequency must be between %d and %d.",
                MIN_MAINS_FREQ, MAX_MAINS_FREQ);
    }
    /* make sure the sampling frequency is within pre-defined limits */
    if (sampFreq < mainsFreq || sampFreq > mainsFreq * LPF_BUFLEN) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_detect_qrs_int:badSampFreq",
                "Sampling frequency must be between %d and %d.",
                mainsFreq, mainsFreq * LPF_BUFLEN);
    }
    /* make sure the MBD filter order is within pre-defined limits */
    if (mbdOrder < MIN_MBD_ORDER || mbdOrder > MAX_MBD_ORDER) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_detect_qrs_int:badMbdOrder",
                "MBD filter order must be between %d and %d.",
                MIN_MBD_ORDER, MAX_MBD_ORDER);
    }
}

/*=========================================================================
 * Hanlde output arguments 
 *=======================================================================*/
void handleOutputs( int nlhs, mxArray *plhs[],
                    mwSize nrows, mwSize ncols)
{
    int *outVectors[MAX_OUTPUTS-1] = {NULL};
    
    /* create the output vectors */
    for (mwSize i = 0; i < min(nlhs,MAX_OUTPUTS-1); i++) {
        plhs[i] = mxCreateNumericMatrix(nrows,ncols,mxINT32_CLASS,mxREAL);
        outVectors[i] = (int *)mxGetData(plhs[i]);
    }
    
    /* get pointers to the output vectors */
    filtSig = outVectors[0];
    qrsHist = outVectors[1];
    qrs2Hist = outVectors[2];
    thr1Hist = outVectors[3];
    thr2Hist = outVectors[4];
    rrintHist = outVectors[5];
}

/*=========================================================================
 * Adjust the size of an mxArray 
 *=======================================================================*/
void adjustSize( int *vector, mxArray *mxarray,
                 mwSize nrows, mwSize ncols, mwSize newlen)
{
    if (vector != NULL && mxarray != NULL) {
        if (ncols == 1) {
            mxSetM(mxarray, newlen);
        }
        else {
            mxSetN(mxarray, newlen);
        }
        mxSetPr(mxarray, mxRealloc(vector, newlen * sizeof(int)));
    }
}

/*=========================================================================
 * Finalization routine 
 *=======================================================================*/
void finalize( int nlhs, mxArray *plhs[],
               mwSize nrows, mwSize ncols)
{
    /* adjust the size of qrs history output */
    if (nlhs > 1) {
        adjustSize(qrsHist,plhs[1],nrows,ncols,qrsCount);
    }
    
    /* adjust the size of qrs2 history output */
    if (nlhs > 2) {
        adjustSize(qrs2Hist,plhs[2],nrows,ncols,qrsCount2);
    }
    
    /* create the last output (preprocessing delay) */
    if (nlhs > 6) {
        plhs[6] = mxCreateDoubleScalar(delay);
    }
    
    /* deallocate memory */
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
    mwSize inLen;
    double time;
    
    /* check argument correctness */
    checkArgs(nlhs, plhs, nrhs, prhs);
    
    /* handle input arguments */
    handleInputs(nrhs, prhs, &nrows, &ncols);
    
    /* handle output arguments */
    handleOutputs(nlhs, plhs, nrows, ncols);
    
    /* calculate the length of the signal */
    inLen = max(nrows,ncols);
    
    /* design the preprocessing filters */
    designPreprocessingFilters();
    
    /* initialize detection variables */
    initDetectionVariables();
    
    /* process one input sample at a time */
    tic();
    for (mwSize i = 0; i < inLen; i++) {
        onNewSample(inputSig[i]);
    }
    time = toc();
    mexPrintf("Total processing time: %.2f ms\n", 1000*time);
    mexPrintf("Average time per sample: %.2f ns\n", 1000000000*time/inLen);
    
    /* perform some final adjustments */
    finalize(nlhs, plhs, nrows, ncols);
}
