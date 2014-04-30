/*=========================================================================
 * c_detect_qrs_double.c
 * 
 *  Title: detect QRS complex using floating-point arithmetic
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
#include "c_ecg_utils.h"
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
#define MAX_MBD_ORDER   4
#define DEF_MBD_ORDER   3

#define ACG_LEN_SEC     0.05
#define MAF_LEN_SEC     0.10

#define LPF_BUFLEN      256     // supports sampling frequencies of up to (LPF_BUFLEN/2*mainsFreq) Hz
#define WMF_BUFLEN      512     // 2^nextpow2(MAX_WMF_LIST_SIZE) [determined experimentaly]
#define ACG_BUFLEN      512     // 2^nextpow2(LPF_BUFLEN/2 * MAX_MAINS_FREQ * ACG_LEN_SEC)
#define MBD_BUFLEN      8       // 2^nextpow2(LPF_BUFLEN/2 * MAX_MBD_ORDER)
#define MAF_BUFLEN      1024    // 2^nextpow2(LPF_BUFLEN/2 * MAX_MAINS_FREQ * MAF_LEN_SEC)

/*=========================================================================
 * Input variables
 *=======================================================================*/
static double *inputSig;    // input signal
static int sampFreq;        // sampling frequency
static int mainsFreq;       // mains frequency
static int mbdOrder;        // MBD filter order

/*=========================================================================
 * Output variables
 *=======================================================================*/
static double *filtSig;     // filtered signal
static double *qrsHist;     // indices of detected QRS complex
static double *qrs2Hist;    // indices of QRS detected by searchback
static double *thr1Hist;    // history of signal threshold
static double *thr2Hist;    // history of secondary signal threshold
static double *rrintHist;   // history of average RR interval
static double delay;        // overall preprocessing delay

/*=========================================================================
 * Preprocessing variables
 *=======================================================================*/
static mwSize lpfWin;                       // LPF window length
static mwSize wmfWin;                       // WMF window length
static mwSize agcWin;                       // AGC window length
static mwSize mbdWin;                       // MBD window length
static mwSize mafWin;                       // MAF window length
static double agcMin;                       // AGC lower limit
static double agcMax;                       // AGC upper limit
static double mbdGain;                      // MBD gain
static double lpfBuf[LPF_BUFLEN] = {0};     // LPF buffer
static double wmfBuf[WMF_BUFLEN] = {0};     // WMF buffer
static unsigned int wmfAux[WMF_BUFLEN];     // WMF aux buffer
static double agcBuf[ACG_BUFLEN] = {0};     // AGC buffer
static double mbdBuf[MBD_BUFLEN] = {1};     // MBD buffer
static double mafBuf[MAF_BUFLEN] = {0};     // MAF buffer
static unsigned int ci;                     // current sample index

/*=========================================================================
 * Fast lookup for filter buffers
 *=======================================================================*/
#define lpfb(I) (lpfBuf[(ci+(I))&(LPF_BUFLEN-1)])
#define wmfb(I) (wmfBuf[(I)&(WMF_BUFLEN-1)])
#define wmfa(I) (wmfAux[(I)&(WMF_BUFLEN-1)])
#define agcb(I) (agcBuf[(ci+(I))&(ACG_BUFLEN-1)])
#define mbdb(I) (mbdBuf[(ci+(I))&(MBD_BUFLEN-1)])
#define mafb(I) (mafBuf[(ci+(I))&(MAF_BUFLEN-1)])

/*=========================================================================
 * QRS Detection variables
 *=======================================================================*/
static double sigThreshold;     // signal threshold
static double signalLevel;      // signal level
static double noiseLevel;       // noise level
static double estRatio;         // ratio of signal/noise level estimation
static double peakAmp;          // currently detected peak amplitude
static double rrIntMean;        // running average of RR intervals
static double rrIntMiss;        // interval for qrs to be assumed as missed
static double rrIntLow;         // lower bound for acceptance of new RR
static double rrIntHigh;        // upper bound for acceptance of new RR
static double rrIntMin;         // lowest bound applied to estimated RR
static mwSize trainingPeriod;   // length of training period
static mwSize qrsHalfLength;    // half the length of a QRS complex
static mwSize twaveTolerance;   // tolerance for T-wave detection
static mwSize refractoryPeriod; // length of refractory period for QRS detection
static mwSize qrsCount;         // count of QRS complex in main QRS buffer
static mwSize qrsCount2;        // current position in second QRS buffer
static unsigned int lastQrsIdx;     // index of last detected QRS complex
static unsigned int searchBackIdx;  // index of searchback starting point
static unsigned int peakIdx;        // currently detected peak index
static bool isSignalRising;     // flag to indicate a rise in the signal

/*=========================================================================
 * Low-pass filter and second-order backward difference
 *=======================================================================*/
double lpf(double sample)
{
    // update filter output: y[n] = x[n] - 2*x[n-M/2] + x[n-M]
    double y0 = sample - 2 * lpfb(-lpfWin/2) + lpfb(-lpfWin);
    
    // update filter memory
    lpfb(0) = sample;
    
    // compute result
    return y0;// / 4
}

/*=========================================================================
 * Windowed-maximum
 *=======================================================================*/
double wmf(double sample)
{
    static mwSize first = 0;
    static mwSize count = 0;
    mwSize j,k;
    
    // search for the first element greater than the current sample
    j = count;
    k = first + count - 1;
    while (j > 0 && sample >= wmfb(k)) {
        j--;
        k--;
    }
    // put the sample next to the element found and adjust the length
    wmfb(k + 1) = sample;
    wmfa(k + 1) = ci;
    count = j + 1;
    // check if the first in line has gone out of the window length
    if (count > wmfWin || wmfa(first) == ci - wmfWin) {
        first++;
        count--;
    }
    // return the max in the window
    return wmfb(first);
}

/*=========================================================================
 * Automatic gain control
 *=======================================================================*/
double agc(double sample)
{
    double gain = wmf(sample);
    double y0 = agcb(-agcWin);
    
    if (gain >= agcMin) {
        // update filter output: y[n] = x[n-M] * maxGain / gain
        y0 *= agcMax / gain;
    }
    
    // update filter memory
    agcb(0) = sample;
    
    // compute result
    return fmin(y0, agcMax-1);
}

/*=========================================================================
 * Multiplication of absolute backward differences
 *=======================================================================*/
double mdb(double sample)
{
    double y0 = sample;
    
    // update filter output: y[n] = x[n] * ... * x[n-M+1]
    for (mwSize k = 1; k < mbdWin; k++) {
        y0 *= mbdb(-k);
    }
    
    // update filter memory
    mbdb(0) = sample;
    
    // compute result
    return y0 / mbdGain;
}

/*=========================================================================
 * Moving-average 
 *=======================================================================*/
double maf(double sample)
{
    static double y0 = 0;
    
    // update filter output: y[n] = y[n-1] + x[n] - x[n-M]
    y0 += sample - mafb(-mafWin);
    
    // update filter memory
    mafb(0) = sample;
    
    // compute result
    return y0;// / mafWin;
}

/*=========================================================================
 * Preprocess one input sample 
 *=======================================================================*/
double preprocessSample(double sample)
{
    sample = lpf(sample);
    sample = fabs(sample);
    sample = agc(sample);
    sample = mdb(sample);
    return maf(sample);
}

/*=========================================================================
 * Check for T wave 
 *=======================================================================*/
bool istwave(unsigned int candQrs, unsigned int lastQrs)
{
    // check if the canditate QRS occurs near the previous one
    if ((int)(candQrs - lastQrs) < twaveTolerance) {
        // calculate starting indices
        unsigned int start1 = candQrs - qrsHalfLength + 1;
        unsigned int start2 = lastQrs - qrsHalfLength + 1;
        
        // max slope of waveforms
        double slope1 = maxdiff(filtSig + start1, qrsHalfLength);
        double slope2 = maxdiff(filtSig + start2, qrsHalfLength);
        
        // check condition for T wave
        return (slope1 < 0.5 * slope2);
    }
    else return false;
}

/*=========================================================================
 * Peak detection 
 *=======================================================================*/
bool detectPeak(double newAmp)
{
    static unsigned int lastPeakIdx = 0;    // index of the last peak
    static double lastPeakAmp = 0.0;        // amplitude of the last peak
    
    // check if the new amplitude is greater than that of the last peak
    if (newAmp > lastPeakAmp) {
        // signalize beginning or continuation of positive slope
        isSignalRising = true;
        // update peak info
        lastPeakAmp = newAmp;
        lastPeakIdx = ci;
        return false;
    }
    else if (!isSignalRising) {
        // update current amplitude
        lastPeakAmp = newAmp;
        return false;
    }
    else if (newAmp < 0.5 * lastPeakAmp) {
        // report the new peak amplitude and location
        peakAmp = lastPeakAmp;
        peakIdx = lastPeakIdx;
        // reset state of positive slope
        isSignalRising = false;
        lastPeakAmp = newAmp;
        // signalize detection
        return true;
    }
    else return false;
}

/*=========================================================================
 * Peak evaluation 
 *=======================================================================*/
bool evaluatePeak()
{
    if (peakAmp >= sigThreshold) {
        signalLevel = estimate(signalLevel, estRatio, peakAmp);
        
        if (ci <= (unsigned int)trainingPeriod) {
            lastQrsIdx = peakIdx;
            searchBackIdx = peakIdx;
            return false;
        }
        else if ((int)(peakIdx - lastQrsIdx) > refractoryPeriod) {
            return !istwave(peakIdx, lastQrsIdx);
        }
        else return false;
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
    if ((!isSignalRising) && ci > (unsigned int)trainingPeriod &&
            (int)(ci - searchBackIdx) >= (int)rrIntMiss) {
        // search back and locate the max in this interval
        mwSize len = (int)rrIntMean;
        unsigned int begin = ci - len + 1;
        peakIdx = begin + findmax(filtSig + begin, len);
        peakAmp = filtSig[peakIdx];
        
        // check if candidate peak is from qrs
        if (peakAmp < sigThreshold &&
                peakAmp >= 0.5 * sigThreshold &&
                !istwave(peakIdx, lastQrsIdx)) {
            // adjust signal level
            signalLevel = estimate(signalLevel, 0.25, peakAmp);
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
    mwSize newRR = peakIdx - lastQrsIdx;
    if (rrIntLow < newRR && newRR < rrIntHigh) {
        rrIntMean = fmax(rrIntMin, estimate(rrIntMean, 0.125, newRR));
    }
    
    // calculate RR missed limit
    rrIntLow = 0.5 * rrIntMean;
    rrIntHigh = 1.5 * rrIntMean;
    rrIntMiss = 1.75 * rrIntMean;
    
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
    wasPeakDetected = detectPeak(filtSig[ci]);
    
    // detect qrs
    wasQrsDetected = wasPeakDetected && evaluatePeak();
    
    // search back
    *detectedBySearchback = (!wasQrsDetected) && searchBack();
    
    // update threshold
    sigThreshold = noiseLevel + 0.25*fabs(signalLevel - noiseLevel);
    
    // update qrs info
    if (wasQrsDetected || *detectedBySearchback) {
        updateQrsInfo();
    }
    
    // decrease estimation ratio for times beyond the training period
    if (ci == trainingPeriod) {
        estRatio = estRatio / 2.0;
    }
    
    return wasQrsDetected || *detectedBySearchback;
}

/*=========================================================================
 * Event triggered for an arriving sample
 *=======================================================================*/
void onNewSample(double sample)
{
    bool detectedBySearchback;
    
    // preprocess sample
    filtSig[ci] = preprocessSample(sample);
    
    // detect qrs and update qrs history
    if (detectQrs(&detectedBySearchback) && qrsHist != NULL) {
        qrsHist[qrsCount++] = peakIdx + 1;
    }
    
    // update qrs2 history
    if (detectedBySearchback && qrs2Hist != NULL) {
        qrs2Hist[qrsCount2++] = peakIdx + 1;
    }
    
    // update main threshold history
    if (thr1Hist != NULL) {
        thr1Hist[ci] = sigThreshold;
    }
    
    // update secondary threshold history
    if (thr2Hist != NULL) {
        thr2Hist[ci] = 0.5*sigThreshold;
    }
    
    // update RR interval history
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
    int n;
    
    // low-pass and second-order backward difference
    n = (int)round(sampFreq / mainsFreq);
    lpfWin = 2 * n;
    
    // windowed-maximum
    wmfWin = 2 * sampFreq;
    
    // automatic gain control
    agcWin = (int)(sampFreq * (double)ACG_LEN_SEC);
    agcMax = pow(2, 15);
    agcMin = pow(2, 3);
    
    // multiplication of absolute backward differences
    mbdWin = mbdOrder;
    mbdGain = pow(2, (mbdOrder - 1) * 15);
    
    // moving-average
    mafWin = (int)(sampFreq * (double)MAF_LEN_SEC);
    
    // calculate overall preprocessing delay (LPF + WMF + AGC + MBD + MAF)
    delay = n + 0 + agcWin + (mbdWin - 1) / 2.0 + (mafWin - 1) / 2.0;
}

/*=========================================================================
 * Initialize detection variables 
 *=======================================================================*/
void initDetectionVariables()
{
    // levels and threshold
    sigThreshold = 0.0;
    signalLevel = 0.0;
    noiseLevel = 0.0;
    estRatio = 0.125;
    
    // peak
    peakIdx = 0;
    peakAmp = 0.0;
    
    // QRS
    qrsCount = 0;
    qrsCount2 = 0;
    lastQrsIdx = 0;
    searchBackIdx = 0;
    
    // RR interval
    rrIntMean = sampFreq;
    rrIntLow = 0.5 * rrIntMean;
    rrIntHigh = 1.5 * rrIntMean;
    rrIntMiss = 1.75 * rrIntMean;
    rrIntMin = 0.2 * sampFreq;
    
    // flags
    isSignalRising = false;
    
    // limits
    trainingPeriod = 2 * sampFreq;
    qrsHalfLength = (int)round(0.05 * sampFreq);
    twaveTolerance = (int)round(0.36 * sampFreq);
    refractoryPeriod = (int)round(0.20 * sampFreq);
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
                "EcgToolbox:c_detect_qrs_double:nrhs",
                "%d input(s) required.", MIN_INPUTS);
    }
    /* check for proper number of output arguments */
    if (nlhs < MIN_OUTPUTS || nlhs > MAX_OUTPUTS) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_detect_qrs_double:nlhs",
                "%d output(s) required.", MIN_OUTPUTS);
    }
    /* make sure the first input argument is a vector */
    if (mxGetM(prhs[0]) != 1 && mxGetN(prhs[0]) != 1) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_detect_qrs_double:notVector",
                "First input must be a Vector.");
    }
    /* make sure the first input argument is of type double */
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_detect_qrs_double:notDouble",
                "First input must be of type Double.");
    }
    /* make sure the remaining arguments are all scalars */
    for (mwSize i = 1; i < nrhs; i++) {
        if (mxGetNumberOfElements(prhs[i]) != 1) {
            mexErrMsgIdAndTxt(
                    "EcgToolbox:c_detect_qrs_double:notScalar",
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
    inputSig = mxGetPr(prhs[0]);
    
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
        mainsFreq = DEF_MAINS_FREQ;     // default mains frequency
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
                "EcgToolbox:c_detect_qrs_double:badMainsFreq",
                "Mains frequency must be between %d and %d.",
                MIN_MAINS_FREQ, MAX_MAINS_FREQ);
    }
    /* make sure the sampling frequency is within pre-defined limits */
    if (sampFreq < mainsFreq || sampFreq > mainsFreq * LPF_BUFLEN/2) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_detect_qrs_double:badSampFreq",
                "Sampling frequency must be between %d and %d.",
                mainsFreq, mainsFreq * LPF_BUFLEN/2);
    }
    /* make sure the MBD filter order is within pre-defined limits */
    if (mbdOrder < MIN_MBD_ORDER || mbdOrder > MAX_MBD_ORDER) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_detect_qrs_double:badMbdOrder",
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
    double *outVectors[MAX_OUTPUTS-1] = {NULL};
    
    /* create the output vectors */
    for (mwSize i = 0; i < min(nlhs,MAX_OUTPUTS-1); i++) {
        plhs[i] = mxCreateDoubleMatrix(nrows,ncols,mxREAL);
        outVectors[i] = mxGetPr(plhs[i]);
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
