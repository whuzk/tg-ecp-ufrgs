/*=========================================================================
 * c_prod_detect_qrs_int.c
 * 
 *  Title: real-time detection of QRS complex using integer arithmetic
 *         (production code)
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
 *      1. location of detected QRS (in samples)
 *      2. history of R-R intervals (in samples)
 *      3. overall preprocessing delay (in samples)
 *
 *  *required
 *
 *=======================================================================*/
#include <math.h>
#include "mex.h"

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
#define MIN_OUTPUTS     0
#define MAX_OUTPUTS     3

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
#define MBD_BUFLEN      8       // 2^nextpow2(MAX_MBD_ORDER)
#define MAF_BUFLEN      1024    // 2^nextpow2(LPF_BUFLEN/2 * MAX_MAINS_FREQ * MAF_LEN_SEC)

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
static int *qrsHist;        // indices of detected QRS complex
static int *rrintHist;      // history of average RR interval
static double delay;        // overall preprocessing delay

/*=========================================================================
 * Preprocessing variables
 *=======================================================================*/
static mwSize lpfWin;                       // LPF window length
static mwSize wmfWin;                       // WMF window length
static mwSize agcWin;                       // AGC window length
static mwSize mbdWin;                       // MBD window length
static mwSize mafWin;                       // MAF window length
static int agcMin;                          // AGC lower limit
static int agcMax;                          // AGC upper limit
static int mbdGainLog2;                     // log2 of MDB gain
static short lpfBuf[LPF_BUFLEN] = {0};      // LPF buffer
static int wmfBuf[WMF_BUFLEN] = {0};        // WMF buffer
static unsigned int wmfAux[WMF_BUFLEN];     // WMF aux buffer
static int agcBuf[ACG_BUFLEN] = {0};        // AGC buffer
static int mbdBuf[MBD_BUFLEN] = {1};        // MBD buffer
static int mafBuf[MAF_BUFLEN] = {0};        // MAF buffer
static int *detBuf;                         // detection buffer
static mwSize detBufLen;                    // detection buffer length
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
#define detb(I) (detBuf[(ci+(I))&(detBufLen-1)])

/*=========================================================================
 * QRS Detection variables
 *=======================================================================*/
static int sigThreshold;        // signal threshold
static int signalLevel;         // signal level
static int noiseLevel;          // noise level
static int estRatio;            // ratio of signal/noise level estimation
static int peakAmp;             // currently detected peak amplitude
static mwSize rrIntMean;        // running average of RR intervals
static mwSize rrIntMean2;       // secondary running average of RR intervals
static mwSize rrIntMiss;        // interval for qrs to be assumed as missed
static mwSize rrIntLow;         // lower bound for acceptance of new RR
static mwSize rrIntHigh;        // upper bound for acceptance of new RR
static mwSize rrIntMin;         // lowest bound applied to estimated RR
static mwSize trainCountDown;   // training period countdown
static mwSize qrsHalfLength;    // half the length of a QRS complex
static mwSize twaveTolerance;   // tolerance for T-wave detection
static mwSize refractoryPeriod; // length of refractory period for QRS detection
static mwSize qrsCount;         // count of QRS complex in main QRS buffer
static mwSize lastQrsIdx;       // index of last detected QRS complex
static mwSize searchBackIdx;    // index of searchback starting point
static mwSize peakIdx;          // currently detected peak index
static bool isSignalRising;     // flag to indicate a rise in the signal

/*=========================================================================
 * Low-pass filter and second-order backward difference
 *=======================================================================*/
int lpf(short sample)
{
    // update filter output: y[n] = x[n] - 2*x[n-M/2] + x[n-M]
    int y0 = sample - (lpfb(-(lpfWin >> 1)) << 1) + lpfb(-lpfWin);
    
    // update filter memory
    lpfb(0) = sample;
    
    // compute result
    return y0;// >> 2;
}

/*=========================================================================
 * Windowed-maximum
 *=======================================================================*/
int wmf(int sample)
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
int agc(int sample)
{
    int gain = wmf(sample);
    int y0 = agcb(-agcWin);
    
    if (gain >= agcMin) {
        // update filter output: y[n] = x[n-M] * maxGain / gain
        y0 *= agcMax / gain;
    }
    
    // update filter memory
    agcb(0) = sample;
    
    // compute result
    return min(y0, agcMax-1);
}

/*=========================================================================
 * Multiplication of absolute backward differences
 *=======================================================================*/
int mdb(int sample)
{
    long long y0 = sample;
    
    // update filter output: y[n] = x[n] * ... * x[n-M+1]
    for (mwSize k = 1; k < mbdWin; k++) {
        y0 *= mbdb(-k);
    }
    
    // update filter memory
    mbdb(0) = sample;
    
    // compute result
    return (int)(y0 >> mbdGainLog2);
}

/*=========================================================================
 * Moving-average 
 *=======================================================================*/
int maf(int sample)
{
    static int y0 = 0;
    
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
int preprocessSample(short sample)
{
    int sample2 = lpf(sample);
    sample2 = abs(sample2);
    sample2 = agc(sample2);
    sample2 = mdb(sample2);
    return maf(sample2);
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
    int y, newy;
    
    y = detb(start);
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
            int slope1 = maxdiff(start1, qrsHalfLength);
            int slope2 = maxdiff(start2, qrsHalfLength);

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
bool detectPeak(int newAmp)
{
    static mwSize lastPeakIdx = 0;      // index of the last peak
    static int lastPeakAmp = 0;         // amplitude of the last peak
    
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
    if (peakAmp >= sigThreshold) {
        signalLevel = estimate(signalLevel, estRatio, peakAmp);
        
        if (trainCountDown > 0) {
            lastQrsIdx = peakIdx;
            searchBackIdx = peakIdx;
            return false;
        }
        else if (peakIdx - lastQrsIdx > refractoryPeriod) {
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
    if ((!isSignalRising) && trainCountDown == 0 &&
            0 >= searchBackIdx + rrIntMiss) {
        mwSize begin = 0 - rrIntMean + 1;
        
        if (isvalidindex(begin)) {
            // search back and locate the max in this interval
            peakIdx = begin + findmax(begin, rrIntMean);
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
                //signalLevel = 0.5*signalLevel;
                //noiseLevel = 0.5*noiseLevel;
                // postpone searchback
                searchBackIdx += rrIntMean;
                return false;
            }
        }
        else return false;
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
        rrIntMean = max(rrIntMin, estimate(rrIntMean, 3, newRR));
    }
    rrIntMean2 = estimate(rrIntMean2, 3, newRR);
    
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
bool detectQrs(int sample)
{
    bool wasPeakDetected;   // flag to indicate that a peak was detected
    bool wasQrsDetected;    // flag to indicate that a qrs was detected
    
    // update detection buffer
    detb(0) = sample;
    
    // detect peak
    wasPeakDetected = detectPeak(sample);
    
    // detect qrs
    wasQrsDetected = wasPeakDetected && evaluatePeak();
    
    // search back
    wasQrsDetected = wasQrsDetected || searchBack();
    
    // update threshold
    sigThreshold = noiseLevel + (abs(signalLevel - noiseLevel) >> 2);
    
    // update qrs info
    if (wasQrsDetected) {
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
    lastQrsIdx--;
    searchBackIdx--;
    
    return wasQrsDetected;
}

/*=========================================================================
 * Event triggered for an arriving sample
 *=======================================================================*/
void onNewSample(short sample)
{
    // preprocess sample
    int sample2 = preprocessSample(sample);
    
    // detect qrs and update outputs
    if (detectQrs(sample2)) {
        if (qrsHist != NULL) {
            qrsHist[qrsCount] = ci + peakIdx + 1;
        }
        if (rrintHist != NULL) {
            rrintHist[qrsCount] = rrIntMean2;
        }
        qrsCount++;
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
    lpfWin = n << 1;
    
    // windowed-maximum
    wmfWin = sampFreq << 1;
    
    // automatic gain control
    agcWin = (int)(sampFreq * (double)ACG_LEN_SEC);
    agcMax = 1 << 15;
    agcMin = 1 << 3;
    
    // multiplication of absolute backward differences
    mbdWin = mbdOrder;
    mbdGainLog2 = (mbdOrder - 1) * 15;
    
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
    lastQrsIdx = 0;
    searchBackIdx = 0;
    
    // RR interval
    rrIntMean = sampFreq;
    rrIntMean2 = sampFreq;
    rrIntLow = rrIntMean >> 1;
    rrIntHigh = rrIntMean + rrIntLow;
    rrIntMiss = rrIntHigh + (rrIntMean >> 2);
    rrIntMin = sampFreq / 5;
    
    // flags
    isSignalRising = false;
    
    // periods
    trainCountDown = sampFreq << 1;
    qrsHalfLength = sampFreq / 20;
    twaveTolerance = (sampFreq * 9) / 25;
    refractoryPeriod = sampFreq / 5;
    
    // create buffer for the filtered signal
    detBufLen = 1 << (int)ceil(log2(sampFreq << 1));
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
    if (sampFreq < mainsFreq || sampFreq > mainsFreq * LPF_BUFLEN/2) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_detect_qrs_int:badSampFreq",
                "Sampling frequency must be between %d and %d.",
                mainsFreq, mainsFreq * LPF_BUFLEN/2);
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
    /* create output for the QRS history */
    if (nlhs > 0) {
        plhs[0] = mxCreateNumericMatrix(nrows,ncols,mxINT32_CLASS,mxREAL);
        qrsHist = (int *)mxGetData(plhs[0]);
    }
    
    /* create output for the RR interval history */
    if (nlhs > 1) {
        plhs[1] = mxCreateNumericMatrix(nrows,ncols,mxINT32_CLASS,mxREAL);
        rrintHist = (int *)mxGetData(plhs[1]);
    }
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
    /* adjust the size of QRS history output */
    if (nlhs > 0) {
        adjustSize(qrsHist,plhs[0],nrows,ncols,qrsCount);
    }
    
    /* adjust the size of RR interval history output */
    if (nlhs > 1) {
        adjustSize(rrintHist,plhs[1],nrows,ncols,qrsCount);
    }
    
    /* create the last output (preprocessing delay) */
    if (nlhs > MAX_OUTPUTS-1) {
        plhs[MAX_OUTPUTS-1] = mxCreateDoubleScalar(delay);
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
    for (mwSize i = 0; i < inLen; i++) {
        onNewSample(inputSig[i]);
    }
    
    /* perform some final adjustments */
    finalize(nlhs, plhs, nrows, ncols);
}
