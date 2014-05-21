/*=========================================================================
 * c_preprocess_detect.c
 * 
 *  Title: real-time detection of QRS complex, fiducial marks, etc.
 *  Author:     Diego Sogari
 *  Modified:   13/May/2014
 *
 *  Intputs:
 *      1. ECG signal*
 *      2. Derivative of ECG signal*
 *      3. Morphological derivative of ECG signal*
 *      4. Denoised ECG signal*
 *      5. Sampling frequency (in Hz)*
 *      6. Number of beats for the adaptive template construction
 *
 *  Outputs:
 *      1. List of R peak locations (in samples)*
 *      2. List of averaged R-R intervals (in samples)*
 *      3. List of fiducial mark locations (in samples)*
 *      4. List of extracted beats*
 *      5. The template of normal beat*
 *
 *  *required
 *
 *=======================================================================*/
#include <math.h>
#include "mex.h"
#include "c_mathutils.h"
#include "c_timeutils.h"

/*=========================================================================
 * Constants
 *=======================================================================*/
#define MIN_INPUTS      5       // minimum number of input arguments
#define MAX_INPUTS      6       // maximum number of input arguments
#define MIN_OUTPUTS     5       // minimum number of output arguments
#define MAX_OUTPUTS     5       // maximum number of output arguments

#define MIN_SAMP_FREQ   50*2    // minimum sampling frequency
#define MAX_SAMP_FREQ   60*40   // maximum sampling frequency

#define QRS_MEM_STEP    10000   // step size for increasing the memory
                                //      allocation of the outputs vectors
#define RRI_BUFLEN      8       // eight most recent RR intervals
#define SEARCH_L1_SEC   0.10    // search limit for R-wave points
#define SEARCH_L2_SEC   0.02    // search limit for Ronset and Roffset
#define SEARCH_L3_SEC   0.15    // search limit for T- and P-wave points
#define SEARCH_L4_SEC   0.08    // search limit for R onset and offset
#define NUM_BL_SAMPLES  5       // number of samples for baseline removal
#define NUM_FDP_ROWS    7       // number of rows in the FDP matrix
#define HALF_FRAME      0.6     // half the size of the beat frame
#define DEF_NUM_BEATS   30      // default number of beats for the adaptive
                                //      template construction

/*=========================================================================
 * Input and output variables
 *=======================================================================*/
static double *inputSig1;       // input signal 1
static double *inputSig2;       // input signal 2
static double *inputSig3;       // input signal 3
static double *inputSig4;       // input signal 4
static mwSize inputLen;         // input signal length
static double sampFreq;         // sampling frequency
static double *qrsHist;         // history of QRS location
static double *rrHist;          // history of averaged RR
static double *fdpList;         // list of fiducial points
static double *beatList;        // list of extracted beats
static double *template;        // template of normal beat
static mwSize numBeats;         // number of beats for adaptation
static mwSize totalQrsLen;      // total number of columns in the output
                                //      vectors

/*=========================================================================
 * Buffer variables
 *=======================================================================*/
static unsigned int ci;             // current sample index
static unsigned int bi;             // current beat index
static mwSize rr1Buf[RRI_BUFLEN];   // primary RR buffer
static mwSize rr2Buf[RRI_BUFLEN];   // secondary RR buffer
static int *tpkBuf;                 // buffer for the tompkins signal
static int *deBuf;                  // buffer for the derivative signal
static int *mdBuf;                  // buffer for the MD signal
static double *noisBuf;             // buffer for the denoised signal
static double *beatBuf;             // buffer for the current beat
static double *tempBuf;             // buffer for the current template
static double *auxBuf1;             // auxiliary buffer for math operations
static double *auxBuf2;             // second auxiliary buffer
static mwSize bufLen;               // length of the buffers
static mwSize frameSize;            // length of beat frame

/*=========================================================================
 * Fast lookup for buffers
 *=======================================================================*/
#define rr1b(I)     (rr1Buf[(I)&(RRI_BUFLEN-1)])
#define rr2b(I)     (rr2Buf[(I)&(RRI_BUFLEN-1)])
#define tpkb(I)     (tpkBuf[(ci+(I))&(bufLen-1)])
#define deb(I)      (deBuf[(ci+(I))&(bufLen-1)])
#define mdb(I)      (mdBuf[(ci+(I))&(bufLen-1)])
#define noisb(I)    (noisBuf[(ci+(I))&(bufLen-1)])
#define fdp(I,J)    (fdpList[(I)*NUM_FDP_ROWS+(J)])
#define beat(I,J)   (beatList[(I)*frameSize+(J)])
#define temp(I,J)   (tempList[(I)*frameSize+(J)])

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
static mwSize lastQrsIdx;       // location of last detected QRS
static mwSize searchBackIdx;    // location of searchback starting point
static mwSize peakIdx;          // currently detected peak location
static bool isSignalRising;     // flag to indicate a rise in the signal
static mwSize Ponset;           // current P-wave onset
static mwSize Ppeak;            // current P-wave peak
static mwSize Ronset;           // current R-wave onset
static mwSize Rpeak;            // current R-wave peak
static mwSize Roffset;          // current R-wave offset
static mwSize Tpeak;            // current T-wave peak
static mwSize Toffset;          // current T-wave offset
static mwSize searchL1;         // search limit 1
static mwSize searchL2;         // search limit 2
static mwSize searchL3;         // search limit 3
static mwSize searchL4;         // search limit 4
static double artThresh[2];     // thresolds for artifact detection
static double tempRatio;        // ratio of adaptation for template
static double rmsdVal;          // value of RMSD (beat and template)
static mwSize savedRRInterval;  // copy of current RR interval average

/*=========================================================================
 * Check if index is valid in detection buffer 
 *=======================================================================*/
bool isvalidindex(mwSize idx)
{
    return (idx <= 0 && 0 < bufLen + idx);
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
    mwSize i;
    
    for (i = 1; i < len; i++) {
        newd = tpkb(start + i) - tpkb(start + i - 1);
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
    mwSize i, pos = 0;
    int y, newy;
    
    y = tpkb(start);
    for (i = 1; i < len; i++) {
        newy = tpkb(start + i);
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
    if ((!isSignalRising || tpkb(0) < (sigThreshold >> 1)) &&
            searchBackIdx != 0 && 0 >= searchBackIdx + rrIntMiss) {
        // calculate starting point for the search
        mwSize begin = 0 - rrIntMean2 + 1;
        
        if (isvalidindex(begin)) {
            // search back and locate the max in this interval
            peakIdx = findmax(begin, rrIntMean2);
            peakAmp = tpkb(peakIdx);
            
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
        mwSize i;
        for (i = 0; i < RRI_BUFLEN; i++) {
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
bool detectQrs()
{
    bool peakDetected;      // flag to indicate that a peak was detected
    bool qrsDetected;       // flag to indicate that a qrs was detected
    
    // detect peak
    peakDetected = detectPeak(tpkb(0));
    
    // detect qrs
    qrsDetected = peakDetected && evaluatePeak();
    
    // search back
    qrsDetected = qrsDetected || searchBack();
    
    // update threshold
    sigThreshold = noiseLevel + (abs(signalLevel - noiseLevel) >> 2);
    
    // update qrs info
    if (qrsDetected) {
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
    
    return qrsDetected;
}

/*=========================================================================
 * Check that the indices are correct for detection of fiducial points
 *=======================================================================*/
bool check_indices(mwSize start, mwSize end, mwSize *inc)
{
    if (abs(end - start) <= 1 || !isvalidindex(start) || !isvalidindex(end)) {
        return false;
    }
    else if (start < end) {
        *inc = 1;
        return true;
    }
    else {
        *inc = -1;
        return true;
    }
}

/*=========================================================================
 * Perform the detection of fiducial points that are peaks
 *=======================================================================*/
mwSize search_peak_abs(mwSize start, mwSize end, mwSize def)
{
    mwSize i, imax, idx[2];
    int y, ymax, val[2];
    mwSize left, right, inc;
    
    if (!check_indices(start, end, &inc)) {
        return def;
    }
    
    val[0] = val[1] = 0;
    idx[0] = idx[1] = def;
    for (i = start + inc; i != end - inc; i += inc) {
        y = deb(i);
        if (deb(i - 1) < y && y >= deb(i + 1)) {
            if (y > val[0]) {
                idx[0] = i;
                val[0] = y;
            }
        }
        else if (deb(i - 1) > y && y <= deb(i + 1)) {
            if (y < val[1]) {
                idx[1] = i;
                val[1] = y;
            }
        }
    }
    
    if (idx[0] < idx[1]) {
        left = idx[0];
        right = idx[1];
    }
    else {
        left = idx[1];
        right = idx[0];
    }
    
    ymax = 0;
    imax = left;
    for (i = left + 1; i < right - 1; i++) {
        y = abs(mdb(i));
        if (abs(mdb(i - 1)) < y && y >= abs(mdb(i + 1))) {
            if (y > ymax) {
                ymax = y;
                imax = i;
            }
        }
    }
    
    return imax;
}

/*=========================================================================
 * Perform the detection of the first maximum peak on the MD signal
 *=======================================================================*/
mwSize search_first_mark_max(mwSize start, mwSize end, mwSize def)
{
    mwSize i, imax;
    int y, ymax;
    mwSize inc;
    bool found;
    
    if (!check_indices(start, end, &inc)) {
        return def;
    }
    
    ymax = 0;
    imax = start;
    found = false;
    i = start + inc;
    while (!found && i != end - inc) {
        y = mdb(i);
        if (mdb(i - 1) < y && y >= mdb(i + 1)) {
            found = true;
        }
        else {
            /*if (y > ymax) {
                ymax = y;
                imax = i;
            }*/
            i += inc;
        }
    }
    
    if (!found) {
        return imax;
    }
    else return i;
}

/*=========================================================================
 * Perform the detection of the first minimum peak on the MD signal
 *=======================================================================*/
mwSize search_first_mark_min(mwSize start, mwSize end, mwSize def)
{
    mwSize i, imin;
    int y, ymin;
    mwSize inc;
    bool found;
    
    if (!check_indices(start, end, &inc)) {
        return def;
    }
    
    ymin = 0;
    imin = start;
    found = false;
    i = start + inc;
    while (!found && i != end - inc) {
        y = mdb(i);
        if (mdb(i - 1) > y && y <= mdb(i + 1)) {
            found = true;
        }
        else {
            /*if (y < ymin) {
                ymin = y;
                imin = i;
            }*/
            i += inc;
        }
    }
    
    if (!found) {
        return imin;
    }
    else return i;
}

/*=========================================================================
 * Perform the detection of the best maximum peak on the absolute value of
 * the MD signal
 *=======================================================================*/
mwSize search_best_mark_abs(mwSize start, mwSize end, mwSize def)
{
    mwSize i, imax, imaxpk;
    int y, ymax, ymaxpk;
    mwSize inc;
    bool found;
    
    if (!check_indices(start, end, &inc)) {
        return def;
    }
    
    found = false;
    ymax = ymaxpk = 0;
    imax = imaxpk = start;
    for (i = start + inc; i != end - inc; i += inc) {
        y = abs(mdb(i));
        if (abs(mdb(i - 1)) < y && y >= abs(mdb(i + 1))) {
            if (y > ymaxpk) {
                ymaxpk = y;
                imaxpk = i;
                found = true;
            }
        }
        else if (y > ymax) {
            ymax = y;
            imax = i;
        }
    }
    
    if (!found) {
        return imax;
    }
    else return imaxpk;
}

/*=========================================================================
 * Perform the detection of the best maximum peak on the MD signal
 *=======================================================================*/
mwSize search_best_mark_max(mwSize start, mwSize end, mwSize def)
{
    mwSize i, imax, imaxpk;
    int y, ymax, ymaxpk;
    mwSize inc;
    bool found;
    
    if (!check_indices(start, end, &inc)) {
        return def;
    }
    
    found = false;
    ymax = ymaxpk = 0;
    imax = imaxpk = start;
    for (i = start + inc; i != end - inc; i += inc) {
        y = mdb(i);
        if (mdb(i - 1) < y && y >= mdb(i + 1)) {
            if (y > ymaxpk) {
                ymaxpk = y;
                imaxpk = i;
                found = true;
            }
        }
        else if (y > ymax) {
            ymax = y;
            imax = i;
        }
    }
    
    if (!found) {
        return imax;
    }
    else return imaxpk;
}

/*=========================================================================
 * Perform the detection of the best minimum peak on the MD signal
 *=======================================================================*/
mwSize search_best_mark_min(mwSize start, mwSize end, mwSize def)
{
    mwSize i, imin, iminpk;
    int y, ymin, yminpk;
    mwSize inc;
    bool found;
    
    if (!check_indices(start, end, &inc)) {
        return def;
    }
    
    found = false;
    ymin = yminpk = 0;
    imin = iminpk = start;
    for (i = start + inc; i != end - inc; i += inc) {
        y = mdb(i);
        if (mdb(i - 1) > y && y <= mdb(i + 1)) {
            if (y < yminpk) {
                yminpk = y;
                iminpk = i;
                found = true;
            }
        }
        else if (y < ymin) {
            ymin = y;
            imin = i;
        }
    }
    
    if (!found) {
        return imin;
    }
    else return iminpk;
}

/*=========================================================================
 * Special search
 *=======================================================================*/
mwSize special_search(mwSize test, mwSize def)
{
    int y1 = mdb(test);
    int y2 = mdb(def);
    if (SIGN(y1) != SIGN(y2) && abs(y1) >= abs(y2 >> 2)) {
        return test;
    }
    else return def;
}

/*=========================================================================
 * Perform the detection of fiducial marks
 *=======================================================================*/
void searchMarks(mwSize r, mwSize rr1, mwSize rr2)
{
    mwSize len, temp1, temp2;
    
    // R-wave
    Rpeak = search_peak_abs(r - searchL1, r + searchL1, r);
    if (mdb(Rpeak) > 0) {
        // inverted
        Ronset = search_first_mark_min(Rpeak - searchL2, Rpeak - searchL1, Rpeak);
        Roffset = search_first_mark_min(Rpeak + searchL2, Rpeak + searchL1, Rpeak);
        temp1 = search_first_mark_max(Ronset - 1, Ronset - searchL4, Ronset);
        temp2 = search_first_mark_max(Roffset + 1, Roffset + searchL4, Roffset);
    }
    else {
        // normal
        Ronset = search_first_mark_max(Rpeak - searchL2, Rpeak - searchL1, Rpeak);
        Roffset = search_first_mark_max(Rpeak + searchL2, Rpeak + searchL1, Rpeak);
        temp1 = search_first_mark_min(Ronset - 1, Ronset - searchL4, Ronset);
        temp2 = search_first_mark_min(Roffset + 1, Roffset + searchL4, Roffset);
    }
    Ronset = special_search(temp1, Ronset);
    Roffset = special_search(temp2, Roffset);
    
    // P-wave
    len = (rr1 + (rr1 << 1)) >> 3;
    Ppeak = search_peak_abs(Ronset - 1, Ronset - len, Ronset);
    if (mdb(Ppeak) > 0) {
        // inverted
        Ponset = search_best_mark_min(Ppeak - 1, Ppeak - searchL3, Ppeak);
    }
    else {
        // normal
        Ponset = search_best_mark_max(Ppeak - 1, Ppeak - searchL3, Ppeak);
    }
    
    // T-wave
    len = rr2 >> 1;
    Tpeak = search_peak_abs(Roffset + 1, Roffset + len, Roffset);
    if (mdb(Tpeak) > 0) {
        // inverted
        Toffset = search_best_mark_min(Tpeak + 1, Tpeak + searchL3, Tpeak);
    }
    else {
        // normal
        Toffset = search_best_mark_max(Tpeak + 1, Tpeak + searchL3, Tpeak);
    }
}

/*=========================================================================
 * Detect beat
 *=======================================================================*/
bool detectBeat(bool qrsDetected)
{
    static mwSize countDown = -1;
    static mwSize rrInt = -1;
    static mwSize qrsIdx;
    bool beatDetected = false;
    mwSize rrInt2;
    
    // simulate QRS detection
    if ((qrsDetected && countDown > 0) || countDown == 0 || ci == inputLen - 1) {
        if (ci == inputLen - 1) {
            rrInt2 = -qrsIdx;
        }
        else if (countDown == 0) {
            rrInt2 = rrInt;
        }
        else {
            rrInt2 = peakIdx - qrsIdx;
        }
        // perform fiducial mark detection
        searchMarks(qrsIdx, rrInt, rrInt2);
        // signalize detection
        beatDetected = true;
    }
    
    if (qrsDetected) {
        // update rr
        if (rrInt == -1) {
            rrInt = (mwSize)sampFreq;
        }
        else rrInt = peakIdx - qrsIdx;
        // update qrs
        qrsIdx = peakIdx;
        // update countdown
        countDown = rrInt;
    }
    else if (countDown >= 0) {
        // update countdown
        countDown--;
    }
    
    qrsIdx--;
    
    return beatDetected;
}

/*=========================================================================
 * Perform the baseline removal procedure
 *=======================================================================*/
void removeBaseline()
{
    double x[2], y[2], p[2], mean;
    mwSize i, startIdx, endIdx;
    
    // calculate start and end points
    startIdx = max(1 - bufLen, Ponset - (NUM_BL_SAMPLES >> 1));
    endIdx = min(0, Toffset + (NUM_BL_SAMPLES >> 1));
    
    // average of first BL_NUM_PTS samples
    mean = 0.0;
    for (i = startIdx; i < startIdx + NUM_BL_SAMPLES; i++) {
        mean += noisb(i);
    }
    x[0] = 0.0;
    y[0] = mean / NUM_BL_SAMPLES;
    
    // average of last BL_NUM_PTS samples
    mean = 0.0;
    for (i = endIdx - NUM_BL_SAMPLES + 1; i <= endIdx; i++) {
        mean += noisb(i);
    }
    x[1] = endIdx - startIdx;
    y[1] = mean / NUM_BL_SAMPLES;
    
    // create polynomial fit
    first_order_fit(x, y, p);
    
    // calculate new amplitudes
    for (i = startIdx; i <= endIdx; i++) {
        noisb(i) -= p[0] * (i - startIdx) + p[1];
    }
}

/*=========================================================================
 * Frame one beat
 *=======================================================================*/
void frameBeat()
{
    mwSize startIdx, endIdx, absIdx1, absIdx2, len1, len2;
    mwSize r, half, beatsize, pad, cut1, cut2;
    
    // calculate padding and trimming sizes
    half = frameSize >> 1;
    beatsize = Toffset - Ponset + 1;
    r = Rpeak - Ponset + 1;
    pad = max(0, half - r + 1);
    cut1 = max(0, r - 1 - half);
    cut2 = max(0, beatsize - r - half);
    
    // update start and end indices
    startIdx = Ponset + cut1;
    endIdx = Toffset - cut2;
    if (startIdx > endIdx || !isvalidindex(startIdx) ||
            !isvalidindex(endIdx)) {
        return;
    }
    
    // calculate absolute indices in the buffer
    absIdx1 = (ci + startIdx) & (bufLen - 1);
    absIdx2 = (ci + endIdx) & (bufLen - 1);
    
    // frame the beat into the buffer
    memset(beatBuf, 0, frameSize * sizeof(double));
    if (absIdx1 <= absIdx2) {
        len1 = (absIdx2 - absIdx1 + 1) * sizeof(double);
        memcpy(beatBuf + pad, noisBuf + absIdx1, len1);
    }
    else {
        len1 = (bufLen - absIdx1) * sizeof(double);
        len2 = (absIdx2 + 1) * sizeof(double);
        memcpy(beatBuf + pad, noisBuf + absIdx1, len1);
        memcpy(beatBuf + pad + bufLen - absIdx1, noisBuf, len2);
    }
}

/*=========================================================================
 * Extract beat
 *=======================================================================*/
void extractBeat()
{
    // remove baseline
    removeBaseline();
    
    // frame the beat
    frameBeat();
}

/*=========================================================================
 * Detect artifact and update template beat
 *=======================================================================*/
bool detectArtifact()
{
    // save difference in auxBuf1 and calculate rmsd
    double_subtract(beatBuf, tempBuf, auxBuf1, frameSize);
    double_array_multiply(auxBuf1, auxBuf1, auxBuf2, frameSize);
    rmsdVal = sqrt(double_sum(auxBuf2, frameSize) / frameSize);
    
    if (bi == 0) {
        // first beat
        memcpy(tempBuf, beatBuf, frameSize * sizeof(double));
    }
    else if (bi < (unsigned int)numBeats || rmsdVal < artThresh[0]) {
        // very good beat
        double_scalar_multiply(auxBuf1, tempRatio, auxBuf1, frameSize);
        double_add(tempBuf, auxBuf1, tempBuf, frameSize);
        artThresh[0] += tempRatio * (2 * rmsdVal - artThresh[0]);
        artThresh[1] += tempRatio * (4 * rmsdVal - artThresh[1]);
    }
    else if (rmsdVal > artThresh[1]) {
        // artifact
        return true;
    }
    return false;
}

/*=========================================================================
 * Adjust the output dimensions
 *=======================================================================*/
void adjustOutputSizes()
{
    if (bi == totalQrsLen) {
        totalQrsLen += QRS_MEM_STEP;
        qrsHist = mxRealloc(qrsHist, totalQrsLen * sizeof(double));
        rrHist = mxRealloc(rrHist, totalQrsLen * sizeof(double));
        fdpList = mxRealloc(fdpList, NUM_FDP_ROWS * totalQrsLen * sizeof(double));
        beatList = mxRealloc(beatList, frameSize * totalQrsLen * sizeof(double));
    }
}

/*=========================================================================
 * Update output buffers
 *=======================================================================*/
void updateOutputs()
{
    mwSize center;

    // adjust the outputs
    adjustOutputSizes();

    // save R peak
    qrsHist[bi] = ci + Rpeak + 1;

    // save RR interval
    rrHist[bi] = savedRRInterval;

    // save fiducial marks
    center = (frameSize >> 1) + 1;
    fdp(bi,0) = max(1, center - (Rpeak - Ponset));
    fdp(bi,1) = max(1, center - (Rpeak - Ppeak));
    fdp(bi,2) = max(1, center - (Rpeak - Ronset));
    fdp(bi,3) = center;
    fdp(bi,4) = min(frameSize, center + (Roffset - Rpeak));
    fdp(bi,5) = min(frameSize, center + (Tpeak - Rpeak));
    fdp(bi,6) = min(frameSize, center + (Toffset - Rpeak));
    
    // save beat
    memcpy(&beat(bi, 0), beatBuf, frameSize * sizeof(double));

    // save template
    if (bi == numBeats - 1) {
        memcpy(template, tempBuf, frameSize * sizeof(double));
    }

    // increment beat index
    bi++;
}

/*=========================================================================
 * Event triggered for an arriving sample
 *=======================================================================*/
void onNewSample(double sample1, double sample2, double sample3, double sample4)
{
    bool qrsDetected;
    
    // update buffers
    tpkb(0) = (int)sample1;
    deb(0) = (int)sample2;
    mdb(0) = (int)sample3;
    noisb(0) = sample4;
    
    // detect qrs
    qrsDetected = detectQrs();
    
    // detect beat
    if (detectBeat(qrsDetected)) {
        // extract beat
        extractBeat();
        
        // update template and detect artifact
        //if (!detectArtifact()) {
        detectArtifact();
            // update outputs
            updateOutputs();
        //}
    }
    
    // save RR interval average
    if (qrsDetected) {
        savedRRInterval = rrIntMean1;
    }
    
    // increment global index
    ci++;
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
            "EcgToolbox:c_preprocess_detect:nrhs",
            "%d input(s) required and %d optional",
            MIN_INPUTS, MAX_INPUTS - MIN_INPUTS);
    }
    // check for proper number of output arguments
    if (nlhs < MIN_OUTPUTS || nlhs > MAX_OUTPUTS) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_preprocess_detect:nlhs",
            "%d output(s) required and %d optional",
            MIN_OUTPUTS, MAX_OUTPUTS - MIN_OUTPUTS);
    }
    // make sure all input arguments are of type double
    for (i = 0; i < nrhs; i++) {
        if (!mxIsDouble(prhs[i]) || mxIsComplex(prhs[i])) {
            mexErrMsgIdAndTxt(
                "EcgToolbox:c_preprocess_detect:notDouble",
                "Input #%d must be of type double.", i + 1);
        }
    }
    // make sure the first 4 input arguments are vectors
    for (i = 0; i < 4; i++) {
        if (mxGetM(prhs[i]) != 1 && mxGetN(prhs[i]) != 1) {
            mexErrMsgIdAndTxt(
                "EcgToolbox:c_preprocess_detect:notVector",
                "Input #%d must be a vector.", i + 1);
        }
    }
    // make sure the remaining input arguments are all scalars
    for (i = 4; i < nrhs; i++) {
        if (mxGetNumberOfElements(prhs[i]) != 1) {
            mexErrMsgIdAndTxt(
                "EcgToolbox:c_preprocess_detect:notScalar",
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
    
    // get the dimensions of the input vector
    *nrows = (mwSize)mxGetM(prhs[0]);
    *ncols = (mwSize)mxGetN(prhs[0]);
    
    // make sure the 1-4 input arguments have compatible dimensions
    for (i = 1; i < 4; i++) {
        if (mxGetM(prhs[i]) != *nrows || mxGetN(prhs[i]) != *ncols) {
            mexErrMsgIdAndTxt(
                "EcgToolbox:c_preprocess_detect:badDimensions",
                "Inputs 1-4 must have the same dimensions.");
        }
    }
    
    // get a pointer to the data in the input vectors
    inputSig1 = mxGetPr(prhs[0]);
    inputSig2 = mxGetPr(prhs[1]);
    inputSig3 = mxGetPr(prhs[2]);
    inputSig4 = mxGetPr(prhs[3]);
    
    // get the sampling frequency
    sampFreq = mxGetScalar(prhs[4]);
    
    // make sure the sampling frequency is within pre-defined limits
    if (sampFreq < MIN_SAMP_FREQ || sampFreq > MAX_SAMP_FREQ) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_preprocess_detect:badSampFreq",
            "Sampling frequency must be between %d and %d.",
            MIN_SAMP_FREQ, MAX_SAMP_FREQ);
    }
    
    // get the number of beats for template adaptation
    if (nrhs > 5) {
        numBeats = (int)mxGetScalar(prhs[5]);
    }
    else numBeats = DEF_NUM_BEATS;
}

/*=========================================================================
 * Hanlde output arguments 
 *=======================================================================*/
void handleOutputs(int nlhs, mxArray *plhs[])
{
    totalQrsLen = QRS_MEM_STEP;
    
    // get a pointer to the list of R peaks
    plhs[0] = mxCreateDoubleMatrix(totalQrsLen, 1, mxREAL);
    qrsHist = mxGetPr(plhs[0]);
    
    // get a pointer to the list of RR intervals
    plhs[1] = mxCreateDoubleMatrix(totalQrsLen, 1, mxREAL);
    rrHist = mxGetPr(plhs[1]);
    
    // get a pointer to the list of fiducial points
    plhs[2] = mxCreateDoubleMatrix(NUM_FDP_ROWS, totalQrsLen, mxREAL);
    fdpList = mxGetPr(plhs[2]);
    
    // get a pointer to the list of beats
    plhs[3] = mxCreateDoubleMatrix(frameSize, totalQrsLen, mxREAL);
    beatList = mxGetPr(plhs[3]);
    
    // get a pointer to the template vector
    plhs[4] = mxCreateDoubleMatrix(frameSize, 1, mxREAL);
    template = mxGetPr(plhs[4]);
}

/*=========================================================================
 * The initializarion routine
 *=======================================================================*/
void init()
{
    // initialize global indices
    ci = 0;
    bi = 0;
    
    // qrs detection variables
    sigThreshold = 0;
    signalLevel = 0;
    noiseLevel = 0;
    estRatio = 2;
    peakIdx = 0;
    peakAmp = 0;
    lastQrsIdx = 0;
    searchBackIdx = 0;
    isSignalRising = false;
    
    // qrs detection limits
    trainCountDown = (mwSize)(2 * sampFreq);
    qrsHalfLength = (mwSize)(0.10 * sampFreq);
    twaveTolerance = (mwSize)(0.36 * sampFreq);
    refractoryPeriod = (mwSize)(0.20 * sampFreq);
    
    // fiducial mark search limits
    searchL1 = (mwSize)(SEARCH_L1_SEC * sampFreq);
    searchL2 = (mwSize)(SEARCH_L2_SEC * sampFreq);
    searchL3 = (mwSize)(SEARCH_L3_SEC * sampFreq);
    searchL4 = (mwSize)(SEARCH_L4_SEC * sampFreq);
    
    // calculate the ratio for template adaptation
    tempRatio = 1.0 / numBeats;
    
    // create buffer for the filtered signal
    bufLen = 1 << (1 + ILOG2(4 * sampFreq - 1));
    tpkBuf = (int *)mxMalloc(bufLen * sizeof(int));
    deBuf = (int *)mxMalloc(bufLen * sizeof(int));
    mdBuf = (int *)mxMalloc(bufLen * sizeof(int));
    noisBuf = (double *)mxMalloc(bufLen * sizeof(double));
    beatBuf = (double *)mxMalloc(frameSize * sizeof(double));
    tempBuf = (double *)mxCalloc(frameSize, sizeof(double));
    auxBuf1 = (double *)mxMalloc(frameSize * sizeof(double));
    auxBuf2 = (double *)mxMalloc(frameSize * sizeof(double));
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
    for (i = 0; i < inputLen; i++) {
        onNewSample(inputSig1[i], inputSig2[i], inputSig3[i], inputSig4[i]);
    }
    
    // stop time counter
    time = toc();
    
    // display time statistics
    mexPrintf("Total processing time: %.2f ms\n", 1000 * time);
    mexPrintf("Average time per sample: %.2f ns\n", 1000000000 * time / inputLen);
}

/*=========================================================================
 * Adjust the size of an mxArray
 *=======================================================================*/
void adjustSize(mxArray *mxarray, double *vector, mwSize newm, mwSize newn)
{
    mxSetM(mxarray, newm);
    mxSetN(mxarray, newn);
    mxSetPr(mxarray, vector);
}

/*=========================================================================
 * The finalization routine 
 *=======================================================================*/
void finalize(int nlhs, mxArray *plhs[])
{
    // adjust the size of the list of R peaks
    qrsHist = mxRealloc(qrsHist, bi * sizeof(double));
    adjustSize(plhs[0], qrsHist, bi, 1);
    
    // adjust the size of the list of RR intervals
    rrHist = mxRealloc(rrHist, bi * sizeof(double));
    adjustSize(plhs[1], rrHist, bi, 1);
    
    // adjust the size of the list of fiducial points
    fdpList = mxRealloc(fdpList, NUM_FDP_ROWS * bi * sizeof(double));
    adjustSize(plhs[2], fdpList, NUM_FDP_ROWS, bi);
    
    // adjust the size of the list of Beats
    beatList = mxRealloc(beatList, frameSize * bi * sizeof(double));
    adjustSize(plhs[3], beatList, frameSize, bi);
    
    // deallocate memory
    mxFree(tpkBuf);
    mxFree(deBuf);
    mxFree(mdBuf);
    mxFree(noisBuf);
    mxFree(beatBuf);
    mxFree(tempBuf);
    mxFree(auxBuf1);
    mxFree(auxBuf2);
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
    
    // calculate the length of the signal
    inputLen = max(nrows,ncols);
    
    // calculate the frame size
    frameSize = 2 * (int)(HALF_FRAME * sampFreq) + 1;
    
    // handle output arguments
    handleOutputs(nlhs, plhs);
    
    // make some initializations
    init();
    
    // do the actual processing
    doTheJob();
    
    // perform final adjustments
    finalize(nlhs, plhs);
}
