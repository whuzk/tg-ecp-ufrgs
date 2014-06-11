/*=========================================================================
 * qrsdetector.h
 * 
 *  Title: real-time detection of QRS complex
 *  Author:     Diego Sogari
 *  Modified:   11/June/2014
 *
 *=======================================================================*/
#ifndef QRSDETECTOR
#define QRSDETECTOR

#include <stdlib.h>
#include <math.h>

/*=========================================================================
 * Constants
 *=======================================================================*/
#define RR_LEN  8   // eight most recent RR intervals

/*=========================================================================
 * Type definitions
 *=======================================================================*/
typedef struct {
    unsigned int ci;            // current sample index
    int *buffer;                // detection buffer
    int sigThreshold;           // signal threshold
    int signalLevel;            // signal level
    int noiseLevel;             // noise level
    int estRatio;               // ratio of signal/noise level estimation
    int peakAmp;                // currently detected peak amplitude
    int lastPeakAmp;            // amplitude of the last peak
    mwSize lastPeakIdx;         // index of the last peak
    mwSize rrIntMean1;          // primary running average of RR
    mwSize rrIntMean2;          // secondary running average of RR
    mwSize rrIntMiss;           // interval for QRS to be assumed as missed
    mwSize trainCountDown;      // countdown of training period
    mwSize qrsHalfLength;       // half the length of a QRS (predefined)
    mwSize twaveTolerance;      // tolerance for T-wave detection (predefined)
    mwSize refractoryPeriod;    // length of refractory period (predefined)
    mwSize lastQrsIdx;          // location of last detected QRS
    mwSize searchBackIdx;       // location of searchback starting point
    mwSize peakIdx;             // currently detected peak location
    mwSize bufLen;              // length of detection buffer
    mwSize rr1Buf[RR_LEN];      // primary RR buffer
    mwSize rr2Buf[RR_LEN];      // secondary RR buffer
    bool isSignalRising;        // flag to indicate a rise in the signal
    bool rrInitialized;         // flag to tell if RR measures are init'd
    mwSize rrIntLow;            // lower limit for the RR interval
    mwSize rrIntHigh;           // upper limit for the RR interval
    mwSize rrIdx1;              // index in the primary RR buffer
    mwSize rrIdx2;              // index in the secondary RR buffer
    mwSize rrY1;                // internal variable for RR estimation
    mwSize rrY2;                // internal variable for RR estimation
} qrsdetobject;

/*=========================================================================
 * Macros
 *=======================================================================*/
#define RR1VAL(o,k)     ((o).rr1Buf[(k)&(RR_LEN-1)])
#define RR2VAL(o,k)     ((o).rr2Buf[(k)&(RR_LEN-1)])
#define DETVAL(o,k)     ((o).buffer[((o).ci+(k))&((o).bufLen-1)])
#define IDXOK(o,k)      ((k)<=0 && 0<(o).bufLen+(k))
#define ESTIM(x,r,v)    ((((1<<(r))-1)*(x)+(v))>>(r))
#define ESTSIG(o,r)     (ESTIM((o).signalLevel,r,(o).peakAmp))
#define ESTNOIS(o,r)    (ESTIM((o).noiseLevel,r,(o).peakAmp))
#define ESTTHR(o,r)     ((o).noiseLevel+(abs((o).signalLevel-(o).noiseLevel)>>r))

#define initqrsdetector(o)  \
    (o).buffer = NULL;      \
    (o).lastPeakAmp = 0;    \
    (o).lastPeakIdx = 0;    \
    (o).rrInitialized = false;
    
#define endqrsdetector(o)   \
    free((o).buffer);

/*=========================================================================
 * Get the maximum slope on the detection signal 
 *=======================================================================*/
static int maxdiff(qrsdetobject *object, mwSize start, mwSize len)
{
    int d = 0;
    int newd;
    mwSize i;
    
    for (i = 1; i < len; i++) {
        newd = DETVAL(*object, start + i) - DETVAL(*object, start + i - 1);
        if (newd > d) {
            d = newd;
        }
    }
    return d;
}

/*=========================================================================
 * Find the position of the maximum value on the detection signal 
 *=======================================================================*/
static mwSize findmax(qrsdetobject *object, mwSize start, mwSize len)
{
    mwSize i, pos = 0;
    int y, newy;
    
    y = DETVAL(*object, start);
    for (i = 1; i < len; i++) {
        newy = DETVAL(*object, start + i);
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
static bool istwave(qrsdetobject *object)
{
    mwSize candQrs = object->peakIdx;
    mwSize lastQrs = object->lastQrsIdx;
    mwSize rr = candQrs - lastQrs;
    
    if (rr <= object->refractoryPeriod) {
        return true;
    }
    else if (rr > object->twaveTolerance) {
        return false;
    }
    else {
        // calculate starting indices
        mwSize start1 = candQrs - object->qrsHalfLength + 1;
        mwSize start2 = lastQrs - object->qrsHalfLength + 1;
        
        // check index validity
        if (IDXOK(*object, start1) && IDXOK(*object, start2)) {
            // check condition for T wave
            int slope1 = maxdiff(object, start1, object->qrsHalfLength);
            int slope2 = maxdiff(object, start2, object->qrsHalfLength);
            return (slope1 < (slope2 >> 1));
        }
        else return false;
    }
}

/*=========================================================================
 * Peak detection 
 *=======================================================================*/
static bool detectPeak(qrsdetobject *object)
{
    int newAmp = DETVAL(*object, 0);
    
    // check if the new amplitude is greater than that of the last peak
    if (newAmp > object->lastPeakAmp) {
        // signalize beginning or continuation of positive slope
        object->isSignalRising = true;
        // update peak info
        object->lastPeakAmp = newAmp;
        object->lastPeakIdx = 0;
        return false;
    }
    else if (!object->isSignalRising) {
        // update current amplitude
        object->lastPeakAmp = newAmp;
        return false;
    }
    else if (newAmp < (object->lastPeakAmp >> 1)) {
        // report the new peak amplitude and location
        object->peakAmp = object->lastPeakAmp;
        object->peakIdx = object->lastPeakIdx - 1;
        // reset state of positive slope
        object->isSignalRising = false;
        // reset peak info
        object->lastPeakAmp = newAmp;
        object->lastPeakIdx = 0;
        // signalize detection
        return true;
    }
    else {
        object->lastPeakIdx--;
        return false;
    }
}

/*=========================================================================
 * Peak evaluation 
 *=======================================================================*/
static bool evaluatePeak(qrsdetobject *object)
{
    if (object->trainCountDown > 0) {
        if (object->peakAmp >= object->sigThreshold) {
            object->signalLevel = max(object->signalLevel,
                    object->peakAmp);
            object->lastQrsIdx = object->peakIdx;
        }
        else {
            object->noiseLevel = max(object->noiseLevel, object->peakAmp);
        }
        return false;
    }
    else if (object->peakAmp >= object->sigThreshold) {
        object->signalLevel = ESTSIG(*object, object->estRatio);
        return !istwave(object);
    }
    else {
        object->noiseLevel = ESTNOIS(*object, object->estRatio);
        return false;
    }
}

/*=========================================================================
 * Search back procedure 
 *=======================================================================*/
static bool searchBack(qrsdetobject *object)
{
    if ((!object->isSignalRising ||
                DETVAL(*object, 0) < (object->sigThreshold >> 1)) &&
            object->searchBackIdx != 0 &&
            object->searchBackIdx + object->rrIntMiss <= 0) {
        // calculate starting point for the search
        mwSize begin = 0 - object->rrIntMean2 + 1;
        
        if (IDXOK(*object, begin)) {
            // search back and locate the max in this interval
            object->peakIdx = findmax(object, begin, object->rrIntMean2);
            object->peakAmp = DETVAL(*object, object->peakIdx);
            
            // check if candidate peak is from qrs
            if (object->peakAmp < object->sigThreshold &&
                    object->peakAmp >= (object->sigThreshold >> 1) &&
                    !istwave(object)) {
                // adjust signal level
                object->signalLevel = ESTSIG(*object, 2);
                // signalize qrs detection
                return true;
            }
            else {
                // reduce levels by half
                object->signalLevel >>= 1;
                object->noiseLevel >>= 1;
                // postpone searchback
                object->searchBackIdx += object->rrIntMean2;
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
static void updateRRInfo(qrsdetobject *object)
{
    mwSize mean1, mean2;
    mwSize newRR = object->peakIdx - object->lastQrsIdx;
    
    if (!object->rrInitialized) {
        mwSize i;
        for (i = 0; i < RR_LEN; i++) {
            object->rr1Buf[i] = newRR;
            object->rr2Buf[i] = newRR;
        }
        object->rrY1 = object->rrY2 = newRR << 3;
        object->rrIdx1 = object->rrIdx2 = 0;
        object->rrInitialized = true;
    }
    else {
        object->rrY1 += newRR - RR1VAL(*object, object->rrIdx1);
        RR1VAL(*object, object->rrIdx1++) = newRR;
        if (object->rrIntLow <= newRR && newRR <= object->rrIntHigh) {
            object->rrY2 += newRR - RR2VAL(*object, object->rrIdx2);
            RR2VAL(*object, object->rrIdx2++) = newRR;
        }
    }
    
    // update RR averages
    mean1 = object->rrY1 >> 3;
    mean2 = object->rrY2 >> 3;
    object->rrIntLow = object->rrIntMean2 - ((object->rrY2 - mean2) >> 6);
    object->rrIntHigh = object->rrIntLow + ((object->rrY2 + mean2) >> 5);
    object->rrIntMiss = object->rrIntHigh + (mean2 >> 2);
    object->rrIntMean1 = mean1;
    object->rrIntMean2 = mean2;
}

/*=========================================================================
 * Update QRS detector with an incoming sample, and return the output.
 *=======================================================================*/
static bool qrsdetnewx(qrsdetobject *object, int sample, mwSize *rpeak)
{
    bool peakDetected, qrsDetected;
    
    // store new sample
    DETVAL(*object, 0) = sample;
    
    // detect peak
    peakDetected = detectPeak(object);
    
    // detect qrs
    qrsDetected = peakDetected && evaluatePeak(object);
    
    // search back
    qrsDetected = qrsDetected || searchBack(object);
    
    // update threshold
    object->sigThreshold = ESTTHR(*object, 2);
    
    // update qrs info
    if (qrsDetected) {
        updateRRInfo(object);
        object->lastQrsIdx = object->peakIdx;
        object->searchBackIdx = object->peakIdx;
    }
    
    // update training countdown
    if (object->trainCountDown > 0) {
        object->trainCountDown--;
        // decrease estimation ratio for times beyond the training period
        if (object->trainCountDown == 0) {
            object->estRatio++;
        }
    }
    
    // update indices
    object->lastQrsIdx--;
    object->searchBackIdx--;
    
    // increment sample index
    object->ci++;
    
    // save result
    *rpeak = object->lastQrsIdx;
    
    return qrsDetected;
}

/*=========================================================================
 * Initialize a QRS detector object and allocate memory for its properties
 *=======================================================================*/
static void create_qrsdet(qrsdetobject *object, double Fs)
{
    //TODO
}

#endif