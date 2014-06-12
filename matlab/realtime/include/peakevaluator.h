/*=========================================================================
 * peakevaluator.h
 * 
 *  Title: real-time evaluation of peaks for the QRS detector
 *  Author:     Diego Sogari
 *  Modified:   12/June/2014
 *
 *=======================================================================*/
#ifndef PEAKEVALUATOR
#define PEAKEVALUATOR

#include <math.h>
#include "mex.h"

/*=========================================================================
 * Type definitions
 *=======================================================================*/
template <class type>
class PeakEvaluator {
protected:
    type signalLevel;
    type noiseLevel;
    type estimRatio;
    mwSize bufferLen;
    mwSize lastPeakIdx;
    bool qrsDetected;
    type estimate(type oldVal, type newVal);
public:
    PeakEvaluator(type ratio, mwSize len);
    ~PeakEvaluator();
    void newx(const type *buffer, mwSize peakIdx, type thresh,
            bool searchback, bool training);
    type outputSignalLevel();
    type outputNoiseLevel();
    bool outputQrsDetected();
};

/*=========================================================================
 * Constructor
 *=======================================================================*/
template <class type>
PeakEvaluator<type>::PeakEvaluator(type ratio, mwSize len)
{
    this->signalLevel = (type)0;
    this->noiseLevel = (type)0;
    this->estimRatio = ratio;
    this->bufferLen = len;
    this->lastPeakIdx = 0;
    this->qrsDetected = false;
}

/*=========================================================================
 * Destructor
 *=======================================================================*/
template <class type>
PeakEvaluator<type>::~PeakEvaluator()
{
}

/*=========================================================================
 * Iterative estimator (default for real numbers)
 *=======================================================================*/
template <class type>
type PeakEvaluator<type>::estimate(type oldVal, type newVal)
{
    return (1 - estimRatio) * oldVal + estimRatio * newVal;
}

/*=========================================================================
 * Iterative estimator (for integers)
 *=======================================================================*/
int PeakEvaluator<int>::estimate(int oldVal, int newVal)
{
    return (((1 << estimRatio) - 1) * oldVal + newVal) >> estimRatio;
}

/*=========================================================================
 * Update filter memory with an incoming sample
 *=======================================================================*/
template <class type>
void PeakEvaluator<type>::newx(const type *buffer, mwSize peakIdx,
        type thresh, bool searchback, bool training)
{
    type peakAmp = buffer[bufferLen - 1 + peakIdx];
    
    if (searchback) {
        type save = estimRatio;
        if (estimRatio > 1) {
            estimRatio -= 1;
        }
        else estimRatio *= 2;
        signalLevel = estimate(signalLevel, peakAmp);
        estimRatio = save;
        qrsDetected = true;
    }
    else if (lastPeakIdx != peakIdx) {
        if (peakAmp >= thresh) {
            if (training) {
                signalLevel = (type)fmax((double)signalLevel, (double)peakAmp);
            }
            else signalLevel = estimate(signalLevel, peakAmp);
            qrsDetected = true;
        }
        else {
            if (training) {
                noiseLevel = (type)fmax((double)noiseLevel, (double)peakAmp);
            }
            else noiseLevel = estimate(noiseLevel, peakAmp);
            qrsDetected = false;
        }
        lastPeakIdx = peakIdx;
    }
    else qrsDetected = false;
}

/*=========================================================================
 * Return the output
 *=======================================================================*/
template <class type>
type PeakEvaluator<type>::outputSignalLevel()
{
    return this->signalLevel;
}

/*=========================================================================
 * Return the output
 *=======================================================================*/
template <class type>
type PeakEvaluator<type>::outputNoiseLevel()
{
    return this->noiseLevel;
}

/*=========================================================================
 * Return the output
 *=======================================================================*/
template <class type>
bool PeakEvaluator<type>::outputQrsDetected()
{
    return this->qrsDetected;
}

#endif