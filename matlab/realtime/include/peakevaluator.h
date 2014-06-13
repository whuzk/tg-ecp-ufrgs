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

#define IDXOK(k)    (bufferLen+(k) > 0)
#define BUFVAL(k)   (buffer[bufferLen-1+(k)])
#define MAX(a,b)    ((a) > (b) ? (a) : (b))

/*=========================================================================
 * Type definitions
 *=======================================================================*/
template <class type>
class PeakEvaluator {
protected:
    mwSize bufferLen;
    type signalLevel;
    type noiseLevel;
    type estimFactor;
    bool qrsDetected;
    type estimate(type oldVal, type newVal);
public:
    PeakEvaluator(mwSize buflen, type factor);
    ~PeakEvaluator();
    void newx(const type *buffer, int peakIdx, type thresh,
            bool normalpeak, bool searchback, bool training);
    type outputSignalLevel();
    type outputNoiseLevel();
    bool outputQrsDetected();
};

/*=========================================================================
 * Constructor
 *=======================================================================*/
template <class type>
PeakEvaluator<type>::PeakEvaluator(mwSize buflen, type factor)
{
    this->bufferLen = buflen;
    this->signalLevel = (type)0;
    this->noiseLevel = (type)0;
    this->estimFactor = factor;
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
 * Iterative estimator
 *=======================================================================*/
template <class type>
type PeakEvaluator<type>::estimate(type oldVal, type newVal)
{
    return ((estimFactor - 1) * oldVal + newVal) / estimFactor;
}

/*=========================================================================
 * Update filter memory with an incoming sample
 *=======================================================================*/
template <class type>
void PeakEvaluator<type>::newx(const type *buffer, int peakIdx,
        type thresh, bool normalpeak, bool searchback, bool training)
{
    qrsDetected = false;
    
    if (IDXOK(peakIdx)) {
        type peakAmp = BUFVAL(peakIdx);
        
        if (searchback) {
            type save = estimFactor;
            estimFactor /= (type)2;
            if (peakAmp >= thresh / (type)2) {
                signalLevel = estimate(signalLevel, peakAmp);
                qrsDetected = true;
            }
            else noiseLevel = estimate(noiseLevel, peakAmp);
            estimFactor = save;
        }
        else if (normalpeak) {
            if (peakAmp >= thresh) {
                if (training) {
                    signalLevel = MAX(signalLevel, peakAmp);
                }
                else signalLevel = estimate(signalLevel, peakAmp);
                qrsDetected = true;
            }
            else if (training) {
                noiseLevel = MAX(noiseLevel, peakAmp);
            }
            else noiseLevel = estimate(noiseLevel, peakAmp);
        }
    } 
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