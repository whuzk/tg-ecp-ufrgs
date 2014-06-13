/*=========================================================================
 * peakdetector.h
 * 
 *  Title: real-time detection of peaks for the QRS detector
 *  Author:     Diego Sogari
 *  Modified:   12/June/2014
 *
 *=======================================================================*/
#ifndef PEAKDETECTOR
#define PEAKDETECTOR

#include "mex.h"

/*=========================================================================
 * Type definitions
 *=======================================================================*/
template <class type>
class PeakDetector {
protected:
    int peakIdx;            // index of the detected peak
    int lastPeakIdx;        // index of the last peak
    type lastPeakAmp;       // amplitude of the last peak
    bool wasPeakDetected;   // flag to indicate a new detected peak
    bool isSignalRising;    // flag to indicate a rise in the signal
public:
    PeakDetector();
    ~PeakDetector();
    void newx(type x);
    int outputPeakIndex();
    bool outputPeakDetected();
    bool outputSignalRising();
};

/*=========================================================================
 * Constructor
 *=======================================================================*/
template <class type>
PeakDetector<type>::PeakDetector()
{
    this->peakIdx = 0;
    this->lastPeakIdx = 0;
    this->lastPeakAmp = (type)0;
    this->wasPeakDetected = false;
    this->isSignalRising = false;
}

/*=========================================================================
 * Destructor
 *=======================================================================*/
template <class type>
PeakDetector<type>::~PeakDetector()
{
}

/*=========================================================================
 * Update filter memory with an incoming sample
 *=======================================================================*/
template <class type>
void PeakDetector<type>::newx(type x)
{
    wasPeakDetected = false;
    
    // check if the new amplitude is greater than that of the last peak
    if (x > lastPeakAmp) {
        // update state
        isSignalRising = true;
        lastPeakAmp = x;
        lastPeakIdx = 0;
    }
    else if (!isSignalRising) {
        // update current amplitude
        lastPeakAmp = x;
    }
    else if (x < lastPeakAmp / (type)2) {
        // report the new peak
        peakIdx = lastPeakIdx - 1;
        wasPeakDetected = true;
        // reset state
        isSignalRising = false;
        lastPeakAmp = x;
        lastPeakIdx = 0;
    }
    else lastPeakIdx--;
}

/*=========================================================================
 * Return the output
 *=======================================================================*/
template <class type>
int PeakDetector<type>::outputPeakIndex()
{
    return this->peakIdx;
}

/*=========================================================================
 * Return the output
 *=======================================================================*/
template <class type>
bool PeakDetector<type>::outputPeakDetected()
{
    return this->wasPeakDetected;
}

/*=========================================================================
 * Return the output
 *=======================================================================*/
template <class type>
bool PeakDetector<type>::outputSignalRising()
{
    return this->isSignalRising;
}

#endif