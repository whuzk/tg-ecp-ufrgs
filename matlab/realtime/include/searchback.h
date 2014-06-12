/*=========================================================================
 * searchback.h
 * 
 *  Title: real-time searchback procedure for the QRS detector
 *  Author:     Diego Sogari
 *  Modified:   12/June/2014
 *
 *=======================================================================*/
#ifndef SEARCHBACK
#define SEARCHBACK

#include <math.h>
#include "mex.h"

/*=========================================================================
 * Type definitions
 *=======================================================================*/
template <class type>
class SearchBack {
protected:
    type signalLevel;
    type noiseLevel;
    type estimRatio;
    mwSize lastPeakIdx;
    bool qrsDetected;
    type estimate(type oldVal, type newVal);
public:
    SearchBack(type ratio);
    ~SearchBack();
    void newx();
    //type outputSignalLevel();
    //type outputNoiseLevel();
    //bool outputQrsDetected();
};

/*=========================================================================
 * Constructor
 *=======================================================================*/
template <class type>
SearchBack<type>::SearchBack(type ratio)
{
}

/*=========================================================================
 * Destructor
 *=======================================================================*/
template <class type>
SearchBack<type>::~SearchBack()
{
}

/*=========================================================================
 * Iterative estimator (default for real numbers)
 *=======================================================================*/
template <class type>
type SearchBack<type>::estimate(type oldVal, type newVal)
{
    return (1 - estimRatio) * oldVal + estimRatio * newVal;
}

/*=========================================================================
 * Iterative estimator (default for real numbers)
 *=======================================================================*/
int SearchBack<int>::estimate(int oldVal, int newVal)
{
    return (((1 << estimRatio) - 1) * oldVal + newVal) >> estimRatio;
}

/*=========================================================================
 * Update filter memory with an incoming sample
 *=======================================================================*/
template <class type>
void SearchBack<type>::newx(type peakAmp, mwSize peakIdx, type thresh,
        bool training, bool searchback)
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
 * Return the output
 *=======================================================================*/
template <class type>
type SearchBack<type>::outputSignalLevel()
{
    return this->signalLevel;
}

/*=========================================================================
 * Return the output
 *=======================================================================*/
template <class type>
type SearchBack<type>::outputNoiseLevel()
{
    return this->noiseLevel;
}

/*=========================================================================
 * Return the output
 *=======================================================================*/
template <class type>
bool SearchBack<type>::outputQrsDetected()
{
    return this->qrsDetected;
}

#endif