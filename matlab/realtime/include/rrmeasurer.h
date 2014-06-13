/*=========================================================================
 * rrmeasurer.h
 * 
 *  Title: real-time measuring of RR intervals for the QRS detector
 *  Author:     Diego Sogari
 *  Modified:   12/June/2014
 *
 *=======================================================================*/
#ifndef RRMEASURER
#define RRMEASURER

#include <string.h>
#include <math.h>
#include "mex.h"

#define RR_BUFLEN       8
#define LOG2_BUFLEN     3
#define RR1BUFVAL       (rr1buffer[rrIdx1&(RR_BUFLEN-1)])
#define RR2BUFVAL       (rr2buffer[rrIdx2&(RR_BUFLEN-1)])

/*=========================================================================
 * Type definitions
 *=======================================================================*/
class RrMeasurer {
private:
    mwSize rr1buffer[RR_BUFLEN];
    mwSize rr2buffer[RR_BUFLEN];
protected:
    bool initialized;
    mwSize rrMean1;
    mwSize rrMean2;
    mwSize rrLow;
    mwSize rrHigh;
    mwSize rrMiss;
    mwSize rrAcc1;
    mwSize rrAcc2;
    int rrIdx1;
    int rrIdx2;
    int lastPeakIdx;
public:
    RrMeasurer();
    ~RrMeasurer();
    void newx(mwSize newRR, bool qrsDetected, bool qrsConfirmed);
    mwSize outputRrMean();
    mwSize outputRrMiss();
};

/*=========================================================================
 * Constructor
 *=======================================================================*/
RrMeasurer::RrMeasurer()
{
    this->initialized = false;
    this->rrMean1 = 0;
    this->rrMean2 = 0;
    this->rrLow = 0;
    this->rrHigh = 0;
    this->rrMiss = 0;
    this->rrAcc1 = 0;
    this->rrAcc2 = 0;
    this->rrIdx1 = 0;
    this->rrIdx2 = 0;
    this->lastPeakIdx = 0;
}

/*=========================================================================
 * Destructor
 *=======================================================================*/
RrMeasurer::~RrMeasurer()
{
}

/*=========================================================================
 * Update filter memory with an incoming sample
 *=======================================================================*/
void RrMeasurer::newx(int peakIdx, bool qrsDetected, bool qrsConfirmed)
{
    if (qrsConfirmed) {
        // calculate new RR interval
        mwSize newRR = peakIdx - lastPeakIdx;
        
        // update RR buffers
        if (!initialized) {
            int i;
            for (i = 0; i < RR_BUFLEN; i++) {
                rr1buffer[i] = newRR;
                rr2buffer[i] = newRR;
            }
            rrAcc1 = rrAcc2 = newRR << LOG2_BUFLEN;
            initialized = true;
        }
        else {
            rrAcc1 += newRR - RR1BUFVAL;
            RR1BUFVAL = newRR;
            rrIdx1++;
            if (rrLow <= newRR && newRR <= rrHigh) {
                rrAcc2 += newRR - RR2BUFVAL;
                RR2BUFVAL = newRR;
                rrIdx2++;
            }
        }
        
        // update RR averages
        rrMean1 = rrAcc1 >> LOG2_BUFLEN;
        rrMean2 = rrAcc2 >> LOG2_BUFLEN;
        rrLow = rrMean2 - ((rrAcc2 - rrMean2) >> 6);    // 0.89*RRmean
        rrHigh = rrLow + ((rrAcc2 + rrMean2) >> 5);     // 1.17*RRmean
        rrMiss = rrHigh + (rrMean2 >> 2);               // 1.67*RRmean
    }
    
    if (qrsDetected) {
        lastPeakIdx = peakIdx;
    }
    else lastPeakIdx--;
}

/*=========================================================================
 * Return the output
 *=======================================================================*/
mwSize RrMeasurer::outputRrMean()
{
    return this->rrMean2;
}

/*=========================================================================
 * Return the output
 *=======================================================================*/
mwSize RrMeasurer::outputRrMiss()
{
    return this->rrMiss;
}

#endif