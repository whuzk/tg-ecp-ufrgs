/*=========================================================================
 * beatdetector.h
 * 
 *  Title: real-time detection of beats for the FP detector
 *  Author:     Diego Sogari
 *  Modified:   12/June/2014
 *
 *=======================================================================*/
#ifndef BEATDETECTOR
#define BEATDETECTOR

#include <math.h>
#include "mex.h"

/*=========================================================================
 * Type definitions
 *=======================================================================*/
class BeatDetector {
protected:
    mwSize firstRrInt;
    mwSize countDown;
    mwSize rrIntLeft;
    mwSize rrIntRight;
    int lastQrsIdx;
    bool beatDetected;
public:
    BeatDetector(double Fs);
    ~BeatDetector();
    void newx(int qrsIdx, bool qrsDetected);
    int outputLastQrsIdx();
    mwSize outputRrIntLeft();
    mwSize outputRrIntRight();
    bool outputBeatDetected();
};

/*=========================================================================
 * Constructor
 *=======================================================================*/
BeatDetector::BeatDetector(double Fs)
{
    this->firstRrInt = (int)Fs;
    this->countDown = -1;
    this->rrIntLeft = -1;
    this->rrIntRight = 0;
    this->lastQrsIdx = 0;
    this->beatDetected = false;
}

/*=========================================================================
 * Destructor
 *=======================================================================*/
BeatDetector::~BeatDetector()
{
}

/*=========================================================================
 * Update filter memory with an incoming sample
 *=======================================================================*/
void BeatDetector::newx(int qrsIdx, bool qrsDetected)
{
    // simulate QRS detection
    if ((qrsDetected && countDown > 0) || countDown == 0) {
        if (countDown > 0) {
            rrIntRight = qrsIdx - lastQrsIdx;
        }
        else rrIntRight = rrIntLeft;
        beatDetected = true;
    }
    else beatDetected = false;
    
    if (qrsDetected) {
        // update rr
        if (rrIntLeft == -1) {
            rrIntLeft = firstRrInt;
        }
        else rrIntLeft = qrsIdx - lastQrsIdx;
        // update qrs
        lastQrsIdx = qrsIdx;
        // update countdown
        countDown = rrIntLeft;
    }
    else if (countDown >= 0) {
        // update countdown
        countDown--;
    }
    
    lastQrsIdx--;
}

/*=========================================================================
 * Return the output
 *=======================================================================*/
int BeatDetector::outputLastQrsIdx()
{
    return this->lastQrsIdx;
}

/*=========================================================================
 * Return the output
 *=======================================================================*/
mwSize BeatDetector::outputRrIntLeft()
{
    return this->rrIntLeft;
}

/*=========================================================================
 * Return the output
 *=======================================================================*/
mwSize BeatDetector::outputRrIntRight()
{
    return this->rrIntRight;
}

/*=========================================================================
 * Return the output
 *=======================================================================*/
bool BeatDetector::outputBeatDetected()
{
    return this->beatDetected;
}

#endif