/*=========================================================================
 * qrsevaluator.h
 * 
 *  Title: real-time evaluation of QRS for the QRS detector
 *  Author:     Diego Sogari
 *  Modified:   12/June/2014
 *
 *=======================================================================*/
#ifndef QRSEVALUATOR
#define QRSEVALUATOR

#include <math.h>
#include "mex.h"

#define IDXOK(k)    (bufferLen+(k) > 0)
#define BUFVAL(k)   (buffer[bufferLen-1+(k)])

/*=========================================================================
 * Type definitions
 *=======================================================================*/
template <class type>
class QrsEvaluator {
private:
    type maxdiff(const type *buffer, int start, mwSize len);
protected:
    mwSize bufferLen;
    int lastQrsIdx;
    mwSize refracPeriod;
    mwSize qrsHalfLength;
    mwSize twaveTolerance;
    bool qrsConfirmed;
public:
    QrsEvaluator(mwSize buflen, double Fs);
    ~QrsEvaluator();
    void newx(const type *buffer, int qrsIdx, bool qrsDetected);
    bool outputQrsConfirmed();
};

/*=========================================================================
 * Constructor
 *=======================================================================*/
template <class type>
QrsEvaluator<type>::QrsEvaluator(mwSize buflen, double Fs)
{
    this->bufferLen = buflen;
    this->lastQrsIdx = 0;
    this->refracPeriod = (mwSize)(0.20 * Fs);
    this->qrsHalfLength = (mwSize)(0.10 * Fs);
    this->twaveTolerance = (mwSize)(0.36 * Fs);
    this->qrsConfirmed = false;
}

/*=========================================================================
 * Destructor
 *=======================================================================*/
template <class type>
QrsEvaluator<type>::~QrsEvaluator()
{
}

/*=========================================================================
 * Get the maximum slope on the detection signal 
 *=======================================================================*/
template <class type>
type QrsEvaluator<type>::maxdiff(const type *buffer, int start, mwSize len)
{
    type d = 0;
    type newd;
    mwSize i;
    
    for (i = 1; i < len; i++) {
        newd = BUFVAL(start + i) - BUFVAL(start + i - 1);
        if (newd > d) {
            d = newd;
        }
    }
    return d;
}

/*=========================================================================
 * Update filter memory with an incoming sample
 *=======================================================================*/
template <class type>
void QrsEvaluator<type>::newx(const type *buffer, int qrsIdx,
        bool qrsDetected)
{
    qrsConfirmed = false;
    
    if (qrsDetected) {
        mwSize rr = qrsIdx - lastQrsIdx;
        
        if (rr > twaveTolerance) {
            qrsConfirmed = true;
        }
        else if (rr > refracPeriod) {
            // calculate starting indices
            int start1 = qrsIdx - qrsHalfLength + 1;
            int start2 = lastQrsIdx - qrsHalfLength + 1;

            // check index validity
            if (IDXOK(start1) && IDXOK(start2)) {
                // check condition for T wave
                type slope1 = maxdiff(buffer, start1, qrsHalfLength);
                type slope2 = maxdiff(buffer, start2, qrsHalfLength);
                qrsConfirmed = (slope1 >= slope2 / (type)2);
            }
        }
        lastQrsIdx = qrsIdx;
    }
    
    lastQrsIdx--;
}

/*=========================================================================
 * Return the output
 *=======================================================================*/
template <class type>
bool QrsEvaluator<type>::outputQrsConfirmed()
{
    return this->qrsConfirmed;
}

#endif