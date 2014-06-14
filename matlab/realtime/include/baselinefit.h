/*=========================================================================
 * baselinefit.h
 * 
 *  Title: real-time baseline estimation for the segmentation procedure
 *  Author:     Diego Sogari
 *  Modified:   12/June/2014
 *
 *=======================================================================*/
#ifndef BASELINEFIT
#define BASELINEFIT

#include <math.h>
#include "mex.h"

#define BUFVAL(k)   (buffer[bufferLen-1+(k)])
#define MAX(a,b)    ((a) > (b) ? (a) : (b))
#define MIN(a,b)    ((a) < (b) ? (a) : (b))

/*=========================================================================
 * Type definitions
 *=======================================================================*/
template <class type>
class BaselineFit {
protected:
    mwSize bufferLen;
    int numSamples;
    double coeffs[2];
public:
    BaselineFit(mwSize buflen, int numsamp);
    ~BaselineFit();
    void newx(const type *buffer, int ponset, int toffset);
    double outputCoeff(int idx);
};

/*=========================================================================
 * Constructor
 *=======================================================================*/
template <class type>
BaselineFit<type>::BaselineFit(mwSize buflen, int numsamp)
{
    this->bufferLen = buflen;
    this->numSamples = numsamp;
    this->coeffs[0] = 0.0;
    this->coeffs[1] = 0.0;
}

/*=========================================================================
 * Destructor
 *=======================================================================*/
template <class type>
BaselineFit<type>::~BaselineFit()
{
}

/*=========================================================================
 * Update filter memory with an incoming sample
 *=======================================================================*/
template <class type>
void BaselineFit<type>::newx(const type *buffer, int ponset, int toffset)
{
    int i, begin, end, x[2];
    double y[2], mean;
    
    // calculate start and ending points
    begin = MAX(1 - bufferLen, ponset - (numSamples >> 1));
    end = MIN(0, toffset + (numSamples >> 1));
    
    // average of first numSamples
    mean = (type)0;
    for (i = begin; i < begin + numSamples; i++) {
        mean += (double)BUFVAL(i);
    }
    x[0] = 0.0;
    y[0] = mean / (double)numSamples;
    
    // average of last numSamples
    mean = 0.0;
    for (i = end - numSamples + 1; i <= end; i++) {
        mean += (double)BUFVAL(i);
    }
    x[1] = toffset - ponset; // begin - end;
    y[1] = mean / (double)numSamples;
    
    // create polynomial fit
    coeffs[0] = (y[1] - y[0]) / (double)(x[1] - x[0]);
    coeffs[1] = y[0] - (double)x[0] * coeffs[0];
}

/*=========================================================================
 * Return the output
 *=======================================================================*/
template <class type>
double BaselineFit<type>::outputCoeff(int idx)
{
    return this->coeffs[idx];
}

#endif