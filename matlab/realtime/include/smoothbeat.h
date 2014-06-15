/*=========================================================================
 * smoothbeat.h
 * 
 *  Title: smoothing of a cardiac beat
 *  Author:     Diego Sogari
 *  Modified:   12/June/2014
 *
 *=======================================================================*/
#ifndef SMOOTHBEAT
#define SMOOTHBEAT

#include <string.h>
#include <math.h>
#include "mex.h"
#include "rtfilter.h"

/*=========================================================================
 * Type definitions
 *=======================================================================*/
template <class type>
class SmoothBeat {
private:
    type *create_numerator(int len);
protected:
    mwSize bufferLen;
    mwSize delay;
    RtFilter<type> *filter;
public:
    SmoothBeat(mwSize bufferLen, int width);
    ~SmoothBeat();
    void newx(const double *buffer, double *outbuffer);
};

/*=========================================================================
 * Constructor
 *=======================================================================*/
template <class type>
type *SmoothBeat<type>::create_numerator(int len)
{
    type *b = new type[len];
    memset(b, 0, len * sizeof(type));
    b[0] = (type)1;
    b[len - 1] = (type)-1;
    return b;
}

/*=========================================================================
 * Constructor
 *=======================================================================*/
template <class type>
SmoothBeat<type>::SmoothBeat(mwSize buflen, int width)
{
    type a[2] = {(type)1, (type)-1};
    type *b = create_numerator(width + 1);
    this->bufferLen = buflen;
    this->delay = width >> 1;
    this->filter = new RtFilter<type>(b, width + 1, a, 2);
    delete[] b;
}

/*=========================================================================
 * Destructor
 *=======================================================================*/
template <class type>
SmoothBeat<type>::~SmoothBeat()
{
    delete filter;
}

/*=========================================================================
 * Update filter memory with an incoming sample
 *=======================================================================*/
template <class type>
void SmoothBeat<type>::newx(const double *buffer, double *outbuffer)
{
    mwSize i;
    
    for (i = 0; i < delay; i++) {
        filter->newx(buffer[i]);
    }
    for (i = delay; i < bufferLen; i++) {
        outbuffer[i - delay] = filter->newx(buffer[i]);
    }
    for (i = bufferLen; i < bufferLen + delay; i++) {
        outbuffer[i - delay] = filter->newx(0);
    }
}

#endif