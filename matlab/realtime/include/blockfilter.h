/*=========================================================================
 * blockfilter.h
 * 
 *  Title: smoothing of a cardiac beat
 *  Author:     Diego Sogari
 *  Modified:   12/June/2014
 *
 *=======================================================================*/
#ifndef BLOCKFILTER
#define BLOCKFILTER

#include <string.h>
#include <math.h>
#include "mex.h"
#include "rtfilter.h"

/*=========================================================================
 * Type definitions
 *=======================================================================*/
template <class type>
class BlockFilter {
private:
    type *create_numerator(int len);
protected:
    mwSize bufferLen;
    RtFilter<type> *filter;
    int delay;
public:
    BlockFilter(mwSize bufferLen, const type *b, mwSize nb,
            const type *a, mwSize na, int delay);
    ~BlockFilter();
    void newx(const double *buffer, double *outbuffer);
};

/*=========================================================================
 * Constructor
 *=======================================================================*/
template <class type>
BlockFilter<type>::BlockFilter(mwSize buflen, const type *b, mwSize nb,
        const type *a, mwSize na, int delay)
{
    this->bufferLen = buflen;
    this->filter = new RtFilter<type>(b, nb, a, na);
    this->delay = delay;
}

/*=========================================================================
 * Destructor
 *=======================================================================*/
template <class type>
BlockFilter<type>::~BlockFilter()
{
    delete filter;
}

/*=========================================================================
 * Update filter memory with an incoming sample
 *=======================================================================*/
template <class type>
void BlockFilter<type>::newx(const double *buffer, double *outbuffer)
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