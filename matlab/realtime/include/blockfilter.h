/*=========================================================================
 * blockfilter.h
 * 
 *  Title: implementation of filtering on a block without delay
 *  Author:     Diego Sogari
 *  Modified:   12/June/2014
 *
 *=======================================================================*/
#ifndef BLOCKFILTER
#define BLOCKFILTER

#include <math.h>
#include "mex.h"
#include "rtfilter.h"

/*=========================================================================
 * Type definitions
 *=======================================================================*/
template <class type>
class BlockFilter {
protected:
    mwSize bufferLen;
    RtFilter<type> *filter;
    int delay;
public:
    BlockFilter(mwSize bufferLen, const type *b, mwSize nb,
            const type *a, mwSize na, int delay);
    ~BlockFilter();
    void newx(const type *buffer, type *outbuffer);
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
void BlockFilter<type>::newx(const type *buffer, type *outbuffer)
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