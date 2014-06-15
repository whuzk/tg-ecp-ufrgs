/*=========================================================================
 * diffbeat.h
 * 
 *  Title: differentiation of a cardiac beat
 *  Author:     Diego Sogari
 *  Modified:   12/June/2014
 *
 *=======================================================================*/
#ifndef DIFFBEAT
#define DIFFBEAT

#include <math.h>
#include "mex.h"
#include "rtfilter.h"

/*=========================================================================
 * Type definitions
 *=======================================================================*/
template <class type>
class DiffBeat {
protected:
    mwSize bufferLen;
    RtFilter<type> *filter;
public:
    DiffBeat(mwSize bufferLen);
    ~DiffBeat();
    void newx(const double *buffer, double *outbuffer);
};

/*=========================================================================
 * Constructor
 *=======================================================================*/
template <class type>
DiffBeat<type>::DiffBeat(mwSize buflen)
{
    type a = (type)1;
    type b[2] = {(type)1, (type)-1};
    this->bufferLen = buflen;
    this->filter = new RtFilter<type>(b, 2, &a, 1);
}

/*=========================================================================
 * Destructor
 *=======================================================================*/
template <class type>
DiffBeat<type>::~DiffBeat()
{
    delete filter;
}

/*=========================================================================
 * Update filter memory with an incoming sample
 *=======================================================================*/
template <class type>
void DiffBeat<type>::newx(const double *buffer, double *outbuffer)
{
    mwSize i;
    
    for (i = 0; i < bufferLen; i++) {
        outbuffer[i] = filter->newx(buffer[i]);
    }
}

#endif