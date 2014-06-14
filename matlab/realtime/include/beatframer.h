/*=========================================================================
 * beatframer.h
 * 
 *  Title: real-time framing for the segmentation procedure
 *  Author:     Diego Sogari
 *  Modified:   12/June/2014
 *
 *=======================================================================*/
#ifndef BEATFRAMER
#define BEATFRAMER

#include <string.h>
#include <math.h>
#include "mex.h"

#define BUFVAL(k)   (buffer[bufferLen-1+(k)])
#define MAX(a,b)    ((a) > (b) ? (a) : (b))
#define MIN(a,b)    ((a) < (b) ? (a) : (b))

/*=========================================================================
 * Type definitions
 *=======================================================================*/
template <class type>
class BeatFramer {
protected:
    mwSize bufferLen;
    mwSize frameLen;
public:
    BeatFramer(mwSize buflen, mwSize framelen);
    ~BeatFramer();
    void newx(const type *buffer, type *outbuffer, int rpeak, int ponset,
            int toffset, const double *p);
};

/*=========================================================================
 * Constructor
 *=======================================================================*/
template <class type>
BeatFramer<type>::BeatFramer(mwSize buflen, mwSize framelen)
{
    this->bufferLen = buflen;
    this->frameLen = framelen;
}

/*=========================================================================
 * Destructor
 *=======================================================================*/
template <class type>
BeatFramer<type>::~BeatFramer()
{
}

/*=========================================================================
 * Update filter memory with an incoming sample
 *=======================================================================*/
template <class type>
void BeatFramer<type>::newx(const type *buffer, type *frame, int rpeak,
        int ponset, int toffset, const double *p)
{
    mwSize len1, len2, halflen, beatlen, pad, cut1, cut2;
    int i, j, r, start, end;
    type baseline;
    
    // calculate padding and cutting sizes
    halflen = frameLen >> 1;
    beatlen = toffset - ponset + 1;
    r = rpeak - ponset + 1;
    pad = MAX(0, halflen - r + 1);
    cut1 = MAX(0, r - 1 - halflen);
    cut2 = MAX(0, beatlen - r - halflen);
    
    // update start and end indices
    start = ponset + cut1;
    end = toffset - cut2;
    
    // calculate new amplitudes
    memset(frame, 0, frameLen * sizeof(type));
    for (i = start, j = pad; i <= end; i++, j++) {
        baseline = (p[0] * (i - ponset) + p[1]);
        frame[j] = BUFVAL(i) - baseline;
    }
}

#endif