/*=========================================================================
 * edgedetector.h
 * 
 *  Title: detection of edge
 *  Author:     Diego Sogari
 *  Modified:   12/June/2014
 *
 *=======================================================================*/
#ifndef EDGEDETECTOR
#define EDGEDETECTOR

#include <math.h>
#include "mex.h"

/*=========================================================================
 * Type definitions
 *=======================================================================*/
template <class type>
class EdgeDetector {
protected:
    mwSize bufferLen;
    int windowLen;
    type threshold;
    int edgeIdx;
public:
    EdgeDetector(mwSize bufferLen, int width, type thresh);
    ~EdgeDetector();
    void newx(const type *buffer, int start, int end, int def);
    int outputEdgeIdx();
};

/*=========================================================================
 * Constructor
 *=======================================================================*/
template <class type>
EdgeDetector<type>::EdgeDetector(mwSize buflen, int width, type thresh)
{
    this->bufferLen = buflen;
    this->windowLen = width;
    this->threshold = thresh;
    this->edgeIdx = 0;
}

/*=========================================================================
 * Destructor
 *=======================================================================*/
template <class type>
EdgeDetector<type>::~EdgeDetector()
{
}

/*=========================================================================
 * Update filter memory with an incoming sample
 *=======================================================================*/
template <class type>
void EdgeDetector<type>::newx(const type *buffer, int start, int end,
        int def)
{
    int half = windowLen >> 1;
    int firstAbove = 0;
    int inc = 1;
    int i;
    
    if (start > end) {
        inc = -inc;
        half = -half;
    }
    
    edgeIdx = def;
    for (i = start; i != end + inc; i += inc) {
        if (abs(buffer[i]) > threshold) {
            firstAbove = 1;
        }
        else if (firstAbove == windowLen) {
            edgeIdx = i - half;
            break;
        }
        else firstAbove++;
    }
}

/*=========================================================================
 * Return the output
 *=======================================================================*/
template <class type>
int EdgeDetector<type>::outputEdgeIdx()
{
    return this->edgeIdx;
}

#endif