/*=========================================================================
 * searchback.h
 * 
 *  Title: real-time searchback procedure for the QRS detector
 *  Author:     Diego Sogari
 *  Modified:   12/June/2014
 *
 *=======================================================================*/
#ifndef SEARCHBACK
#define SEARCHBACK

#include <math.h>
#include "mex.h"

#define IDXOK(k)    (bufferLen+(k) > 0)
#define BUFVAL(k)   (buffer[bufferLen-1+(k)])

/*=========================================================================
 * Type definitions
 *=======================================================================*/
template <class type>
class SearchBack {
private:
    int findmax(const type *buffer, int start, mwSize len);
protected:
    mwSize bufferLen;
    int peakIdx;
    int searchbackIdx;
public:
    SearchBack(mwSize buflen);
    ~SearchBack();
    void newx(const type *buffer, int qrsIdx, type thresh, mwSize rrMean,
        mwSize rrMiss, bool sigRising, bool qrsConfirmed);
    int outputPeakIdx();
};

/*=========================================================================
 * Constructor
 *=======================================================================*/
template <class type>
SearchBack<type>::SearchBack(mwSize buflen)
{
    this->bufferLen = buflen;
    this->peakIdx = 0;
    this->searchbackIdx = 0;
}

/*=========================================================================
 * Destructor
 *=======================================================================*/
template <class type>
SearchBack<type>::~SearchBack()
{
}

/*=========================================================================
 * Find the position of the maximum value on the detection signal 
 *=======================================================================*/
template <class type>
int SearchBack<type>::findmax(const type *buffer, int start, mwSize len)
{
    mwSize i;
    type y, newy;
    int pos = 0;
    
    y = BUFVAL(start);
    for (i = 1; i < len; i++) {
        newy = BUFVAL(start + i);
        if (newy > y) {
            y = newy;
            pos = i;
        }
    }
    return start + pos;
}

/*=========================================================================
 * Update filter memory with an incoming sample
 *=======================================================================*/
template <class type>
void SearchBack<type>::newx(const type *buffer, int qrsIdx, type thresh,
        mwSize rrMean, mwSize rrMiss, bool sigRising, bool qrsConfirmed)
{
    peakIdx = 0;
    
    if (qrsConfirmed) {
        searchbackIdx = qrsIdx;
    }
    else if ((!sigRising || BUFVAL(0) < (thresh / (type)2)) &&
            (searchbackIdx != 0 && 0 >= searchbackIdx + rrMiss) &&
            IDXOK(1 - rrMean)) {
        // search back and locate the max in this interval
        peakIdx = findmax(buffer, 1 - rrMean, rrMean);
    }
    
    if (searchbackIdx != 0) {
        searchbackIdx--;
    }
}

/*=========================================================================
 * Return the output
 *=======================================================================*/
template <class type>
int SearchBack<type>::outputPeakIdx()
{
    return this->peakIdx;
}

#endif