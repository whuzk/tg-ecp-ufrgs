/*=========================================================================
 * searchfirst.h
 * 
 *  Title: real-time search of first max/min for the FP detector
 *  Author:     Diego Sogari
 *  Modified:   12/June/2014
 *
 *=======================================================================*/
#ifndef SEARCHFIRST
#define SEARCHFIRST

#include <math.h>
#include "mex.h"

#define IDXOK(k)    (bufferLen+(k) > 0)
#define BUFMVAL(k)  (bufferM[bufferLen-1+(k)])

/*=========================================================================
 * Type definitions
 *=======================================================================*/
template <class type>
class SearchFirst {
private:
    bool check_indices(int start, int end, int *inc);
protected:
    mwSize bufferLen;
    int peakIdx;
    bool ismax;
public:
    SearchFirst(mwSize buflen, bool ismax);
    ~SearchFirst();
    void newx(const type *bufferM, int start, int end, int def,
            bool beatDetected);
    int outputPeakIdx();
};

/*=========================================================================
 * Constructor
 *=======================================================================*/
template <class type>
SearchFirst<type>::SearchFirst(mwSize buflen, bool ismax)
{
    this->bufferLen = buflen;
    this->peakIdx = 0;
    this->ismax = ismax;
}

/*=========================================================================
 * Destructor
 *=======================================================================*/
template <class type>
SearchFirst<type>::~SearchFirst()
{
}

/*=========================================================================
 * Check that the indices are correct for detection of fiducial points
 *=======================================================================*/
template <class type>
bool SearchFirst<type>::check_indices(int start, int end, int *inc)
{
    if (abs(end - start) <= 1 || !IDXOK(start) || !IDXOK(end)) {
        return false;
    }
    else if (start < end) {
        *inc = 1;
        return true;
    }
    else {
        *inc = -1;
        return true;
    }
}

/*=========================================================================
 * Update filter memory with an incoming sample
 *=======================================================================*/
template <class type>
void SearchFirst<type>::newx(const type *bufferM, int start, int end,
        int def, bool beatDetected)
{
    if (beatDetected) {
        int i, inc;
        type y;

        if (!check_indices(start, end, &inc)) {
            peakIdx = def;
            return;
        }
        if (ismax) {
            for (i = start + inc; i != end - inc; i += inc) {
                y = BUFMVAL(i);
                if (BUFMVAL(i - 1) < y && y >= BUFMVAL(i + 1)) {
                    break;
                }
            }
        }
        else {
            for (i = start + inc; i != end - inc; i += inc) {
                y = BUFMVAL(i);
                if (BUFMVAL(i - 1) > y && y <= BUFMVAL(i + 1)) {
                    break;
                }
            }
        }
        peakIdx = i;
    }
}

/*=========================================================================
 * Return the output
 *=======================================================================*/
template <class type>
int SearchFirst<type>::outputPeakIdx()
{
    return this->peakIdx;
}

#endif