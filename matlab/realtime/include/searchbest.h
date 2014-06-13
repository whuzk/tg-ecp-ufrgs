/*=========================================================================
 * searchbest.h
 * 
 *  Title: real-time search of best max/min for the FP detector
 *  Author:     Diego Sogari
 *  Modified:   12/June/2014
 *
 *=======================================================================*/
#ifndef SEARCHBEST
#define SEARCHBEST

#include <math.h>
#include "mex.h"

#define IDXOK(k)    (bufferLen+(k) > 0)
#define BUFVAL(k)  (buffer[bufferLen-1+(k)])

/*=========================================================================
 * Type definitions
 *=======================================================================*/
template <class type>
class SearchBest {
private:
    bool check_indices(int start, int end, int *inc);
protected:
    mwSize bufferLen;
    int peakIdx;
    bool ismax;
public:
    SearchBest(mwSize buflen, bool ismax);
    ~SearchBest();
    void newx(const type *buffer, int start, int end, int def);
    int outputPeakIdx();
};

/*=========================================================================
 * Constructor
 *=======================================================================*/
template <class type>
SearchBest<type>::SearchBest(mwSize buflen, bool ismax)
{
    this->bufferLen = buflen;
    this->peakIdx = 0;
    this->ismax = ismax;
}

/*=========================================================================
 * Destructor
 *=======================================================================*/
template <class type>
SearchBest<type>::~SearchBest()
{
}

/*=========================================================================
 * Check that the indices are correct for detection of fiducial points
 *=======================================================================*/
template <class type>
bool SearchBest<type>::check_indices(int start, int end, int *inc)
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
void SearchBest<type>::newx(const type *buffer, int start, int end, int def)
{
    int i, isave, isavepk, inc;
    type y, ysave, ysavepk;
    bool found;

    if (!check_indices(start, end, &inc)) {
        peakIdx = def;
        return;
    }

    found = false;
    ysave = ysavepk = 0;
    isave = isavepk = start;
    if (ismax) {
        for (i = start + inc; i != end - inc; i += inc) {
            y = BUFVAL(i);
            if (BUFVAL(i - 1) < y && y >= BUFVAL(i + 1)) {
                if (y > ysavepk) {
                    ysavepk = y;
                    isavepk = i;
                    found = true;
                }
            }
            else if (y > ysave) {
                ysave = y;
                isave = i;
            }
        }
    }
    else {
        for (i = start + inc; i != end - inc; i += inc) {
            y = BUFVAL(i);
            if (BUFVAL(i - 1) > y && y <= BUFVAL(i + 1)) {
                if (y < ysavepk) {
                    ysavepk = y;
                    isavepk = i;
                    found = true;
                }
            }
            else if (y < ysave) {
                ysave = y;
                isave = i;
            }
        }
    }

    if (!found) {
        peakIdx = isave;
    }
    else peakIdx = isavepk;
}

/*=========================================================================
 * Return the output
 *=======================================================================*/
template <class type>
int SearchBest<type>::outputPeakIdx()
{
    return this->peakIdx;
}

#endif