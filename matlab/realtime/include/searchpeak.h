/*=========================================================================
 * searchpeak.h
 * 
 *  Title: real-time search of peaks for the FP detector
 *  Author:     Diego Sogari
 *  Modified:   12/June/2014
 *
 *=======================================================================*/
#ifndef SEARCHPEAK
#define SEARCHPEAK

#include <math.h>
#include "mex.h"

#define IDXOK(k)    (bufferLen+(k) > 0)
#define BUFDVAL(k)  (bufferD[bufferLen-1+(k)])
#define BUFMVAL(k)  (bufferM[bufferLen-1+(k)])

/*=========================================================================
 * Type definitions
 *=======================================================================*/
template <class type>
class SearchPeak {
private:
    bool check_indices(int start, int end, int *inc);
protected:
    mwSize bufferLen;
    int peakIdx;
public:
    SearchPeak(mwSize buflen);
    ~SearchPeak();
    void newx(const type *bufferD, const type *bufferM, int start,
            int end, int def, bool beatDetected);
    int outputPeakIdx();
};

/*=========================================================================
 * Constructor
 *=======================================================================*/
template <class type>
SearchPeak<type>::SearchPeak(mwSize buflen)
{
    this->bufferLen = buflen;
    this->peakIdx = 0;
}

/*=========================================================================
 * Destructor
 *=======================================================================*/
template <class type>
SearchPeak<type>::~SearchPeak()
{
}

/*=========================================================================
 * Check that the indices are correct for detection of fiducial points
 *=======================================================================*/
template <class type>
bool SearchPeak<type>::check_indices(int start, int end, int *inc)
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
void SearchPeak<type>::newx(const type *bufferD, const type *bufferM,
        int start, int end, int def, bool beatDetected)
{
    if (beatDetected) {
        int i, imax, idx[2];
        type y, ymax, val[2];
        int left, right, inc;

        if (!check_indices(start, end, &inc)) {
            peakIdx = def;
            return;
        }

        val[0] = val[1] = 0;
        idx[0] = idx[1] = def;
        for (i = start + inc; i != end - inc; i += inc) {
            y = BUFDVAL(i);
            if (BUFDVAL(i - 1) < y && y >= BUFDVAL(i + 1)) {
                if (y > val[0]) {
                    idx[0] = i;
                    val[0] = y;
                }
            }
            else if (BUFDVAL(i - 1) > y && y <= BUFDVAL(i + 1)) {
                if (y < val[1]) {
                    idx[1] = i;
                    val[1] = y;
                }
            }
        }

        if (idx[0] < idx[1]) {
            left = idx[0];
            right = idx[1];
        }
        else {
            left = idx[1];
            right = idx[0];
        }

        ymax = 0;
        imax = left;
        for (i = left + 1; i < right - 1; i++) {
            y = abs(BUFMVAL(i));
            if (abs(BUFMVAL(i - 1)) < y && y >= abs(BUFMVAL(i + 1))) {
                if (y > ymax) {
                    ymax = y;
                    imax = i;
                }
            }
        }

        peakIdx = imax;
    }
}

/*=========================================================================
 * Return the output
 *=======================================================================*/
template <class type>
int SearchPeak<type>::outputPeakIdx()
{
    return this->peakIdx;
}

#endif