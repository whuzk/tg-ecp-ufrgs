/*=========================================================================
 * c_timeutils.h
 *=======================================================================*/
#ifndef C_TIMEUTILS
#define C_TIMEUTILS

#include <windows.h>
#include "mex.h"

/*=========================================================================
 * Global variables
 *=======================================================================*/
static double pcfreq = 0.0;

/*=========================================================================
 * Start timer
 *=======================================================================*/
long long tic()
{
    LARGE_INTEGER li;
    if (QueryPerformanceFrequency(&li)) {
        pcfreq = (double)li.QuadPart;
        if (QueryPerformanceCounter(&li)) {
            return li.QuadPart;
        }
        else {
            mexErrMsgIdAndTxt(
                    "EcgToolbox:c_timeutils:pcCounterNotAccessible",
                    "Could not read PC counter.");
            return -1;
        }
    }
    else {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_timeutils:pcFreqNotAccessible",
                "Could not read PC frequency.");
        pcfreq = 0.0;
    }
    
}

/*=========================================================================
 * Stop timer
 *=======================================================================*/
double toc(long long start)
{
    LARGE_INTEGER li;
    if (pcfreq == 0.0) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_timeutils:pcFreqNotSet",
                "PC frequency was not set.");
        return -1.0;
    }
    else if (start < 0) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_timeutils:startCountNotSet",
                "Start count was not set.");
        return -1.0;
    }
    else if (!QueryPerformanceCounter(&li)) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_timeutils:pcCounterNotAccessible",
                "Could not read PC counter.");
        return -1.0;
    }
    else return (li.QuadPart - start) / pcfreq;
}

#endif