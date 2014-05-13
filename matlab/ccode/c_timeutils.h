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
static __int64 start;
static double pcfreq = 0.0;

/*=========================================================================
 * Start timer
 *=======================================================================*/
void tic()
{
    LARGE_INTEGER li;
    if (QueryPerformanceFrequency(&li)) {
        pcfreq = (double)li.QuadPart;
        if (QueryPerformanceCounter(&li)) {
            start = li.QuadPart;
        }
        else {
            mexErrMsgIdAndTxt(
                    "EcgToolbox:c_ecg_utils:pcCounterNotAccessible",
                    "Could not read PC counter.");
            start = -1;
        }
    }
    else {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_ecg_utils:pcFreqNotAccessible",
                "Could not read PC frequency.");
        pcfreq = 0.0;
    }
    
}

/*=========================================================================
 * Stop timer
 *=======================================================================*/
double toc()
{
    LARGE_INTEGER li;
    if (pcfreq == 0.0) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_ecg_utils:pcFreqNotSet",
                "PC frequency was not set.");
        return -1.0;
    }
    else if (start < 0) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_ecg_utils:startCountNotSet",
                "Start count was not set.");
        return -1.0;
    }
    else if (!QueryPerformanceCounter(&li)) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_ecg_utils:pcCounterNotAccessible",
                "Could not read PC counter.");
        return -1.0;
    }
    else return (li.QuadPart - start) / pcfreq;
}

#endif