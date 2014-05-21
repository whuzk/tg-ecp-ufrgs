/*=========================================================================
 * c_fiducial_marks.c
 * 
 *  Title: detection of fiducial marks
 *  Author:     Diego Sogari
 *  Modified:   07/May/2014
 *
 *  Intputs:
 *      1. Derivative of ECG signal
 *      2. Morphological derivative of ECG signal
 *      3. Estimated location of R peaks
 *      4. sampling frequency
 *
 *  Outputs:
 *      1. location of P-wave onset
 *      2. location of P-wave peak
 *      3. location of R-wave onset
 *      4. location of R-wave peak
 *      5. location of R-wave offset
 *      6. location of T-wave peak
 *      7. location of T-wave offset
 *
 *  *all inputs and outputs required
 *
 *=======================================================================*/
#include <math.h>
#include "mex.h"
#include "c_mathutils.h"
#include "c_timeutils.h"

/*=========================================================================
 * Abbreviations
 *  LEN     length
 *  BUF     buffer
 *  MIN     minimum
 *  MAX     maximum
 *  DEF     default
 *  SEC     seconds
 *  FREQ    frequency
 *  SIG     signal
 *  QRS     QRS complex
 *  RR      R-R interval
 *  IDX     index
 *  DE      derivative
 *  MD      morphological derivative
 *=======================================================================*/

/*=========================================================================
 * Constants
 *=======================================================================*/
#define MIN_INPUTS      4       // minimum number of input arguments
#define MAX_INPUTS      4       // maximum number of input arguments
#define MIN_OUTPUTS     7       // minimum number of output arguments
#define MAX_OUTPUTS     7       // maximum number of output arguments

#define SEARCH_L1_SEC   0.10    // search limit for R-wave points
#define SEARCH_L2_SEC   0.02    // search limit for Ronset and Roffset
#define SEARCH_L3_SEC   0.15    // search limit for T- and P-wave points
#define SEARCH_L4_SEC   0.08    // search limit for R onset and offset

/*=========================================================================
 * Input and output variables
 *=======================================================================*/
static double *inputSig1;       // input signal
static double *inputSig2;       // second input signal
static mwSize inputLen;         // input signal length
static double *qrsHist;         // estimated R peaks
static mwSize qrsLen;           // R peaks history length
static double sampFreq;         // sampling frequency
static double *ponHist;         // history of P-wave onsets
static double *ppkHist;         // history of P-wave peaks
static double *ronHist;         // history of R-wave onsets
static double *rpkHist;         // history of R-wave peaks
static double *rofHist;         // history of R-wave offsets
static double *tpkHist;         // history of T-wave peaks
static double *tofHist;         // history of T-wave offsets

/*=========================================================================
 * Buffer variables
 *=======================================================================*/
static unsigned int ci;             // current sample index
static mwSize bi;                   // current beat index
static int *deBuf;                  // buffer for the derivative signal
static int *mdBuf;                  // buffer for the MD signal
static mwSize bufLen;               // length of both buffers
static mwSize searchL1;             // search limit 1
static mwSize searchL2;             // search limit 2
static mwSize searchL3;             // search limit 3
static mwSize searchL4;             // search limit 4

/*=========================================================================
 * Fast lookup for detection buffers
 *=======================================================================*/
#define deb(I) (deBuf[(ci+(I))&(bufLen-1)])
#define mdb(I) (mdBuf[(ci+(I))&(bufLen-1)])

/*=========================================================================
 * Detection variables
 *=======================================================================*/
static mwSize Ponset;           // current P-wave onset
static mwSize Ppeak;            // current P-wave peak
static mwSize Ronset;           // current R-wave onset
static mwSize Rpeak;            // current R-wave peak
static mwSize Roffset;          // current R-wave offset
static mwSize Tpeak;            // current T-wave peak
static mwSize Toffset;          // current T-wave offset

/*=========================================================================
 * Check if index is valid in detection buffers
 *=======================================================================*/
bool isvalidindex(mwSize idx)
{
    return (idx <= 0 && 0 < bufLen + idx);
}

/*=========================================================================
 * Update outputs
 *=======================================================================*/
void updateOutputs()
{
    if (bi < qrsLen) {
        ponHist[bi] = ci + Ponset + 1;
        ppkHist[bi] = ci + Ppeak + 1;
        ronHist[bi] = ci + Ronset + 1;
        rpkHist[bi] = ci + Rpeak + 1;
        rofHist[bi] = ci + Roffset + 1;
        tpkHist[bi] = ci + Tpeak + 1;
        tofHist[bi] = ci + Toffset + 1;
        bi++;
    }
}

/*=========================================================================
 * Check that the indices are correct for detection of peaks
 *=======================================================================*/
bool check_indices(mwSize start, mwSize end, mwSize *inc)
{
    if (abs(end - start) <= 1 || !isvalidindex(start) || !isvalidindex(end)) {
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
 * Perform the detection of the Rpeak
 *=======================================================================*/
mwSize search_peak_abs(mwSize start, mwSize end, mwSize def)
{
    mwSize i, imax, idx[2];
    int y, ymax, val[2];
    mwSize left, right, inc;
    
    if (!check_indices(start, end, &inc)) {
        return def;
    }
    
    val[0] = val[1] = 0;
    idx[0] = idx[1] = def;
    for (i = start + inc; i != end - inc; i += inc) {
        y = deb(i);
        if (deb(i - 1) < y && y >= deb(i + 1)) {
            if (y > val[0]) {
                idx[0] = i;
                val[0] = y;
            }
        }
        else if (deb(i - 1) > y && y <= deb(i + 1)) {
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
        y = abs(mdb(i));
        if (abs(mdb(i - 1)) < y && y >= abs(mdb(i + 1))) {
            if (y > ymax) {
                ymax = y;
                imax = i;
            }
        }
    }
    
    return imax;
}

/*=========================================================================
 * Perform the detection of the first maximum peak on the MD signal
 *=======================================================================*/
mwSize search_first_mark_max(mwSize start, mwSize end, mwSize def)
{
    mwSize i, imax;
    int y, ymax;
    mwSize inc;
    bool found;
    
    if (!check_indices(start, end, &inc)) {
        return def;
    }
    
    ymax = 0;
    imax = start;
    found = false;
    i = start + inc;
    while (!found && i != end - inc) {
        y = mdb(i);
        if (mdb(i - 1) < y && y >= mdb(i + 1)) {
            found = true;
        }
        else {
            /*if (y > ymax) {
                ymax = y;
                imax = i;
            }*/
            i += inc;
        }
    }
    
    if (!found) {
        return def;
    }
    else return i;
}

/*=========================================================================
 * Perform the detection of the first minimum peak on the MD signal
 *=======================================================================*/
mwSize search_first_mark_min(mwSize start, mwSize end, mwSize def)
{
    mwSize i, imin;
    int y, ymin;
    mwSize inc;
    bool found;
    
    if (!check_indices(start, end, &inc)) {
        return def;
    }
    
    ymin = 0;
    imin = start;
    found = false;
    i = start + inc;
    while (!found && i != end - inc) {
        y = mdb(i);
        if (mdb(i - 1) > y && y <= mdb(i + 1)) {
            found = true;
        }
        else {
            /*if (y < ymin) {
                ymin = y;
                imin = i;
            }*/
            i += inc;
        }
    }
    
    if (!found) {
        return def;
    }
    else return i;
}

/*=========================================================================
 * Perform the detection of the best maximum peak on the absolute value of
 * the MD signal
 *=======================================================================*/
mwSize search_best_mark_abs(mwSize start, mwSize end, mwSize def)
{
    mwSize i, imax, imaxpk;
    int y, ymax, ymaxpk;
    mwSize inc;
    bool found;
    
    if (!check_indices(start, end, &inc)) {
        return def;
    }
    
    found = false;
    ymax = ymaxpk = 0;
    imax = imaxpk = start;
    for (i = start + inc; i != end - inc; i += inc) {
        y = abs(mdb(i));
        if (abs(mdb(i - 1)) < y && y >= abs(mdb(i + 1))) {
            if (y > ymaxpk) {
                ymaxpk = y;
                imaxpk = i;
                found = true;
            }
        }
        else if (y > ymax) {
            ymax = y;
            imax = i;
        }
    }
    
    if (!found) {
        return imax;
    }
    else return imaxpk;
}

/*=========================================================================
 * Perform the detection of the best maximum peak on the MD signal
 *=======================================================================*/
mwSize search_best_mark_max(mwSize start, mwSize end, mwSize def)
{
    mwSize i, imax, imaxpk;
    int y, ymax, ymaxpk;
    mwSize inc;
    bool found;
    
    if (!check_indices(start, end, &inc)) {
        return def;
    }
    
    found = false;
    ymax = ymaxpk = 0;
    imax = imaxpk = start;
    for (i = start + inc; i != end - inc; i += inc) {
        y = mdb(i);
        if (mdb(i - 1) < y && y >= mdb(i + 1)) {
            if (y > ymaxpk) {
                ymaxpk = y;
                imaxpk = i;
                found = true;
            }
        }
        else if (y > ymax) {
            ymax = y;
            imax = i;
        }
    }
    
    if (!found) {
        return imax;
    }
    else return imaxpk;
}

/*=========================================================================
 * Perform the detection of the best minimum peak on the MD signal
 *=======================================================================*/
mwSize search_best_mark_min(mwSize start, mwSize end, mwSize def)
{
    mwSize i, imin, iminpk;
    int y, ymin, yminpk;
    mwSize inc;
    bool found;
    
    if (!check_indices(start, end, &inc)) {
        return def;
    }
    
    found = false;
    ymin = yminpk = 0;
    imin = iminpk = start;
    for (i = start + inc; i != end - inc; i += inc) {
        y = mdb(i);
        if (mdb(i - 1) > y && y <= mdb(i + 1)) {
            if (y < yminpk) {
                yminpk = y;
                iminpk = i;
                found = true;
            }
        }
        else if (y < ymin) {
            ymin = y;
            imin = i;
        }
    }
    
    if (!found) {
        return imin;
    }
    else return iminpk;
}

/*=========================================================================
 * Special search
 *=======================================================================*/
mwSize special_search(mwSize test, mwSize def)
{
    int y1 = mdb(test);
    int y2 = mdb(def);
    if (SIGN(y1) != SIGN(y2) && abs(y1) >= abs(y2 >> 2)) {
        return test;
    }
    else return def;
}

/*=========================================================================
 * Perform the detection of fiducial marks
 *=======================================================================*/
void detectFiducialMarks(mwSize r, mwSize rr1, mwSize rr2)
{
    mwSize len, temp1, temp2;
    
    // R-wave
    Rpeak = search_peak_abs(r - searchL1, r + searchL1, r);
    if (mdb(Rpeak) > 0) {
        // inverted
        Ronset = search_first_mark_min(Rpeak - searchL2, Rpeak - searchL1, Rpeak);
        Roffset = search_first_mark_min(Rpeak + searchL2, Rpeak + searchL1, Rpeak);
        temp1 = search_first_mark_max(Ronset - 1, Ronset - searchL4, Ronset);
        temp2 = search_first_mark_max(Roffset + 1, Roffset + searchL4, Roffset);
    }
    else {
        // normal
        Ronset = search_first_mark_max(Rpeak - searchL2, Rpeak - searchL1, Rpeak);
        Roffset = search_first_mark_max(Rpeak + searchL2, Rpeak + searchL1, Rpeak);
        temp1 = search_first_mark_min(Ronset - 1, Ronset - searchL4, Ronset);
        temp2 = search_first_mark_min(Roffset + 1, Roffset + searchL4, Roffset);
    }
    Ronset = special_search(temp1, Ronset);
    Roffset = special_search(temp2, Roffset);
    
    // P-wave
    len = (rr1 + (rr1 << 1)) >> 3;
    Ppeak = search_peak_abs(Ronset - 1, Ronset - len, Ronset);
    if (mdb(Ppeak) > 0) {
        // inverted
        Ponset = search_best_mark_min(Ppeak - 1, Ppeak - searchL3, Ppeak);
    }
    else {
        // normal
        Ponset = search_best_mark_max(Ppeak - 1, Ppeak - searchL3, Ppeak);
    }
    
    // T-wave
    len = rr2 >> 1;
    Tpeak = search_peak_abs(Roffset + 1, Roffset + len, Roffset);
    if (mdb(Tpeak) > 0) {
        // inverted
        Toffset = search_best_mark_min(Tpeak + 1, Tpeak + searchL3, Tpeak);
    }
    else {
        // normal
        Toffset = search_best_mark_max(Tpeak + 1, Tpeak + searchL3, Tpeak);
    }
}

/*=========================================================================
 * Event triggered for an arriving sample
 *=======================================================================*/
void onNewSample(double sample1, double sample2, bool newQrs)
{
    static mwSize countDown = -1;
    static mwSize rr = -1;
    static mwSize qrs;
    
    // update detection buffers
    deb(0) = (int)sample1;
    mdb(0) = (int)sample2;
    
    // simulate QRS detection
    if ((newQrs && countDown > 0) || countDown == 0 || ci == inputLen - 1) {
        
        // perform fiducial mark detection
        detectFiducialMarks(qrs, rr, 0 - qrs);
        // update the outputs
        updateOutputs();
    }
    
    if (newQrs) {
        // update rr
        if (rr == -1) {
            rr = (mwSize)sampFreq;
        }
        else rr = 0 - qrs;
        // update qrs
        qrs = 0;
        // update countdown
        countDown = rr;
    }
    else if (countDown >= 0) {
        // update countdown
        countDown--;
    }
    
    // decrement qrs index
    qrs--;
    
    // increment sample index
    ci++;
}

/*=========================================================================
 * Check correct number of arguments and their types 
 *=======================================================================*/
void checkArgs( int nlhs, mxArray *plhs[],
                int nrhs, const mxArray *prhs[])
{
    mwSize i;
    
    // check for proper number of input arguments
    if (nrhs < MIN_INPUTS || nrhs > MAX_INPUTS) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_fiducial_marks:nrhs",
            "%d input(s) required and %d optional",
            MIN_INPUTS, MAX_INPUTS - MIN_INPUTS);
    }
    // check for proper number of output arguments
    if (nlhs < MIN_OUTPUTS || nlhs > MAX_OUTPUTS) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_fiducial_marks:nlhs",
            "%d output(s) required and %d optional",
            MIN_OUTPUTS, MAX_OUTPUTS - MIN_OUTPUTS);
    }
    // make sure all input arguments are of type double
    for (i = 0; i < nrhs; i++) {
        if (!mxIsDouble(prhs[i]) || mxIsComplex(prhs[i])) {
            mexErrMsgIdAndTxt(
                "EcgToolbox:c_fiducial_marks:notDouble",
                "Input #%d must be of type double.", i + 1);
        }
    }
    // make sure the first four arguments are vectors
    for (i = 0; i < 3; i++) {
        if (mxGetM(prhs[i]) != 1 && mxGetN(prhs[i]) != 1) {
            mexErrMsgIdAndTxt(
                "EcgToolbox:c_fiducial_marks:notVector",
                "Input #%d must be a vector.", i + 1);
        }
    }
    // make sure the remaining input arguments are all scalars
    for (i = 3; i < nrhs; i++) {
        if (mxGetNumberOfElements(prhs[i]) != 1) {
            mexErrMsgIdAndTxt(
                "EcgToolbox:c_fiducial_marks:notScalar",
                "Input #%d must be a scalar.", i + 1);
        }
    }
}

/*=========================================================================
 * Hanlde input arguments 
 *=======================================================================*/
void handleInputs( int nrhs, const mxArray *prhs[],
                   mwSize *nrows1, mwSize *ncols1,
                   mwSize *nrows2, mwSize *ncols2)
{
    // get pointers to the data in the input vectors
    inputSig1 = mxGetPr(prhs[0]);
    inputSig2 = mxGetPr(prhs[1]);
    qrsHist = mxGetPr(prhs[2]);
    
    // get the dimensions of the first input vector
    *nrows1 = (mwSize)mxGetM(prhs[0]);
    *ncols1 = (mwSize)mxGetN(prhs[0]);
    
    // make sure the dimensions of the second input vector is compatible
    if ((mwSize)mxGetM(prhs[1]) != *nrows1 || (mwSize)mxGetN(prhs[1]) != *ncols1) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_fiducial_marks:badDimensions",
            "The first and second inputs must have the same dimensions.");
    }
    
    // get the dimensions of the third input vector
    *nrows2 = (mwSize)mxGetM(prhs[2]);
    *ncols2 = (mwSize)mxGetN(prhs[2]);
    
    // get the sampling frequency
    sampFreq = mxGetScalar(prhs[3]);
}

/*=========================================================================
 * Hanlde output arguments 
 *=======================================================================*/
void handleOutputs( int nlhs, mxArray *plhs[],
                    mwSize nrows, mwSize ncols)
{
    double *outVectors[MAX_OUTPUTS];
    mwSize i;
    
    // create the output vectors
    for (i = 0; i < nlhs; i++) {
        plhs[i] = mxCreateDoubleMatrix(nrows, ncols, mxREAL);
        outVectors[i] = mxGetPr(plhs[i]);
    }
    
    // get pointers to the output vectors
    ponHist = outVectors[0];
    ppkHist = outVectors[1];
    ronHist = outVectors[2];
    rpkHist = outVectors[3];
    rofHist = outVectors[4];
    tpkHist = outVectors[5];
    tofHist = outVectors[6];
}

/*=========================================================================
 * The initializarion routine
 *=======================================================================*/
void init()
{
    // initialize indices
    ci = 0;
    bi = 0;
    
    // predefined limits
    searchL1 = (mwSize)(SEARCH_L1_SEC * sampFreq);
    searchL2 = (mwSize)(SEARCH_L2_SEC * sampFreq);
    searchL3 = (mwSize)(SEARCH_L3_SEC * sampFreq);
    searchL4 = (mwSize)(SEARCH_L4_SEC * sampFreq);
    
    // create buffers
    bufLen = 1 << (1 + ILOG2(4 * sampFreq - 1));
    deBuf = (int *)mxMalloc(bufLen * sizeof(int));
    mdBuf = (int *)mxMalloc(bufLen * sizeof(int));
}

/*=========================================================================
 * The main routine 
 *=======================================================================*/
void doTheJob()
{
    double time;
    mwSize i, newqrs;
    mwSize beatCount = 0;
    
    // start time counter
    tic();
    
    // initialize qrs variables
    newqrs = (mwSize)qrsHist[beatCount] - 1;
    
    // process one input sample at a time
    for (i = 0; i < inputLen; i++) {
        onNewSample(inputSig1[i], inputSig2[i], newqrs == i);
        if (newqrs == i && beatCount < qrsLen - 1) {
            newqrs = (mwSize)qrsHist[++beatCount] - 1;
        }
    }
    
    // stop time counter
    time = toc();
    
    // display time statistics
    mexPrintf("Total processing time: %.2f ms\n", 1000 * time);
    mexPrintf("Average time per sample: %.2f ns\n", 1000000000 * time / inputLen);
}

/*=========================================================================
 * The finalization routine 
 *=======================================================================*/
void finalize()
{
    // deallocate memory
    mxFree(deBuf);
    mxFree(mdBuf);
}

/*=========================================================================
 * The gateway function 
 *=======================================================================*/
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    mwSize nrows1;
    mwSize ncols1;
    mwSize nrows2;
    mwSize ncols2;
    
    // check argument correctness
    checkArgs(nlhs, plhs, nrhs, prhs);
    
    // handle input arguments
    handleInputs(nrhs, prhs, &nrows1, &ncols1, &nrows2, &ncols2);
    
    // handle output arguments
    handleOutputs(nlhs, plhs, nrows2, ncols2);
    
    // calculate the length of the signal
    inputLen = max(nrows1,ncols1);
    
    // calculate the number of qrs
    qrsLen = max(nrows2,ncols2);
    
    // make some initializations
    init();
    
    // do the actual processing
    doTheJob();
    
    // perform final adjustments
    finalize();
}
