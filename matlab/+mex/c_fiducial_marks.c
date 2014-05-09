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
 *      4. Estimated RR intervals
 *      2. sampling frequency
 *
 *  Outputs:
 *      1. location of P-wave onset
 *      2. location of P-wave peak
 *      3. location of R-wave onset
 *      4. location of R-wave peak
 *      5. location of R-wave offset
 *      6. location of T-wave peak
 *      7. location of T-wave offset
 *      8. location of isoeletric point
 *      9. location of J point
 *
 *  *all inputs and outputs required
 *
 *=======================================================================*/
#include <math.h>
#include "mex.h"
#include "c_mexutils.h"
#include "c_timeutils.h"

/*=========================================================================
 * Abbreviations
 *  LEN     length
 *  BUF     buffer
 *  MIN     minimum
 *  MAX     maximum
 *  DEF     default
 *  SEC     seconds
 *  DEL     delay
 *  WIN     window
 *  FREQ    frequency
 *  SIG     signal
 *  QRS     QRS complex
 *  RR      R-R interval
 *  AMP     amplitude
 *  IDX     index
 *  DE      derivative
 *  MD      morphological derivative
 *=======================================================================*/

/*=========================================================================
 * Constants
 *=======================================================================*/
#define MIN_INPUTS      5       // minimum number of input arguments
#define MAX_INPUTS      5       // maximum number of input arguments
#define MIN_OUTPUTS     9       // minimum number of output arguments
#define MAX_OUTPUTS     9       // maximum number of output arguments

/*=========================================================================
 * Input and output variables
 *=======================================================================*/
static double *inputSig1;       // input signal
static double *inputSig2;       // second input signal
static mwSize inputLen;         // input signal length
static double *qrsHist;         // estimated R peaks
static double *rrHist;          // estimated RR intervals
static mwSize qrsLen;           // R peaks history length
static double sampFreq;         // sampling frequency
static double *ponHist;         // history of P-wave onsets
static double *ppkHist;         // history of P-wave peaks
static double *ronHist;         // history of R-wave onsets
static double *rpkHist;         // history of R-wave peaks
static double *rofHist;         // history of R-wave offsets
static double *tpkHist;         // history of T-wave peaks
static double *tofHist;         // history of T-wave offsets
static double *iptHist;         // history of isoeletrict points
static double *jptHist;         // history of J points

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
static mwSize Ipoint;           // current isoeletrict point
static mwSize Jpoint;           // current J point

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
    ponHist[bi] = ci + Ponset + 1;
    ppkHist[bi] = ci + Ppeak + 1;
    ronHist[bi] = ci + Ronset + 1;
    rpkHist[bi] = ci + Rpeak + 1;
    rofHist[bi] = ci + Roffset + 1;
    tpkHist[bi] = ci + Tpeak + 1;
    tofHist[bi] = ci + Toffset + 1;
    iptHist[bi] = ci + Ipoint + 1;
    jptHist[bi] = ci + Jpoint + 1;
}

/*=========================================================================
 * Perform the detection of the Rpeak
 *=======================================================================*/
mwSize search_rpeak_abs(mwSize start, mwSize end, mwSize def)
{
    mwSize imax, idx[2];
    int y, ymax, val[2];
    mwSize left, right, j;
    
    if (!isvalidindex(start) || !isvalidindex(end)) {
        return def;
    }
    
    j = 1;
    val[0] = val[1] = 0;
    idx[0] = idx[1] = def;
    for (mwSize i = start + 1; i < end - 1; i++) {
        y = abs(deb(i));
        if (abs(deb(i - 1)) < y && y >= abs(deb(i + 1))) {
            if (y > val[j]) {
                idx[j] = i;
                val[j] = y;
                j = (j + 1) & 1;
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
    for (mwSize i = left; i < right; i++) {
        y = abs(mdb(i));
        if (y > ymax) {
            ymax = y;
            imax = i;
        }
    }
    return imax;
}

/*=========================================================================
 * Perform the detection of the ?
 *=======================================================================*/
mwSize search_first_mark_max1(mwSize start, mwSize inc, mwSize end, mwSize def)
{
    mwSize i;
    int y;
    bool found;
    
    if (!isvalidindex(start) || !isvalidindex(end) ||
            abs(end - start - inc) == 0) {
        return def;
    }
    
    found = false;
    i = start + inc;
    while (!found && i != end - inc) {
        y = deb(i);
        if (deb(i - 1) < y && y >= deb(i + 1)) {
            found = true;
        }
        else i += inc;
    }

    if (!found) {
        return def;
    }
    else return i;
}

/*=========================================================================
 * Perform the detection of the ?
 *=======================================================================*/
mwSize search_first_mark_min1(mwSize start, mwSize inc, mwSize end, mwSize def)
{
    mwSize i;
    int y;
    bool found;
    
    if (!isvalidindex(start) || !isvalidindex(end) ||
            abs(end - start - inc) == 0) {
        return def;
    }
    
    found = false;
    i = start + inc;
    while (!found && i != end - inc) {
        y = deb(i);
        if (deb(i - 1) > y && y <= deb(i + 1)) {
            found = true;
        }
        else i += inc;
    }
    
    if (!found) {
        return def;
    }
    else return i;
}

/*=========================================================================
 * Perform the detection of the ?
 *=======================================================================*/
mwSize search_first_mark_max2(mwSize start, mwSize inc, mwSize end, mwSize def)
{
    mwSize i, imax;
    int y, ymax;
    bool found;
    
    if (!isvalidindex(start) || !isvalidindex(end) ||
            abs(end - start - inc) == 0) {
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
            if (y > ymax) {
                ymax = y;
                imax = i;
            }
            i += inc;
        }
    }
    
    if (!found) {
        return imax;
    }
    else return i;
}

/*=========================================================================
 * Perform the detection of the ?
 *=======================================================================*/
mwSize search_first_mark_min2(mwSize start, mwSize inc, mwSize end, mwSize def)
{
    mwSize i, imin;
    int y, ymin;
    bool found;
    
    if (!isvalidindex(start) || !isvalidindex(end) ||
            abs(end - start - inc) == 0) {
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
            if (y < ymin) {
                ymin = y;
                imin = i;
            }
            i += inc;
        }
    }
    
    if (!found) {
        return imin;
    }
    else return i;
}

/*=========================================================================
 * Perform the detection of the ?
 *=======================================================================*/
mwSize search_best_mark_abs(mwSize start, mwSize inc, mwSize end, mwSize def)
{
    mwSize imax, imaxpk;
    int y, ymax, ymaxpk;
    bool found;
    
    if (!isvalidindex(start) || !isvalidindex(end) ||
            abs(end - start - inc) == 0) {
        return def;
    }
    
    found = false;
    ymax = ymaxpk = 0;
    imax = imaxpk = start;
    for (mwSize i = start + inc; i != end - inc; i += inc) {
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
 * Perform the detection of the ?
 *=======================================================================*/
mwSize search_best_mark_max(mwSize start, mwSize inc, mwSize end, mwSize def)
{
    mwSize imax, imaxpk;
    int y, ymax, ymaxpk;
    bool found;
    
    if (!isvalidindex(start) || !isvalidindex(end) ||
            abs(end - start - inc) == 0) {
        return def;
    }
    
    found = false;
    ymax = ymaxpk = 0;
    imax = imaxpk = start;
    for (mwSize i = start + inc; i != end - inc; i += inc) {
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
 * Perform the detection of the ?
 *=======================================================================*/
mwSize search_best_mark_min(mwSize start, mwSize inc, mwSize end, mwSize def)
{
    mwSize imin, iminpk;
    int y, ymin, yminpk;
    bool found;
    
    if (!isvalidindex(start) || !isvalidindex(end) ||
            abs(end - start - inc) == 0) {
        return def;
    }
    
    found = false;
    ymin = yminpk = 0;
    imin = iminpk = start;
    for (mwSize i = start + inc; i != end - inc; i += inc) {
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
 * Perform the detection of fiducial marks
 *=======================================================================*/
void detectFiducialMarks(mwSize r, mwSize rr)
{
    mwSize len;
    
    // R-wave
    Rpeak = search_rpeak_abs(r - searchL1, r + searchL1, r);
    if (mdb(Rpeak) > 0) {
        // inverted
        Ronset = search_first_mark_min2(Rpeak - searchL4, -1, Rpeak - searchL1, Rpeak);
        Roffset = search_first_mark_min2(Rpeak + searchL4, 1, Rpeak + searchL1, Rpeak);
        Ipoint = search_first_mark_max1(Ronset - 1, -1, Ronset - searchL2, Ronset);
        Jpoint = search_first_mark_min1(Roffset + 1, 1, Roffset + searchL2, Roffset);
    }
    else {
        // normal
        Ronset = search_first_mark_max2(Rpeak - searchL4, -1, Rpeak - searchL1, Rpeak);
        Roffset = search_first_mark_max2(Rpeak + searchL4, 1, Rpeak + searchL1, Rpeak);
        Ipoint = search_first_mark_min1(Ronset - 1, -1, Ronset - searchL2, Ronset);
        Jpoint = search_first_mark_max1(Roffset + 1, 1, Roffset + searchL2, Roffset);
    }
    
    // P-wave
    len = rr >> 2;
    Ppeak = search_best_mark_abs(Ronset - searchL2, -1, Ronset - len, Ronset);
    if (mdb(Ppeak) > 0) {
        // inverted
        Ponset = search_best_mark_min(Ppeak - searchL2, -1, Ppeak - searchL3, Ppeak);
    }
    else {
        // normal
        Ponset = search_best_mark_max(Ppeak - searchL2, -1, Ppeak - searchL3, Ppeak);
    }
    
    // T-wave
    len += rr >> 3;
    Tpeak = search_best_mark_abs(Roffset + searchL2, 1, Roffset + len, Roffset);
    if (mdb(Tpeak) > 0) {
        // inverted
        Toffset = search_best_mark_min(Tpeak + searchL2, 1, Tpeak + searchL3, Tpeak);
    }
    else {
        // normal
        Toffset = search_best_mark_max(Tpeak + searchL2, 1, Tpeak + searchL3, Tpeak);
    }
}

/*=========================================================================
 * Event triggered for an arriving sample
 *=======================================================================*/
void onNewSample(double sample1, double sample2)
{
    // update detection buffers
    deb(0) = (int)sample1;
    mdb(0) = (int)sample2;
    
    // simulate QRS detection
    if (bi < qrsLen) {
        mwSize qrs = (mwSize)qrsHist[bi];
        mwSize rr = (mwSize)rrHist[bi];
        if (ci > qrs + rr || ci == inputLen - 1) {
            // perform fiducial mark detection
            detectFiducialMarks(-rr, rr);
            // update the outputs
            updateOutputs();
            // increment beat index
            bi++;
        }
    }
    
    // increment sample index
    ci++;
}

/*=========================================================================
 * Check correct number of arguments and their types 
 *=======================================================================*/
void checkArgs( int nlhs, mxArray *plhs[],
                int nrhs, const mxArray *prhs[])
{
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
    // check the first arguments
    for (mwSize i = 0; i < 4; i++) {
        // make sure the input argument is a vector
        if (mxGetM(prhs[i]) != 1 && mxGetN(prhs[i]) != 1) {
            mexErrMsgIdAndTxt(
                "EcgToolbox:c_fiducial_marks:notVector",
                "Input #%d must be a vector.", i + 1);
        }
        // make sure the input argument is of type double
        if (!mxIsDouble(prhs[i]) || mxIsComplex(prhs[i])) {
            mexErrMsgIdAndTxt(
                "EcgToolbox:c_fiducial_marks:notDouble",
                "Input #%d must be of type double.", i + 1);
        }
    }
    // make sure the remaining arguments are all scalars
    for (mwSize i = 4; i < nrhs; i++) {
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
    rrHist = mxGetPr(prhs[3]);
    
    // get the dimensions of the first input vector
    *nrows1 = (mwSize)mxGetM(prhs[0]);
    *ncols1 = (mwSize)mxGetN(prhs[0]);
    
    // make sure the dimensions of the second input vector is compatible
    if ((mwSize)mxGetM(prhs[0]) != *nrows1 || (mwSize)mxGetN(prhs[0]) != *ncols1) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_fiducial_marks:badDimensions",
            "The first and second inputs must have the same dimensions.");
    }
    
    // get the dimensions of the third input vector
    *nrows2 = (mwSize)mxGetM(prhs[2]);
    *ncols2 = (mwSize)mxGetN(prhs[2]);
    
    // make sure the dimensions of the fourth input vector is compatible
    if ((mwSize)mxGetM(prhs[2]) != *nrows2 || (mwSize)mxGetN(prhs[2]) != *ncols2) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_fiducial_marks:badDimensions",
            "The third and fourth inputs must have the same dimensions.");
    }
    
    // get the sampling frequency
    sampFreq = (int)mxGetScalar(prhs[4]);
}

/*=========================================================================
 * Hanlde output arguments 
 *=======================================================================*/
void handleOutputs( int nlhs, mxArray *plhs[],
                    mwSize nrows, mwSize ncols)
{
    double *outVectors[MAX_OUTPUTS];
    
    // create the output vectors
    for (mwSize i = 0; i < nlhs; i++) {
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
    iptHist = outVectors[7];
    jptHist = outVectors[8];
}

/*=========================================================================
 * The initializarion routine
 *=======================================================================*/
void init()
{
    // initialize sample index
    ci = 0;
    
    // predefined limits
    searchL1 = (mwSize)round(0.10 * sampFreq);
    searchL2 = (mwSize)round(0.04 * sampFreq);
    searchL3 = (mwSize)round(0.20 * sampFreq);
    searchL4 = (mwSize)round(0.02 * sampFreq);
    
    // create buffers
    bufLen = 1 << (1 + ilogb(2 * sampFreq - 1));
    deBuf = (int *)mxMalloc(bufLen * sizeof(int));
    mdBuf = (int *)mxMalloc(bufLen * sizeof(int));
}

/*=========================================================================
 * The main routine 
 *=======================================================================*/
void doTheJob()
{
    double time;
    
    // start time counter
    tic();
    
    // process one input sample at a time
    for (mwSize i = 0; i < inputLen; i++) {
        onNewSample(inputSig1[i], inputSig2[i]);
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
    
    // calculate the length of the signal
    qrsLen = max(nrows2,ncols2);
    
    // make some initializations
    init();
    
    // do the actual processing
    doTheJob();
    
    // perform final adjustments
    finalize();
}
