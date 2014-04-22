/*=========================================================================
 * c_filter_int.c
 *=======================================================================*/
#include <math.h>
#include "mex.h"
#include "c_ecg_utils.h"

#define MIN_INPUTS  2
#define MAX_INPUTS  3
#define MIN_OUTPUTS 1
#define MAX_OUTPUTS 2

#define DEFAULT_M   3
#define GAIN_LEN_MS 50
#define MAFH_LEN_MS 100

// low-pass filter and second-order backward difference
static short *lpfhVec;
static mwSize lpfhLen;
static double lpfhDelay;
static mwSize lpfhHalfLen;

// windowed-max filter
static short *wmfVec;
static unsigned int *wmfAux;
static mwSize wmfLen;
static double wmfDelay;

// automatic gain control
static short *agcVec;
static mwSize agcLen;
static double agcDelay;
static short agcMax;
static short agcMin;
static short agcLog2Max;

// multiplication of absolute backward differences
static int *mobdVec;
static mwSize mobdLen;
static double mobdDelay;

// moving-average filter
static int *mafhVec;
static mwSize mafhLen;
static double mafhDelay;
static mwSize mafhLog2Len;

/* Product of elements */
int prod(int *x, mwSize nx)
{
    int accp = 1;
    for (mwSize k = 0; k < nx; k++) {
        accp *= x[k];
    }
    return accp;
}

/* Low-pass filter and second-order backward difference */
short lpf(short sample)
{
    static mwSize n = 0;
    mwSize n2;
    int y0;
    
    // calculate index at halfway
    n2 = (n + lpfhHalfLen) % lpfhLen;
    
    // update filter output
    y0 = sample - (lpfhVec[n2] << 1) + lpfhVec[n];
    
    // update filter memory
    lpfhVec[n] = sample;
    if (++n == lpfhLen)
        n = 0;
    
    // compute result
    return y0 >> 2;
}

/* Windowed-max filter */
short wmf(short sample)
{
    static mwSize first = 0;
    static mwSize count = 0;
    static unsigned int i = 0;
    mwSize j,k,idx;
    
    // search for the first element greater than the current sample
    j = count;
    k = first + count - 1;
    if (k == wmfLen) {
        k = 0;
    }
    while (j > 0 && sample >= wmfVec[k]) {
        j--;
        if (--k < 0) {
            k = wmfLen - 1;
        }
    }
    // get the index next to element found
    idx = k + 1;
    if (idx == wmfLen) {
        idx = 0;
    }
    // put the sample there and adjust the length
    wmfVec[idx] = sample;
    wmfAux[idx] = i;
    count = j + 1;
    // check if the first in line has gone out of the window length
    if (count > wmfLen || wmfAux[first] == i - wmfLen) {
        if (++first == wmfLen)
            first = 0;
        count--;
    }
    // increment global index
    i++;
    // return the max in the window
    return wmfVec[first];
}

/* Automatic gain control */
short agc(short sample, short gain)
{
    static mwSize n = 0;
    int y0;
    
    // update filter output
    gain = max(gain, agcMin);
    y0 = (agcVec[n] * agcMax) / gain;
    //y0 = agcVec[n] << agcLog2Max;
    //while (gain >>= 1) y0 >>= 1;
    
    // update filter memory
    agcVec[n] = sample;
    if (++n == agcLen)
        n = 0;
    
    // compute result
    return (short)min(y0, agcMax);
}

/* Multiplication of absolute backward differences */
int modb(short sample)
{
    static mwSize n = 0;
    
    // update filter memory
    mobdVec[n] = sample;
    if (++n == mobdLen)
        n = 0;
    
    // compute result
    return prod(mobdVec, mobdLen);
}

/* Moving-average filter of length 32 */
int maf(int sample)
{
    static mwSize n = 0;
    static long long y0 = 0;
    
    // update filter output
    y0 += sample - mafhVec[n];
    
    // update filter memory
    mafhVec[n] = sample;
    if (++n == mafhLen)
        n = 0;
    
    // compute result
    //return (int)min(y0 / mafhLen, INT_MAX);
    return (int)min(y0 >> mafhLog2Len, INT_MAX);
}

/* Process one intput sample */
int processSample(short sample)
{
    sample = abs(lpf(sample));
    sample = agc(sample, wmf(sample));
    return maf(modb(sample));
}

/* Initialize global variables */
void initVariables(int Fs, int m)
{
    // low-pass filter and second-order backward difference
    lpfhLen = 12;
    lpfhVec = (short*)mxCalloc(lpfhLen, sizeof(short));
    lpfhDelay = lpfhLen / 2.0;
    lpfhHalfLen = lpfhLen >> 1;
    
    // windowed-max filter
    wmfLen = 2 * Fs;
    wmfVec = (short*)mxMalloc(wmfLen * sizeof(short));
    wmfAux = (unsigned int*)mxMalloc(wmfLen * sizeof(int));
    wmfDelay = 0;

    // automatic gain control
    agcLen = (int)(Fs * GAIN_LEN_MS * 0.001);
    agcVec = (short*)mxCalloc(agcLen, sizeof(short));
    agcDelay = agcLen;
    agcLog2Max = (short)((sizeof(int) << 3) - 1) / m;
    agcMax = (1 << agcLog2Max) - 1;
    agcMin = 1 << 3;

    // multiplication of absolute backward differences
    mobdLen = m;
    mobdVec = (int*)mxMalloc(mobdLen * sizeof(int));
    mobdDelay = (mobdLen - 1) / 2.0;
    for (mwSize k = 0; k < mobdLen; k++) {
        mobdVec[k] = 1;
    }
    
    // moving-average filter
    mafhLen = (int)(Fs * MAFH_LEN_MS * 0.0005) * 2 + 1;
    mafhVec = (int*)mxCalloc(mafhLen, sizeof(int));
    mafhDelay = (mafhLen - 1) / 2.0;
    mafhLog2Len = (int)ceil(log2(mafhLen));
}

/* Finalize global variables */
void finVariables()
{
    mxFree(lpfhVec);
    mxFree(wmfVec);
    mxFree(wmfAux);
    mxFree(agcVec);
    mxFree(mobdVec);
    mxFree(mafhVec);
}

/* Process input and output arguments */
void processArgs( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[],
                  mwSize *nrows, mwSize *ncols,
                  int *sampFreq, int *mFactor)
{
    /* check for proper number of arguments */
    if (nrhs < MIN_INPUTS || nrhs > MAX_INPUTS) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_filter_int:nrhs",
                "Two inputs required.");
    }
    if (nlhs < MIN_OUTPUTS || nlhs > MAX_OUTPUTS) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_filter_int:nlhs",
                "One output required.");
    }
    /* make sure the first input argument is a vector */
    if (mxGetM(prhs[0]) != 1 && mxGetN(prhs[0]) != 1) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_filter_int:notVector",
                "First input must be a vector.");
    }
    /* make sure the first input argument is of type Int16 */
    if (!mxIsInt16(prhs[0]) || mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_filter_int:notInt16",
                "First input must be of type Int16.");
    }
    /* make sure the second input argument is a scalar */
    if (mxGetNumberOfElements(prhs[1]) != 1) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_filter_int:notScalar",
                "Second input must be a scalar.");
    }
    /* make sure the third input argument, if it exists, is a scalar */
    if (mxGetNumberOfElements(prhs[2]) != 1) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_filter_int:notScalar",
                "Second input must be a scalar.");
    }
    
    /* get dimensions of the input vector */
    *nrows = (mwSize)mxGetM(prhs[0]);
    *ncols = (mwSize)mxGetN(prhs[0]);
    
    /* get the second required input  */
    *sampFreq = (int)mxGetScalar(prhs[1]);
    
    /* get the first optional input  */
    if (nrhs > 2) {
        *mFactor = (int)mxGetScalar(prhs[2]);
    }
    else {
        *mFactor = DEFAULT_M;
    }
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    short *inVector;
    int *outVector;
    int sampFreq;
    int mFactor;
    mwSize nrows;
    mwSize ncols;
    mwSize inLen;
    
    /* process arguments */
    processArgs(nlhs, plhs, nrhs, prhs, &nrows, &ncols,
            &sampFreq, &mFactor);
    
    /* calculate the length of the vectors */
    inLen = max(nrows,ncols);
    
    /* get a pointer to the data in the input vector  */
    inVector = (short*)mxGetData(prhs[0]);
    
    /* create the output vector */
    plhs[0] = mxCreateNumericMatrix(nrows,ncols,mxINT32_CLASS,mxREAL);
    
    /* get a pointer to the data in the output vector */
    outVector = (int*)mxGetData(plhs[0]);
    
    /* Initialize global variables */
    initVariables(sampFreq, mFactor);
    
    /* process each sample of the input vector, and obtain an output sample */
    tic();
    for (mwSize i = 0; i < inLen; i++) {
        outVector[i] = processSample(inVector[i]);
    }
    toc();
    
    /* create the output scalar */
    if (nlhs > 1) {
        double delay = lpfhDelay + wmfDelay + agcDelay + mobdDelay + mafhDelay;
        plhs[1] = mxCreateDoubleScalar(delay);
    }
    
    /* Finalize global variables */
    finVariables();
}
