/*=========================================================================
 * c_detect_qrs_int.c
 *=======================================================================*/
#include <math.h>
#include "mex.h"
#include "c_ecg_utils.h"

#define MIN_INPUTS  2
#define MAX_INPUTS  3
#define MIN_OUTPUTS 1
#define MAX_OUTPUTS 2

#define DEFAULT_M   3
#define ACG_LEN_MS  50
#define MAF_LEN_MS  100

// low-pass filter and second-order backward difference
static short *lpfBuf;
static mwSize lpfLen;
static mwSize lpfDel;
#define lpfb(I) (lpfBuf[(I)&(lpfLen-1)])

// windowed-max filter
static short *wmfBuf;
static unsigned int *wmfAux;
static mwSize wmfLen;
static mwSize wmfDel;
#define wmfb(I) (wmfBuf[(I)&(wmfLen-1)])
#define wmfa(I) (wmfAux[(I)&(wmfLen-1)])

// automatic gain control
static short *agcBuf;
static mwSize agcLen;
static mwSize agcDel;
static short agcMin;
static short agcMax;
#define agcb(I) (agcBuf[(I)&(agcLen-1)])

// multiplication of absolute backward differences
static int *mbdBuf;
static mwSize mbdLen;
static mwSize mbdDel;
static short mdbFac;
#define mbdb(I) (mbdBuf[(I)&(mbdLen-1)])

// moving-average filter
static int *mafBuf;
static mwSize mafLen;
static mwSize mafDel;
static mwSize mafLog2Gain;
#define mafb(I) (mafBuf[(I)&(mafLen-1)])

/* Low-pass filter and second-order backward difference */
short lpf(short sample)
{
    static mwSize n = 0;
    int y0;
    
    // update filter output: y[n] = x[n] - 2*x[n-delay] + x[n-M+1]
    y0 = sample - (lpfb(n - lpfDel) << 1) + lpfb(n - (lpfDel << 1));
    
    // update filter memory
    lpfb(n++) = sample;
    
    // compute result
    return y0 >> 2;
}

/* Windowed-max filter */
short wmf(short sample)
{
    static mwSize first = 0;
    static mwSize count = 0;
    static unsigned int i = 0;
    mwSize j,k;
    
    // search for the first element greater than the current sample
    j = count;
    k = first + count - 1;
    while (j > 0 && sample >= wmfb(k)) {
        j--;
        k--;
    }
    // put the sample next to the element found and adjust the length
    wmfb(k + 1) = sample;
    wmfa(k + 1) = i;
    count = j + 1;
    // check if the first in line has gone out of the window length
    if (count > wmfLen || wmfa(first) == i - wmfLen) {
        first++;
        count--;
    }
    // increment global index
    i++;
    // return the max in the window
    return wmfb(first);
}

/* Automatic gain control */
short agc(short sample, short gain)
{
    static mwSize n = 0;
    int y0;
    
    // update filter output: y[n] = x[n-delay] * maxG / max(minG, G)
    y0 = (agcb(n - agcDel) * agcMax) / max(gain, agcMin);
    
    // update filter memory
    agcb(n++) = sample;
    
    // compute result
    return (short)min(y0, agcMax);
}

/* Multiplication of absolute backward differences */
int modb(short sample)
{
    static mwSize n = 0;
    int y0 = sample;
    
    // update filter output: y[n] = x[n] * ... * x[n-M+1]
    for (mwSize k = 0; k < mdbFac-1; k++)
        y0 *= mbdb(n-k);
    
    // update filter memory
    mbdb(n++) = sample;
    
    // compute result
    return y0;
}

/* Moving-average filter of length 32 */
int maf(int sample)
{
    static mwSize n = 0;
    static long long y0 = 0;
    
    // update filter output: y[n] = y[n-1] + x[n] - x[n-M+1]
    y0 += sample - mafBuf[n - (mafDel << 1)];
    
    // update filter memory
    mafb(n++) = sample;
    
    // compute result
    return (int)min(y0 >> mafLog2Gain, INT_MAX);
}

/* Process one intput sample */
int processSample(short sample)
{
    sample = abs(lpf(sample));
    sample = agc(sample, wmf(sample));
    return maf(modb(sample));
}

/* Initialize global variables */
void initVariables(int Fs, short m)
{
    short log2max;
    
    // low-pass filter and second-order backward difference
    lpfDel = 6;
    lpfLen = (int)ceil(log2((lpfDel << 1) + 1));
    lpfBuf = (short*)mxCalloc(lpfLen, sizeof(short));
    
    // windowed-max filter
    wmfDel = 0;
    wmfLen = (int)ceil(log2(Fs));
    wmfBuf = (short*)mxMalloc(wmfLen * sizeof(short));
    wmfAux = (unsigned int*)mxMalloc(wmfLen * sizeof(int));

    // automatic gain control
    agcDel = (int)(Fs * ACG_LEN_MS * 0.0005);
    agcLen = (int)ceil(log2((agcDel << 1) + 1));
    agcBuf = (short*)mxCalloc(agcLen, sizeof(short));
    log2max = (short)min(15,((sizeof(int) << 3) - 1) / m);
    agcMax = (1 << log2max) - 1;
    agcMin = 1 << 3;
    
    // multiplication of absolute backward differences
    mdbFac = m;
    mbdDel = (mdbFac - 1) >> 1;
    mbdLen = (int)ceil(log2(mdbFac));
    mbdBuf = (int*)mxMalloc(mbdLen * sizeof(int));
    for (mwSize k = 0; k < mbdLen; k++) {
        mbdBuf[k] = 1;
    }
    
    // moving-average filter
    mafDel = (int)(Fs * MAF_LEN_MS * 0.0005);
    mafLen = (int)ceil(log2((mafDel << 1) + 1));
    mafBuf = (int*)mxCalloc(mafLen, sizeof(int));
    mafLog2Gain = (int)ceil(log2(mafLen));
}

/* Finalize global variables */
void finVariables()
{
    mxFree(lpfBuf);
    mxFree(wmfBuf);
    mxFree(wmfAux);
    mxFree(agcBuf);
    mxFree(mbdBuf);
    mxFree(mafBuf);
}

/* Process input and output arguments */
void processArgs( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[],
                  mwSize *nrows, mwSize *ncols,
                  int *sampFreq, short *mFactor)
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
    /* make sure the first input argument is a Buftor */
    if (mxGetM(prhs[0]) != 1 && mxGetN(prhs[0]) != 1) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_filter_int:notBuftor",
                "First input must be a Buftor.");
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
    if (nrhs > 2 && mxGetNumberOfElements(prhs[2]) != 1) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_filter_int:notScalar",
                "Second input must be a scalar.");
    }
    
    /* get dimensions of the input Buftor */
    *nrows = (mwSize)mxGetM(prhs[0]);
    *ncols = (mwSize)mxGetN(prhs[0]);
    
    /* get the second required input  */
    *sampFreq = (int)mxGetScalar(prhs[1]);
    
    /* make sure the sampling frequency is positive */
    if (*sampFreq <= 0) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_filter_int:freqNotPositive",
                "Sampling frequency must be positive.");
    }
    
    /* get the first optional input  */
    if (nrhs > 2) {
        *mFactor = (short)mxGetScalar(prhs[2]);
    }
    else {
        *mFactor = DEFAULT_M;
    }
    
    /* make sure the filter order is within pre-defined limits */
    if (*mFactor < 1 || *mFactor > 8) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_filter_int:badFilterOrder",
                "Filter order must be between 1 and 8.");
    }
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    short *inBuftor;
    int *outBuftor;
    int sampFreq;
    short mFactor;
    mwSize nrows;
    mwSize ncols;
    mwSize inLen;
    
    /* process arguments */
    processArgs(nlhs, plhs, nrhs, prhs, &nrows, &ncols,
            &sampFreq, &mFactor);
    
    /* calculate the length of the Buftors */
    inLen = max(nrows,ncols);
    
    /* get a pointer to the data in the input Buftor  */
    inBuftor = (short*)mxGetData(prhs[0]);
    
    /* create the output Buftor */
    plhs[0] = mxCreateNumericMatrix(nrows,ncols,mxINT32_CLASS,mxREAL);
    
    /* get a pointer to the data in the output Buftor */
    outBuftor = (int*)mxGetData(plhs[0]);
    
    /* Initialize global variables */
    initVariables(sampFreq, mFactor);
    
    /* process each sample of the input Buftor, and obtain an output sample */
    tic();
    for (mwSize i = 0; i < inLen; i++) {
        outBuftor[i] = processSample(inBuftor[i]);
    }
    toc();
    
    /* create the output scalar */
    if (nlhs > 1) {
        int delay = lpfDel + wmfDel + agcDel + mbdDel + mafDel;
        plhs[1] = mxCreateDoubleScalar(delay);
    }
    
    /* Finalize global variables */
    finVariables();
}
