/*=========================================================================
 * c_detect_qrs_double.c
 *=======================================================================*/
#include <math.h>
#include "mex.h"
#include "c_ecg_utils.h"

#define MIN_INPUTS      2
#define MAX_INPUTS      4
#define MIN_OUTPUTS     1
#define MAX_OUTPUTS     2

#define MIN_MAINS_FREQ  50
#define MAX_MAINS_FREQ  60
#define DEF_MAINS_FREQ  60

#define MIN_MBD_ORDER   1
#define MAX_MBD_ORDER   8
#define DEF_MBD_ORDER   3

#define ACG_LEN_SEC     0.05
#define MAF_LEN_SEC     0.10

#define LPF_BUFLEN      256
#define WMF_BUFLEN      1024
#define ACG_BUFLEN      1024    // 2^nextpow2(MAX_MAINS_FREQ * LPF_BUFLEN * ACG_LEN_SEC)
#define MBD_BUFLEN      8       // 2^nextpow2(MAX_MBD_ORDER)
#define MAF_BUFLEN      2048    // 2^nextpow2(MAX_MAINS_FREQ * LPF_BUFLEN * MAF_LEN_SEC)

static double lpfBuf[LPF_BUFLEN] = {0};
static double wmfBuf[WMF_BUFLEN];
static unsigned int wmfAux[WMF_BUFLEN];
static double agcBuf[ACG_BUFLEN] = {0};
static double mbdBuf[MBD_BUFLEN] = {1};
static double mafBuf[MAF_BUFLEN] = {0};

#define lpfb(I) (lpfBuf[(I)&(LPF_BUFLEN-1)])
#define wmfb(I) (wmfBuf[(I)&(WMF_BUFLEN-1)])
#define wmfa(I) (wmfAux[(I)&(WMF_BUFLEN-1)])
#define agcb(I) (agcBuf[(I)&(ACG_BUFLEN-1)])
#define mbdb(I) (mbdBuf[(I)&(MBD_BUFLEN-1)])
#define mafb(I) (mafBuf[(I)&(MAF_BUFLEN-1)])

static int sampFreq;
static int mainsFreq;
static int mbdOrder;

static double lpfDel;
static mwSize lpfWin;
static double wmfDel;
static mwSize wmfWin;
static double agcDel;
static mwSize agcWin;
static short agcMin;
static short agcMax;
static double mbdDel;
static mwSize mbdWin;
static double mafDel;
static mwSize mafWin;


/* Low-pass filter and second-order backward difference */
double lpf(double sample)
{
    static unsigned short n = 0;
    double y0;
    
    // update filter output: y[n] = x[n] - 2*x[n-M/2] + x[n-M+1]
    y0 = sample - 2 * lpfb(n - lpfWin/2) + lpfb(n - lpfWin + 1);
    
    // update filter memory
    lpfb(n++) = sample;
    
    // compute result
    return y0 / 4.0;
}

/* Windowed-max filter */
double wmf(double sample)
{
    static unsigned short first = 0;
    static unsigned short count = 0;
    static unsigned int i = 0;
    unsigned short j,k;
    
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
    if (count > wmfWin || wmfa(first) == i - wmfWin) {
        first++;
        count--;
    }
    // increment global index
    i++;
    // return the max in the window
    return wmfb(first);
}

/* Automatic gain control */
double agc(double sample, double gain)
{
    static unsigned short n = 0;
    double y0;
    
    // update filter output: y[n] = x[n-M] * maxG / max(G, minG)
    y0 = (agcb(n - agcWin) * agcMax) / fmax(gain, agcMin);
    
    // update filter memory
    agcb(n++) = sample;
    
    // compute result
    return fmin(y0, agcMax);
}

/* Multiplication of absolute backward differences */
double mdb(double sample)
{
    static unsigned short n = 0;
    double y0 = sample;
    
    // update filter output: y[n] = x[n] * ... * x[n-M+1]
    for (mwSize k = 1; k < mbdWin; k++) {
        y0 *= mbdb(n-k);
    }
    
    // update filter memory
    mbdb(n++) = sample;
    
    // compute result
    return y0;
}

/* Moving-average filter of length 32 */
double maf(double sample)
{
    static unsigned short n = 0;
    static double y0 = 0;
    
    // update filter output: y[n] = y[n-1] + x[n] - x[n-M]
    y0 += sample - mafb(n - mafWin);
    
    // update filter memory
    mafb(n++) = sample;
    
    // compute result
    return fmin(y0 / mafWin, INT_MAX);
}

/* Process one intput sample */
double processSample(double sample)
{
    double gain;
    sample = lpf(sample);
    sample = fabs(sample);
    gain = wmf(sample);
    sample = agc(sample, gain);
    sample = mdb(sample);
    sample = maf(sample);
    return sample;
}

/* Initialize global variables */
void initVariables()
{
    short maxbits;
    
    // low-pass filter and second-order backward difference
    lpfWin = ((int)round(sampFreq / mainsFreq) << 1) + 1;
    lpfDel = (lpfWin - 1) / (double)2.0;
    //printf("lpd delay = %f\n", lpfDel);
    
    // windowed-max filter
    wmfWin = 2 * sampFreq;
    wmfDel = 0;
    //printf("wmf delay = %f\n", wmfDel);
    
    // automatic gain control
    agcWin = (int)(sampFreq * (double)ACG_LEN_SEC);
    agcDel = (double)agcWin;
    maxbits = (short)((8 * sizeof(int)) - 1) / mbdOrder;
    agcMax = (1 << (short)min(15,maxbits)) - 1;
    agcMin = 1 << 3;
    //printf("agc delay = %f\n", agcDel);
    
    // multiplication of absolute backward differences
    mbdWin = mbdOrder;
    mbdDel = (mbdWin - 1) / (double)2.0;
    //printf("mbd delay = %f\n", mbdDel);
    
    // moving-average filter
    mafWin = (int)(sampFreq * (double)MAF_LEN_SEC);
    mafDel = (mafWin - 1) / (double)2.0;
    //printf("maf delay = %f\n", mafDel);
}

/* Process input and output arguments */
void processArgs( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[],
                  mwSize *nrows, mwSize *ncols)
{
    /* check for proper number of arguments */
    if (nrhs < MIN_INPUTS || nrhs > MAX_INPUTS) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_detect_qrs_double:nrhs",
                "Two inputs required.");
    }
    if (nlhs < MIN_OUTPUTS || nlhs > MAX_OUTPUTS) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_detect_qrs_double:nlhs",
                "One output required.");
    }
    /* make sure the first input argument is a Vector */
    if (mxGetM(prhs[0]) != 1 && mxGetN(prhs[0]) != 1) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_detect_qrs_double:notVector",
                "First input must be a Vector.");
    }
    /* make sure the first input argument is of type Double */
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_detect_qrs_double:notDouble",
                "First input must be of type Double.");
    }
    /* make sure the remaining arguments are all scalars */
    for (mwSize i = 1; i < nrhs; i++) {
        if (mxGetNumberOfElements(prhs[i]) != 1) {
            mexErrMsgIdAndTxt(
                    "EcgToolbox:c_detect_qrs_double:notScalar",
                    "Input #%d must be a scalar.", i+1);
        }
    }
    
    /* get dimensions of the input vector */
    *nrows = (mwSize)mxGetM(prhs[0]);
    *ncols = (mwSize)mxGetN(prhs[0]);
    
    /* get the sampling frequency  */
    sampFreq = (int)mxGetScalar(prhs[1]);
    
    /* get the mains frequency  */
    if (nrhs > 2) {
        mainsFreq = (int)mxGetScalar(prhs[2]);
    }
    else {
        mainsFreq = DEF_MAINS_FREQ;
    }
    
    /* get the mbd filter order  */
    if (nrhs > 3) {
        mbdOrder = (int)mxGetScalar(prhs[3]);
    }
    else {
        mbdOrder = DEF_MBD_ORDER;
    }
    
    /* make sure the mains frequency is within pre-defined limits */
    if (mainsFreq < MIN_MAINS_FREQ || mainsFreq > MAX_MAINS_FREQ) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_detect_qrs_double:badMainsFreq",
                "Mains frequency must be between %d and %d.", MIN_MAINS_FREQ, MAX_MAINS_FREQ);
    }
    /* make sure the sampling frequency is within pre-defined limits */
    if (sampFreq < mainsFreq || sampFreq > mainsFreq * LPF_BUFLEN) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_detect_qrs_double:badSampFreq",
                "Sampling frequency must be between %d and %d.", mainsFreq, mainsFreq * LPF_BUFLEN);
    }
    /* make sure the MBD filter order is within pre-defined limits */
    if (mbdOrder < MIN_MBD_ORDER || mbdOrder > MAX_MBD_ORDER) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_detect_qrs_double:badMbdOrder",
                "MBD filter order must be between %d and %d.", MIN_MBD_ORDER, MAX_MBD_ORDER);
    }
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double *inVector;
    double *outVector;
    mwSize nrows;
    mwSize ncols;
    mwSize inLen;
    
    /* process arguments */
    processArgs(nlhs, plhs, nrhs, prhs, &nrows, &ncols);
    
    /* calculate the length of the Vectors */
    inLen = max(nrows,ncols);
    
    /* get a pointer to the data in the input Vector  */
    inVector = mxGetPr(prhs[0]);
    
    /* create the output Vector */
    plhs[0] = mxCreateDoubleMatrix(nrows,ncols,mxREAL);
    
    /* get a pointer to the data in the output Vector */
    outVector = mxGetPr(plhs[0]);
    
    /* Initialize global variables */
    initVariables();
    
    /* process each sample of the input Vector, and obtain an output sample */
    tic();
    for (mwSize i = 0; i < inLen; i++) {
        outVector[i] = processSample(inVector[i]);
    }
    toc();
    
    /* create the output delay */
    if (nlhs > 1) {
        double delay = lpfDel + wmfDel + agcDel + mbdDel + mafDel;
        plhs[1] = mxCreateDoubleScalar(delay);
    }
}
