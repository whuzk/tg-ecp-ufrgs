/*=========================================================================
 * c_filter_double.c
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

static const double lpfhVector[] = {
    1.0/4,
    0.0/4,
    0.0/4,
    0.0/4,
    0.0/4,
    0.0/4,
   -2.0/4,
    0.0/4,
    0.0/4,
    0.0/4,
    0.0/4,
    0.0/4,
    1.0/4
};

/* Moving-average filter impulse response */
double *createMovingAverageFilter(mwSize nh)
{
    double *h = (double*)mxMalloc(nh*sizeof(double));
    double factor = 1 << (int)ceil(log2(nh));//nh;
    
    for (int i = 0; i < nh; i++) {
        h[i] = 1.0/factor;
    }
    return h;
}

/* Automatic gain control */
void gain_control( const double *x, mwSize nx, double *y,
                   int Fs, double ming, double maxg, mwSize delay)
{
    double *g = (double*)mxMalloc(nx*sizeof(double));
    
    maxfilter(x, nx, 2*Fs, g, -HUGE_VAL);
    
    for (mwSize i = 0; i < delay; i++) {
        y[i] = 0;
    }
    for (mwSize i = delay; i < nx; i++) {
        y[i] = fmin(maxg, x[i-delay] * maxg / fmax(ming, g[i]));
    }
    
    mxFree(g);
}

/* Process the input vector */
double processInput(const double *x, mwSize nx, double *y, int Fs, int m)
{
    double *tempVector;
    double *mafhVector;
    double delay = 0;
    double maxGain;
    double minGain;
    mwSize lpfhLen;
    mwSize mafhLen;
    mwSize gainLen;
    
    /* calculate length of filter impulse responses */
    lpfhLen = sizeof(lpfhVector)/sizeof(double);
    mafhLen = (int)(Fs*MAFH_LEN_MS*0.0005)*2+1;
    gainLen = (int)(Fs*GAIN_LEN_MS*0.001);
    maxGain = (1 << min(15,(8*sizeof(int)-1)/m)) - 1;
    minGain = 1 << 3;
    
    /* create a temporary vector */
    tempVector = (double*)mxMalloc(nx*sizeof(double));
    
    /* create the integration filter impulse response */
    mafhVector = createMovingAverageFilter(mafhLen);
    
    /* call the computational routine 1 */
    fir(x, nx, lpfhVector, lpfhLen, tempVector);
    delay += (lpfhLen-1)/2.0;
    
    /* call the computational routine 2 */
    absolute(tempVector, nx, tempVector);
    gain_control(tempVector, nx, y, Fs, minGain, maxGain, gainLen);
    delay += gainLen;
    
    /* call the computational routine 3 */
    nlfir(y, nx, m, tempVector);
    delay += (m-1)/2.0;
    
    /* call the computational routine 4 */
    fir(tempVector, nx, mafhVector, mafhLen, y);
    delay += (mafhLen-1)/2.0;
    
    /* deallocate memory */
    mxFree(tempVector);
    mxFree(mafhVector);
    
    return delay;
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
                "EcgToolbox:c_filter_double:nrhs",
                "One input required and two optional.");
    }
    if (nlhs < MIN_OUTPUTS || nlhs > MAX_OUTPUTS) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_filter_double:nlhs",
                "One output required.");
    }
    /* make sure the first input argument is a vector */
    if (mxGetM(prhs[0]) != 1 && mxGetN(prhs[0]) != 1) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_filter_double:notVector",
                "First input must be a vector.");
    }
    /* make sure the first input argument is of type Double */
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_filter_double:notDouble",
                "First input must be of type double.");
    }
    /* make sure the second input argument is a scalar */
    if (mxGetNumberOfElements(prhs[1]) != 1) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_filter_double:notScalar",
                "Second input must be a scalar.");
    }
    /* make sure the third input argument, if it exists, is a scalar */
    if (nrhs > 2 && mxGetNumberOfElements(prhs[2]) != 1) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_filter_double:notScalar",
                "Third input must be a scalar.");
    }
    
    /* get dimensions of the input vector */
    *nrows = (mwSize)mxGetM(prhs[0]);
    *ncols = (mwSize)mxGetN(prhs[0]);
    
    /* get the second required input  */
    *sampFreq = (int)mxGetScalar(prhs[1]);
    
    /* make sure the sampling frequency is positive */
    if (*sampFreq <= 0) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_filter_double:freqNotPositive",
                "Sampling frequency must be positive.");
    }
    
    /* get the first optional input  */
    if (nrhs > 2) {
        *mFactor = (int)mxGetScalar(prhs[2]);
    }
    else {
        *mFactor = DEFAULT_M;
    }
    
    /* make sure the filter order is within pre-defined limits */
    if (*mFactor < 1 || *mFactor > 8) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_filter_double:badFilterOrder",
                "Filter order must be between 1 and 8.");
    }
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double *inVector;
    double *outVector;
    double delay;
    int sampFreq;
    int mFactor;
    mwSize nrows;
    mwSize ncols;
    mwSize inLen;
    
    /* process arguments */
    processArgs(nlhs, plhs, nrhs, prhs, &nrows, &ncols,
            &sampFreq, &mFactor);
    
    /* calculate the size of vectors */
    inLen = max(nrows,ncols);
    
    /* get a pointer to the data in the input vector  */
    inVector = mxGetPr(prhs[0]);
    
    /* create the output vector */
    plhs[0] = mxCreateDoubleMatrix(nrows,ncols,mxREAL);
    
    /* get a pointer to the data in the output vector */
    outVector = mxGetPr(plhs[0]);
    
    /* process the input vector, and obtain the delay */
    tic();
    delay = processInput(inVector, inLen, outVector, sampFreq, mFactor);
    toc();
    
    /* create the output scalar */
    if (nlhs > 1) {
        plhs[1] = mxCreateDoubleScalar(delay);
    }
}
