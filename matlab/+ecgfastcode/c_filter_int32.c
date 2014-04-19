/*=========================================================================
 * c_filter_int32.c
 *=======================================================================*/
#include <math.h>
#include "mex.h"

#define FILTER_GAIN_FACTOR  2
#define TARGET_RESOLUTION   7
#define DEFAULT_M_FACTOR    3
#define DEFAULT_RESOLUTION  12
#define INTEGRATION_LEN_MS  100

const int hVector[] = {1,0,0,0,0,0,-2,0,0,0,0,0,1};

/* Multiply and accumulate */
int multacc(const int *x, const int *h, mwSize start, mwSize end)
{
    int acc = 0;
    for (mwSize k = start; k < end; k++) {
        acc += h[k] * x[-k];
    }
    return acc;
}

/* Multiply together */
int multmult(const int *x, mwSize start, mwSize end)
{
    int prod = 1;
    for (mwSize k = start; k < end; k++) {
        prod *= x[-k];
    }
    return prod;
}

/* Normal convolution */
void convolve( const int *x, mwSize nx,
               const int *h, mwSize nh,
               int *y, int divFactor)
{
    mwSize offset = nh/2;
    
    for (mwSize i = offset; i < nh; i++) {
        y[i-offset] = multacc(x+i, h, 0, i+1) >> divFactor;
    }
    for (mwSize i = nh; i < nx; i++) {
        y[i-offset] = multacc(x+i, h, 0, nh) >> divFactor;
    }
    for (mwSize i = nx; i < nx+offset; i++) {
        y[i-offset] = multacc(x+i, h, i-nx+1, nh) >> divFactor;
    }
}

/* Non-linear convolution  */
void nlconvolve(const int *x, mwSize nx, int m, int *y)
{
    mwSize offset = m/2;
    
    for (mwSize i = offset; i < m; i++) {
        y[i-offset] = abs(multmult(x+i, 0, i+1));
    }
    for (mwSize i = m; i < nx; i++) {
        y[i-offset] = abs(multmult(x+i, 0, m));
    }
    for (mwSize i = nx; i < nx+offset; i++) {
        y[i-offset] = abs(multmult(x+i, i-nx+1, m));
    }
}

/* Quantize data to fewer bits */
void quantize(int *x, mwSize nx, int numBits, int targetRes)
{
    if (numBits > targetRes) {
        int factor = numBits-targetRes;
        for (mwSize i = 0; i < nx; i++) {
            x[i] >>= factor;
        }
    }
}

/* Integration filter impulse response */
int *createIntegrationFilter(int Fs, mwSize *hlen, int *factor)
{
    int *h;
    mwSize n;
    
    n = (int)(Fs*INTEGRATION_LEN_MS*0.0005)*2+1;
    h = (int*)mxMalloc(n*sizeof(int));
    for (int i = 0; i < n; i++) {
        h[i] = 1;
    }
    
    *hlen = n;
    *factor = round(log2((double)n));
    return h;
}

/* Process input and output arguments */
void processArgs( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[],
                  mwSize *nrows, mwSize *ncols,
                  int *sampFreq, int *numBits, int *mFactor)
{
    /* check for proper number of arguments */
    if (nrhs < 2 || nrhs > 4) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_filter:nrhs",
                "Two input required and two optional.");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_filter:nlhs",
                "One output required.");
    }
    /* make sure the first input argument is a vector */
    if (mxGetM(prhs[0]) != 1 && mxGetN(prhs[0]) != 1) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_filter:notVector",
                "First input must be a vector.");
    }
    /* make sure the first input argument is of type Int32 */
    if (!mxIsInt32(prhs[0]) || mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_filter:notInt32",
                "First input must be of type Int32.");
    }
    /* make sure the second input argument is a scalar */
    if (mxGetNumberOfElements(prhs[1]) != 1) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_filter:notScalar",
                "Second input must be a scalar.");
    }
    /* make sure the third input argument, if it exists, is a scalar */
    if (nrhs > 2 && mxGetNumberOfElements(prhs[2]) != 1) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_filter:notScalar",
                "Third input must be a scalar.");
    }
    /* make sure the third input argument, if it exists, is a scalar */
    if (nrhs > 3 && mxGetNumberOfElements(prhs[3]) != 1) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_filter:notScalar",
                "Fourth input must be a scalar.");
    }
    
    /* get dimensions of the input vector */
    *nrows = (mwSize)mxGetM(prhs[0]);
    *ncols = (mwSize)mxGetN(prhs[0]);
    
    /* get the second required input  */
    *sampFreq = (int)mxGetScalar(prhs[1]);
    
    /* get the first optional input  */
    if (nrhs > 2) {
        *numBits = (int)mxGetScalar(prhs[2]);
    }
    else {
        *numBits = DEFAULT_RESOLUTION;
    }
    
    /* get the second optional input  */
    if (nrhs > 3) {
        *mFactor = (int)mxGetScalar(prhs[3]);
    }
    else {
        *mFactor = DEFAULT_M_FACTOR;
    }
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    int *inVector;
    int *outVector;
    int sampFreq;
    int numBits;
    int mFactor;
    int *hVector2;
    int *tempVector;
    int filterFactor2;
    
    mwSize nrows;
    mwSize ncols;
    mwSize inLen;
    mwSize hLen;
    mwSize hLen2;
    
    /* process arguments */
    processArgs(nlhs, plhs, nrhs, prhs,
            &nrows, &ncols, &sampFreq, &numBits, &mFactor);
    
    /* calculate the length of the vectors */
    inLen = max(nrows,ncols);
    hLen = sizeof(hVector)/sizeof(int);
    
    /* get a pointer to the data in the input vector  */
    inVector = (int*)mxGetData(prhs[0]);
    
    /* create the output vector */
    plhs[0] = mxCreateNumericMatrix(nrows,ncols,mxINT32_CLASS,mxREAL);
    
    /* get a pointer to the data in the output vector */
    outVector = (int*)mxGetData(plhs[0]);
    
    /* create a temporary vector */
    tempVector = (int*)mxMalloc(inLen*sizeof(int));
    
    /* create the integration filter impulse response */
    hVector2 = createIntegrationFilter(sampFreq,&hLen2,&filterFactor2);
    
    /* call the computational routine 1 */
    convolve(inVector, inLen, hVector, hLen, outVector, FILTER_GAIN_FACTOR);
    
    /* call the computational routine 2 */
    //quantize(tempVector, inLen, numBits, TARGET_RESOLUTION);
    
    /* call the computational routine 3 */
    nlconvolve(outVector, inLen, mFactor, tempVector);
    
    /* call the computational routine 4 */
    convolve(tempVector, inLen, hVector2, hLen2, outVector, filterFactor2);
    
    /* deallocate memory */
    mxFree(tempVector);
    mxFree(hVector2);
}
