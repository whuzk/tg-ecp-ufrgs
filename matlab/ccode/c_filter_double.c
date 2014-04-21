/*=========================================================================
 * c_filter_double.c
 *=======================================================================*/
#include <math.h>
#include "mex.h"
#include "c_ecg_utils.h"

#define DEFAULT_M_FACTOR    3
#define INTEGRATION_LEN_MS  100

static const double hVector[] = {
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

/* Integration filter impulse response */
double *createIntegrationFilter(int Fs, mwSize *hlen)
{
    mwSize n = 32;//(int)(Fs*INTEGRATION_LEN_MS*0.0005)*2+1;
    double *h = (double*)mxMalloc(n*sizeof(double));
    
    for (int i = 0; i < n; i++) {
        h[i] = 1.0/n;
    }
    
    *hlen = n;
    return h;
}

/* Process input and output arguments */
void processArgs( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[],
                  mwSize *nrows, mwSize *ncols,
                  int *sampFreq, int *mFactor)
{
    /* check for proper number of arguments */
    if (nrhs < 2 || nrhs > 3) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_filter_double:nrhs",
                "One input required and two optional.");
    }
    if (nlhs != 1) {
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
    
    /* get the first optional input  */
    if (nrhs > 2) {
        *mFactor = (int)mxGetScalar(prhs[2]);
    }
    else {
        *mFactor = DEFAULT_M_FACTOR;
    }
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double *inVector;
    double *outVector;
    double *hVector2;
    double *tempVector;
    int sampFreq;
    int mFactor;
    
    mwSize nrows;
    mwSize ncols;
    mwSize inLen;
    mwSize hLen;
    mwSize hLen2;
    
    /* process arguments */
    processArgs(nlhs, plhs, nrhs, prhs,
            &nrows, &ncols, &sampFreq, &mFactor);
    
    /* calculate the length of the vectors */
    inLen = max(nrows,ncols);
    hLen = sizeof(hVector)/sizeof(double);
    
    /* get a pointer to the data in the input vector  */
    inVector = mxGetPr(prhs[0]);
    
    /* create the output vector */
    plhs[0] = mxCreateDoubleMatrix(nrows,ncols,mxREAL);
    
    /* get a pointer to the data in the output vector */
    outVector = mxGetPr(plhs[0]);
    
    /* create a temporary vector */
    tempVector = (double*)mxMalloc(inLen*sizeof(double));
    
    /* create the integration filter impulse response */
    hVector2 = createIntegrationFilter(sampFreq,&hLen2);
    
    /* call the computational routine 1 */
    //convolve(inVector, inLen, hVector, hLen, outVector);
    fir(inVector, inLen, hVector, hLen, outVector);
    
    /* call the computational routine 2 */
    absolute(outVector, inLen);
    
    /* call the computational routine 3 */
    //nlconvolve(tempVector, inLen, mFactor, outVector);
    nlfir(outVector, inLen, mFactor, tempVector);
    
    /* call the computational routine 4 */
    //convolve(tempVector, inLen, hVector2, hLen2, outVector);
    fir(tempVector, inLen, hVector2, hLen2, outVector);
    
    /* deallocate memory */
    mxFree(tempVector);
    mxFree(hVector2);
}
