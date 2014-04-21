/*=========================================================================
 * c_filter_int32.c
 *=======================================================================*/
#include <math.h>
#include "mex.h"
#include "c_ecg_utils.h"

#define SAMPL_FREQ  250
#define LPF_LENGTH  13      // must be odd-valued
#define LPF_FACTOR  2       // factor for division of lpf output
#define MVI_LENGTH  32      // must be a power of two
#define MVI_FACTOR  5       // factor for division of mvi output
#define MOBD_ORDER  3       // MOBD_ORDER and TARGET_RES multiplied
#define TARGET_RES  8       // together must be <= sizeof(int)

/* Low-pass filter and second-order backward difference */
short lpfandsbd(short sample)
{
    static short x[LPF_LENGTH*2] = {0};
    static short n = LPF_LENGTH-1;
    int y;
    
    x[n] = x[n+LPF_LENGTH] = sample;
    y = x[n] - (x[n+LPF_LENGTH/2] << 1) + x[n+LPF_LENGTH-1];
    if (--n < 0) {
        n = LPF_LENGTH-1;
    }
    return y >> LPF_FACTOR;
}

/* Multiplication of absolute backward differences */
int modb(short sample)
{
    static int x[MOBD_ORDER] = {1};
    static int n = 0;
    int prod = 1;
    
    x[n++] = abs(sample);
    if (n == MOBD_ORDER) {
        n = 0;
    }
    for (mwSize k = 0; k < MOBD_ORDER; k++) {
        prod *= x[k];
    }
    return prod;
}

/* Moving-average filter of length 32 */
int mvi(int sample)
{
    static int x[MVI_LENGTH*2] = {0};
    static int n = MVI_LENGTH-1;
    static int y0 = 0, y1;
    
    x[n] = x[n+MVI_LENGTH] = sample;
    y1 = y0;
    y0 = x[n] - x[n+MVI_LENGTH-1] + y1;
    if (--n < 0) {
        n = MVI_LENGTH-1;
    }
    return y0 >> MVI_FACTOR;
}

/* Process one sample input */
int processSample(short sample)
{
    sample = lpfandsbd(sample);
    return mvi(modb(sample));
}

/* Process input and output arguments */
void processArgs( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[],
                  mwSize *nrows, mwSize *ncols,
                  int *sampFreq)
{
    /* check for proper number of arguments */
    if (nrhs != 2) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_filter_int:nrhs",
                "Two inputs required.");
    }
    if (nlhs != 1) {
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
    
    /* get dimensions of the input vector */
    *nrows = (mwSize)mxGetM(prhs[0]);
    *ncols = (mwSize)mxGetN(prhs[0]);
    
    /* get the second required input  */
    *sampFreq = (int)mxGetScalar(prhs[1]);
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    short *inVector;
    int *outVector;
    int sampFreq;
    
    mwSize nrows;
    mwSize ncols;
    mwSize inLen;
    
    /* process arguments */
    processArgs(nlhs, plhs, nrhs, prhs, &nrows, &ncols, &sampFreq);
    
    /* calculate the length of the vectors */
    inLen = max(nrows,ncols);
    
    /* get a pointer to the data in the input vector  */
    inVector = (short*)mxGetData(prhs[0]);
    
    /* create the output vector */
    plhs[0] = mxCreateNumericMatrix(nrows,ncols,mxINT32_CLASS,mxREAL);
    
    /* get a pointer to the data in the output vector */
    outVector = (int*)mxGetData(plhs[0]);
    
    /* main algorithm */
    for (mwSize i = 0; i < inLen; i++) {
        outVector[i] = processSample(inVector[i]);
    }
}
