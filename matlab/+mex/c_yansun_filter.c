/*=========================================================================
 * c_yansun_filter.c
 * 
 *  Title: process an input signal using the Yan Sun's filters
 *  Author:     Diego Sogari
 *  Modified:   07/May/2014
 *
 *  Inputs:
 *      1. input signal
 *      2. sampling frequency (in Hz)
 *      3. delay required for outputs (in samples)
 *  Outputs:
 *      1. filtered signal
 *      2. second filtered signal
 *
 *  *all inputs and outputs required
 *
 *=======================================================================*/
#include <math.h>
#include "mex.h"
#include "c_intfilter.h"
#include "c_maxfilter.h"
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
 *  FILT    filtered signal
 *  LPF     low-pass filter
 *  MAF     moving-average filter
 *  DER     derivative filter
 *  MAXF    maximum filter
 *  MINF    minimum filter
 *  MDF     morphological derivative filter
 *  PDF     pure delay filter
 *=======================================================================*/

/*=========================================================================
 * Constants
 *=======================================================================*/
#define MIN_INPUTS      2       // minimum number of input arguments
#define MAX_INPUTS      3       // maximum number of input arguments
#define MIN_OUTPUTS     2       // minimum number of output arguments
#define MAX_OUTPUTS     2       // maximum number of output arguments

#define MIN_SAMP_FREQ   50*2    // minimum sampling frequency
#define MAX_SAMP_FREQ   60*40   // maximum sampling frequency

#define LPF_ORDER       2       // low-pass filter order
#define LPF_CUTOFF      11      // low-pass filter cutoff frequency (in Hz)
#define MAF_WIDTH       0.02    // moving-average filter width (in seconds)
#define DER_ORDERN      1       // derivative filter order N
#define DER_ORDERM      0       // derivative filter order M
#define MDF_WIDTH       0.06    // MD filter width (in seconds)

/*=========================================================================
 * Input and output variables
 *=======================================================================*/
static double *inputSig;        // input signal vector
static mwSize inputLen;         // input signal length
static double sampFreq;         // sampling frequency
static double delay;            // delay required for outputs
static double *filtSig1;        // first filtered signal
static double *filtSig2;        // second filtered signal

/*=========================================================================
 * Filter variables
 *=======================================================================*/
static unsigned int ci;         // current sample index
static intfobject lpFilter;     // low-pass filter object
static intfobject maFilter;     // moving-average filter object
static intfobject deFilter;     // derivative filter object
static intfobject pdFilter1;    // pure delay filter object 1
static maxfobject maxFilter;    // maximum filter object
static maxfobject minFilter;    // minimum filter object
static intfobject mdFilter;     // morphological derivative filter object
static intfobject pdFilter2;    // pure delay filter object 2
static int lpGainLog2;          // log2 of low-pass filter gain
static int maGainLog2;          // log2 of moving-average filter gain
static int deGainLog2;          // log2 of derivative filter gain

/*=========================================================================
 * Low-pass filter
 *=======================================================================*/
int lpf(int sample)
{
    return intfnewx(&lpFilter, ci, sample);// >> lpGainLog2;
}

/*=========================================================================
 * Moving-average filter
 *=======================================================================*/
int maf(int sample)
{
    return intfnewx(&maFilter, ci, sample);// >> maGainLog2;
}

/*=========================================================================
 * Derivative filter
 *=======================================================================*/
int def(int sample)
{
    return intfnewx(&deFilter, ci, sample);// >> deGainLog2;
}

/*=========================================================================
 * Pure delay filter 1
 *=======================================================================*/
int pdf1(int sample)
{
    return intfnewx(&pdFilter1, ci, sample);
}

/*=========================================================================
 * Maximum filter
 *=======================================================================*/
int maxf(int sample)
{
    return maxfnewx(&maxFilter, ci, sample);
}

/*=========================================================================
 * Minimum filter
 *=======================================================================*/
int minf(int sample)
{
    return maxfnewx(&minFilter, ci, sample);
}

/*=========================================================================
 * Morphological derivative filter
 *=======================================================================*/
int mdf(int sample)
{
    return maxf(sample) + minf(sample) -
            (intfnewx(&mdFilter, ci, sample) << 1);
}

/*=========================================================================
 * Pure delay filter 2
 *=======================================================================*/
int pdf2(int sample)
{
    return intfnewx(&pdFilter2, ci, sample);
}

/*=========================================================================
 * Process one input sample, and return an output sample
 *=======================================================================*/
int processSample(int sample, int *mdsamp)
{
    sample = lpf(sample);       // low-pass
    *mdsamp = mdf(sample);      // morphological derivative
    *mdsamp = pdf2(*mdsamp);    // delay of morphological derivative
    
    sample = maf(sample);       // moving-average
    sample = def(sample);       // derivative
    return pdf1(sample);        // delay of derivative
}

/*=========================================================================
 * Event triggered for an arriving sample
 *=======================================================================*/
void onNewSample(double sample)
{
    int mdsamp;
    
    // process sample
    sample = (double)processSample((int)sample, &mdsamp);
    
    // save result
    filtSig1[ci] = sample;
    filtSig2[ci] = (double)mdsamp;
    
    // increment global lindex
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
            "EcgToolbox:c_yansun_filter:nrhs",
            "%d input(s) required and %d optional",
            MIN_INPUTS, MAX_INPUTS - MIN_INPUTS);
    }
    // check for proper number of output arguments
    if (nlhs < MIN_OUTPUTS || nlhs > MAX_OUTPUTS) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_yansun_filter:nlhs",
            "%d output(s) required and %d optional",
            MIN_OUTPUTS, MAX_OUTPUTS - MIN_OUTPUTS);
    }
    // make sure the first input argument is a vector
    if (mxGetM(prhs[0]) != 1 && mxGetN(prhs[0]) != 1) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_yansun_filter:notVector",
            "First input must be a vector.");
    }
    // make sure the first input argument is of type double
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_yansun_filter:notDouble",
            "First input must be of type double.");
    }
    // make sure the remaining arguments are all scalars
    for (mwSize i = 1; i < nrhs; i++) {
        if (mxGetNumberOfElements(prhs[i]) != 1) {
            mexErrMsgIdAndTxt(
                "EcgToolbox:c_yansun_filter:notScalar",
                "Input #%d must be a scalar.", i + 1);
        }
    }
}

/*=========================================================================
 * Hanlde input arguments 
 *=======================================================================*/
void handleInputs( int nrhs, const mxArray *prhs[],
                   mwSize *nrows, mwSize *ncols)
{
    // get a pointer to the data in the input vector
    inputSig = mxGetPr(prhs[0]);
    
    // get the dimensions of the input vector
    *nrows = (mwSize)mxGetM(prhs[0]);
    *ncols = (mwSize)mxGetN(prhs[0]);
    
    // get the sampling frequency 
    sampFreq = mxGetScalar(prhs[1]);
    
    // make sure the sampling frequency is within pre-defined limits
    if (sampFreq < MIN_SAMP_FREQ || sampFreq > MAX_SAMP_FREQ) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_yansun_filter:badSampFreq",
            "Sampling frequency must be between %d and %d.",
            MIN_SAMP_FREQ, MAX_SAMP_FREQ);
    }
    
    // get the required delay
    delay = mxGetScalar(prhs[2]);
}

/*=========================================================================
 * Hanlde output arguments 
 *=======================================================================*/
void handleOutputs( int nlhs, mxArray *plhs[],
                    mwSize nrows, mwSize ncols)
{
    // get a pointer to the first output vector
    plhs[0] = mxCreateDoubleMatrix(nrows, ncols, mxREAL);
    filtSig1 = mxGetPr(plhs[0]);
    
    // get a pointer to the second output vector
    plhs[1] = mxCreateDoubleMatrix(nrows, ncols, mxREAL);
    filtSig2 = mxGetPr(plhs[1]);
}

/*=========================================================================
 * The initialization routine 
 *=======================================================================*/
void init()
{
    int width = (int)round(MDF_WIDTH * sampFreq);
    
    // initialize sample index
    ci = 0;
    
    // initialize filter objects
    initintfilter(lpFilter);
    initintfilter(maFilter);
    initintfilter(deFilter);
    initintfilter(pdFilter1);
    initmaxfilter(maxFilter);
    initmaxfilter(minFilter);
    initintfilter(mdFilter);
    initintfilter(pdFilter2);
    
    // design low-pass filter
    design_lowpass(&lpFilter, sampFreq, "N,Fc,3db", LPF_ORDER,
            (double)LPF_CUTOFF);
    
    // make sure the required delay is sufficiently large
    if (delay < lpFilter.delay + width) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_yansun_filter:badDelay",
            "Required delay must be at least equal to %.2f.",
            lpFilter.delay + width);
        void finalize();
        finalize();
    }
    
    // design moving-average filter
    design_maverage(&maFilter, sampFreq, MAF_WIDTH);
    
    // design derivative filter
    design_derivative(&deFilter, DER_ORDERN, DER_ORDERM);
    
    // create an all-pass filter for delay of the derivative signal
    design_allpass(&pdFilter1, 1, (int)ceil(delay -
            (lpFilter.delay + maFilter.delay + deFilter.delay)));
    
    // create maximum filter
    create_maxfilter(&maxFilter, width * 2 + 1);
    
    // create minimum filter
    create_minfilter(&minFilter, width * 2 + 1);
    
    // create a pure delay filter for the morphological derivative
    design_allpass(&mdFilter, 1, width);
    
    // create an all-pass filter for delay of the morphological derivative
    design_allpass(&pdFilter2, 1, (int)ceil(delay -
            (lpFilter.delay + width)));
    
    // calculate log2 of filter gains (to speedup division)
    lpGainLog2 = 1 + ilogb(lpFilter.gain - 1);
    maGainLog2 = 1 + ilogb(maFilter.gain - 1);
    deGainLog2 = 1 + ilogb(deFilter.gain - 1);
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
        onNewSample(inputSig[i]);
    }
    
    // stop time counter
    time = toc();
    
    // display time statistics
    mexPrintf("Total processing time: %.2f ms\n", 1000 * time);
    mexPrintf("Average time per sample: %.2f ns\n",
            1000000000 * time / inputLen);
}

/*=========================================================================
 * The finalization routine 
 *=======================================================================*/
void finalize()
{
    // deallocate memory of the filter objects
    endintfilter(lpFilter);
    endintfilter(maFilter);
    endintfilter(deFilter);
    endintfilter(pdFilter1);
    endmaxfilter(maxFilter);
    endmaxfilter(minFilter);
    endintfilter(mdFilter);
    endintfilter(pdFilter2);
}

/*=========================================================================
 * The gateway function 
 *=======================================================================*/
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    mwSize nrows;
    mwSize ncols;
    
    // check argument correctness
    checkArgs(nlhs, plhs, nrhs, prhs);
    
    // handle input arguments
    handleInputs(nrhs, prhs, &nrows, &ncols);
    
    // handle output arguments
    handleOutputs(nlhs, plhs, nrows, ncols);
    
    // calculate the length of the signal
    inputLen = max(nrows,ncols);
    
    // make some initializations
    init();
    
    // do the actual processing
    doTheJob();
    
    // perform final adjustments
    finalize();
}
