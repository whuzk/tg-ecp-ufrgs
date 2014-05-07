/*=========================================================================
 * c_tompkins_filter.c
 * 
 *  Title: process an input signal using the Tompkins' filters
 *  Author:     Diego Sogari
 *  Modified:   07/May/2014
 *
 *  Inputs:
 *      1. input signal*
 *      2. sampling frequency*
 *  Outputs:
 *      1. filtered signal*
 *      2. overall preprocessing delay (in samples)
 *
 *  *required arguments
 *
 *=======================================================================*/
#include <math.h>
#include "mex.h"
#include "c_intfilter.h"
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
 *  HPF     high-pass filter
 *  DER     derivative filter
 *  MAF     moving-average filter
 *=======================================================================*/

/*=========================================================================
 * Constants
 *=======================================================================*/
#define MIN_INPUTS      2       // minimum number of input arguments
#define MAX_INPUTS      2       // maximum number of input arguments
#define MIN_OUTPUTS     1       // minimum number of output arguments
#define MAX_OUTPUTS     2       // maximum number of output arguments

#define MIN_SAMP_FREQ   50*2    // minimum sampling frequency
#define MAX_SAMP_FREQ   60*40   // maximum sampling frequency

#define LPF_ORDER       2       // low-pass filter order
#define LPF_CUTOFF      11      // low-pass filter cutoff frequency (in Hz)
#define HPF_ORDER       1       // low-pass filter order
#define HPF_CUTOFF      5       // high-pass filter cutoff frequency (in Hz)
#define DER_ORDERN      1       // derivative filter order N
#define DER_ORDERM      3       // derivative filter order M
#define MAF_WIDTH       0.15    // moving-average width (in seconds)

/*=========================================================================
 * Input and output variables
 *=======================================================================*/
static double *inputSig;        // input signal vector
static mwSize inputLen;         // input signal length
static double sampFreq;         // sampling frequency
static double *filtSig;         // filtered signal
static double delay;            // overall delay

/*=========================================================================
 * Filter variables
 *=======================================================================*/
static unsigned int ci;         // current sample index
static intfobject lpFilter;     // low-pass filter object
static intfobject hpFilter;     // high-pass filter object
static intfobject deFilter;     // derivative filter object
static intfobject maFilter;     // moving-average filter object
static int lpGainLog2;          // log2 of low-pass filter gain
static int hpGainLog2;          // log2 of high-pass filter gain
static int deGainLog2;          // log2 of derivative filter gain
static int maGainLog2;          // log2 of moving-average filter gain

/*=========================================================================
 * Low-pass filter
 *=======================================================================*/
short lpf(short sample)
{
    return intfnewx(&lpFilter, ci, sample) >> lpGainLog2;
}

/*=========================================================================
 * High-pass filter
 *=======================================================================*/
short hpf(short sample)
{
    return intfnewx(&hpFilter, ci, sample) >> hpGainLog2;
}

/*=========================================================================
 * Derivative filter
 *=======================================================================*/
short def(short sample)
{
    return intfnewx(&deFilter, ci, sample);// >> deGainLog2;
}

/*=========================================================================
 * Squaring and hardlimiting
 *=======================================================================*/
short sqr(short sample)
{
    int intsample = (int)sample * sample;
    return min(SHRT_MAX, intsample >> 3);
}

/*=========================================================================
 * Moving-average filter
 *=======================================================================*/
int maf(short sample)
{
    return intfnewx(&maFilter, ci, sample);// >> maGainLog2;
}

/*=========================================================================
 * Process one input sample, and return an output sample
 *=======================================================================*/
int processSample(short sample)
{
    sample = lpf(sample);       // low-pass
    sample = hpf(sample);       // high-pass
    sample = def(sample);       // derivative
    sample = sqr(sample);       // squaring and hardlimiting
    return maf(sample);         // moving-average
}

/*=========================================================================
 * Event triggered for an arriving sample
 *=======================================================================*/
void onNewSample(double sample)
{
    // process sample, and save result
    filtSig[ci++] = (double)processSample((short)sample);
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
            "EcgToolbox:c_tompkins_filter:nrhs",
            "%d input(s) required and %d optional",
            MIN_INPUTS, MAX_INPUTS - MIN_INPUTS);
    }
    // check for proper number of output arguments
    if (nlhs < MIN_OUTPUTS || nlhs > MAX_OUTPUTS) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_tompkins_filter:nlhs",
            "%d output(s) required and %d optional",
            MIN_OUTPUTS, MAX_OUTPUTS - MIN_OUTPUTS);
    }
    // make sure the first input argument is a vector
    if (mxGetM(prhs[0]) != 1 && mxGetN(prhs[0]) != 1) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_tompkins_filter:notVector",
            "First input must be a vector.");
    }
    // make sure the first input argument is of type double
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_tompkins_filter:notDouble",
            "First input must be of type double.");
    }
    // make sure the remaining arguments are all scalars
    for (mwSize i = 1; i < nrhs; i++) {
        if (mxGetNumberOfElements(prhs[i]) != 1) {
            mexErrMsgIdAndTxt(
                "EcgToolbox:c_tompkins_filter:notScalar",
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
            "EcgToolbox:c_tompkins_filter:badSampFreq",
            "Sampling frequency must be between %d and %d.",
            MIN_SAMP_FREQ, MAX_SAMP_FREQ);
    }
}

/*=========================================================================
 * Hanlde output arguments 
 *=======================================================================*/
void handleOutputs( int nlhs, mxArray *plhs[],
                    mwSize nrows, mwSize ncols)
{
    // get a pointer to the output vector
    plhs[0] = mxCreateDoubleMatrix(nrows, ncols, mxREAL);
    filtSig = mxGetPr(plhs[0]);
}

/*=========================================================================
 * The initialization routine 
 *=======================================================================*/
void init()
{
    // initialize sample index
    ci = 0;
    
    // initialize filter objeects
    initintfilter(lpFilter);
    initintfilter(hpFilter);
    initintfilter(deFilter);
    initintfilter(maFilter);
    
    // design low-pass filter
    design_lowpass(&lpFilter, sampFreq, "N,Fc,3db", LPF_ORDER,
            (double)LPF_CUTOFF);
    
    // design high-pass filter
    design_highpass(&hpFilter, sampFreq, "N,Fc,3db", HPF_ORDER,
            (double)HPF_CUTOFF);
    
    // design derivative filter
    design_derivative(&deFilter, DER_ORDERN, DER_ORDERM);
    
    // design moving-average filter
    design_maverage(&maFilter, sampFreq, MAF_WIDTH);
    
    // calculate log2 of filter gains (to speedup division)
    lpGainLog2 = 1 + ilogb(lpFilter.gain - 1);
    hpGainLog2 = 1 + ilogb(hpFilter.gain - 1);
    deGainLog2 = 1 + ilogb(deFilter.gain - 1);
    maGainLog2 = 1 + ilogb(maFilter.gain - 1);

    // calculate overall delay
    delay = lpFilter.delay +
            hpFilter.delay +
            deFilter.delay +
            maFilter.delay;
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
void finalize( int nlhs, mxArray *plhs[],
               mwSize nrows, mwSize ncols)
{
    // assign the overall delay to the corresponding output
    if (nlhs > 1) {
        plhs[1] = mxCreateDoubleScalar(delay);
    }
    
    // deallocate memory of the filter objects
    endintfilter(lpFilter);
    endintfilter(hpFilter);
    endintfilter(deFilter);
    endintfilter(maFilter);
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
    finalize(nlhs, plhs, nrows, ncols);
}
