/*=========================================================================
 * c_noise_filter.c
 * 
 *  Title: process an input signal to remove noise and mains interference
 *  Author:     Diego Sogari
 *  Modified:   10/May/2014
 *
 *  Inputs:
 *      1. input signal*
 *      2. sampling frequency (in Hz)*
 *      3. required output delay (in samples)*
 *      4. mains frequency (in Hz)
 *  Outputs:
 *      1. filtered signal*
 *
 *  *required arguments
 *
 *=======================================================================*/
#include <math.h>
#include "mex.h"
#include "c_fpfilter.h"
#include "c_filtcoeff.h"
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
 *  FREQ    frequency
 *  SIG     signal
 *  FILT    filtered signal
 *  BSF     band-stop filter
 *  LPF     low-pass filter
 *  PDF     pure delay filter
 *=======================================================================*/

/*=========================================================================
 * Constants
 *=======================================================================*/
#define MIN_INPUTS      3       // minimum number of input arguments
#define MAX_INPUTS      4       // maximum number of input arguments
#define MIN_OUTPUTS     1       // minimum number of output arguments
#define MAX_OUTPUTS     1       // maximum number of output arguments

#define MIN_SAMP_FREQ   50*3    // minimum sampling frequency
#define MAX_SAMP_FREQ   60*16   // maximum sampling frequency

/*=========================================================================
 * Input and output variables
 *=======================================================================*/
static double *inputSig;        // input signal vector
static mwSize inputLen;         // input signal length
static double sampFreq;         // sampling frequency
static double mainsFreq;        // mains frequency
static double *filtSig;         // filtered signal
static double delay;            // required output delay

/*=========================================================================
 * Filter variables
 *=======================================================================*/
static unsigned int ci;         // current sample index
static fpfobject bsFilter;      // band-stop filter object
static fpfobject lpFilter;      // low-pass filter object
static fpfobject pdFilter;      // pure delay filter object

/*=========================================================================
 * Low-pass filter
 *=======================================================================*/
double bsf(double sample)
{
    return fpfnewx(&bsFilter, ci, sample);
}

/*=========================================================================
 * High-pass filter
 *=======================================================================*/
double lpf(double sample)
{
    return fpfnewx(&lpFilter, ci, sample);
}

/*=========================================================================
 * Pure delay filter
 *=======================================================================*/
double pdf(double sample)
{
    return fpfnewx(&pdFilter, ci, sample);
}

/*=========================================================================
 * Process one input sample, and return an output sample
 *=======================================================================*/
double processSample(double sample)
{
    sample = bsf(sample);       // band-stop
    sample = lpf(sample);       // low-pass
    return pdf(sample);         // pure delay
}

/*=========================================================================
 * Event triggered for an arriving sample
 *=======================================================================*/
void onNewSample(double sample)
{
    // process sample, and save result
    filtSig[ci++] = processSample(sample);
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
            "EcgToolbox:c_noise_filter:nrhs",
            "%d input(s) required and %d optional",
            MIN_INPUTS, MAX_INPUTS - MIN_INPUTS);
    }
    // check for proper number of output arguments
    if (nlhs < MIN_OUTPUTS || nlhs > MAX_OUTPUTS) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_noise_filter:nlhs",
            "%d output(s) required and %d optional",
            MIN_OUTPUTS, MAX_OUTPUTS - MIN_OUTPUTS);
    }
    // make sure all input arguments are of type double
    for (mwSize i = 0; i < nrhs; i++) {
        if (!mxIsDouble(prhs[i]) || mxIsComplex(prhs[i])) {
            mexErrMsgIdAndTxt(
                "EcgToolbox:c_noise_filter:notDouble",
                "Input #%d must be of type double.", i + 1);
        }
    }
    // make sure the first input argument is a vector
    if (mxGetM(prhs[0]) != 1 && mxGetN(prhs[0]) != 1) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_noise_filter:notVector",
            "First input must be a vector.");
    }
    // make sure the remaining input arguments are all scalars
    for (mwSize i = 1; i < nrhs; i++) {
        if (mxGetNumberOfElements(prhs[i]) != 1) {
            mexErrMsgIdAndTxt(
                "EcgToolbox:c_noise_filter:notScalar",
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
            "EcgToolbox:c_noise_filter:badSampFreq",
            "Sampling frequency must be between %d and %d.",
            MIN_SAMP_FREQ, MAX_SAMP_FREQ);
    }
    
    // get the required delay
    delay = mxGetScalar(prhs[2]);
    
    // get or guess the mains frequency
    if (nrhs > 3) {
        mainsFreq = mxGetScalar(prhs[3]);
    }
    else if ((int)sampFreq % 50 == 0) {
        mainsFreq = 50;
    }
    else {
        mainsFreq = 60;
    }
    
    // make sure the mains frequency is one of the allowed ones
    if ((int)mainsFreq != 50 && (int)mainsFreq != 60) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_noise_filter:badMainsFreq",
            "Mains frequency must be 50 Hz or 60 Hz");
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
    int i = (int)round(sampFreq / mainsFreq) - 3;
    
    // initialize sample index
    ci = 0;
    
    // initialize filter objeects
    initfpfilter(bsFilter);
    initfpfilter(lpFilter);
    initfpfilter(pdFilter);
    
    if ((int)mainsFreq == 50) {
        // create band-stop filter
        fp_create_filter(&bsFilter, chebb50[i], COEFF_COUNT,
                cheba50[i], COEFF_COUNT);

        // design high-pass filter
        fp_create_filter(&lpFilter, buttb50[i], COEFF_COUNT,
                butta50[i], COEFF_COUNT);

        // design derivative filter
        fp_create_allpass(&pdFilter,
                (int)ceil(delay - chebd50[i] - buttd50[i]));
    }
    else {
        // create band-stop filter
        fp_create_filter(&bsFilter, chebb60[i], COEFF_COUNT,
                cheba60[i], COEFF_COUNT);

        // design high-pass filter
        fp_create_filter(&lpFilter, buttb60[i], COEFF_COUNT,
                butta60[i], COEFF_COUNT);

        // design derivative filter
        fp_create_allpass(&pdFilter,
                (int)ceil(delay - chebd60[i] - buttd60[i]));
    }
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
    endfpfilter(bsFilter);
    endfpfilter(lpFilter);
    endfpfilter(pdFilter);
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
