/*=========================================================================
 * c_preprocess_filter.c
 * 
 *  Title: process an input ECG signal
 *  Author:     Diego Sogari
 *  Modified:   13/May/2014
 *
 *  Inputs:
 *      1. input signal*
 *      2. sampling frequency (in Hz)*
 *      3. mains frequency (in Hz)
 *  Outputs:
 *      1. filtered signal 1*
 *      2. filtered signal 2*
 *      3. filtered signal 3*
 *      4. filtered signal 4*
 *      5. overall preprocessing delay (in samples)
 *
 *  *required arguments
 *
 *=======================================================================*/
#include <math.h>
#include "mex.h"
#include "c_fpfilter.h"
#include "c_filtcoeff.h"
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
 *  FREQ    frequency
 *  SIG     signal
 *  FILT    filtered signal
 *  LPF     low-pass filter
 *  HPF     high-pass filter
 *  DER     derivative filter
 *  MAF     moving-average filter
 *  MAXF    maximum filter
 *  MINF    minimum filter
 *  MDF     morphological derivative filter
 *  PDF     pure delay filter
 *  BSF     band-stop filter
 *=======================================================================*/

/*=========================================================================
 * Constants
 *=======================================================================*/
#define MIN_INPUTS      2       // minimum number of input arguments
#define MAX_INPUTS      3       // maximum number of input arguments
#define MIN_OUTPUTS     4       // minimum number of output arguments
#define MAX_OUTPUTS     5       // maximum number of output arguments

#define MIN_SAMP_FREQ   50*3    // minimum sampling frequency
#define MAX_SAMP_FREQ   60*16   // maximum sampling frequency

#define TPK_LPF_ORDER   2       // low-pass filter order
#define TPK_LPF_CUTOFF  11      // low-pass filter cutoff frequency (in Hz)
#define TPK_HPF_ORDER   1       // low-pass filter order
#define TPK_HPF_CUTOFF  5       // high-pass filter cutoff frequency (in Hz)
#define TPK_DER_ORDERN  1       // derivative filter order N
#define TPK_DER_ORDERM  3       // derivative filter order M
#define TPK_MAF_WIDTH   0.15    // moving-average width (in seconds)

#define SUN_LPF_ORDER   2       // low-pass filter order
#define SUN_LPF_CUTOFF  11      // low-pass filter cutoff frequency (in Hz)
#define SUN_MAF_WIDTH   0.05    // moving-average filter width (in seconds)
#define SUN_DER_ORDERN  1       // derivative filter order N
#define SUN_DER_ORDERM  0       // derivative filter order M
#define SUN_MDF_WIDTH   0.06    // MD filter width (in seconds)

/*=========================================================================
 * Input and output variables
 *=======================================================================*/
static double *inputSig;        // input signal vector
static mwSize inputLen;         // input signal length
static double sampFreq;         // sampling frequency
static double mainsFreq;        // mains frequency
static double *filtSig1;        // filtered signal 1
static double *filtSig2;        // filtered signal 2
static double *filtSig3;        // filtered signal 3
static double *filtSig4;        // filtered signal 4
static double delay;            // overall delay

/*=========================================================================
 * Filter variables
 *=======================================================================*/
static unsigned int ci;         // current sample index
static intfobject tpkLpFilter;  // low-pass filter object
static intfobject tpkHpFilter;  // high-pass filter object
static intfobject tpkDeFilter;  // derivative filter object
static intfobject tpkMaFilter;  // moving-average filter object
static intfobject sunLpFilter;  // low-pass filter object
static intfobject sunMaFilter;  // moving-average filter object
static intfobject sunDeFilter;  // derivative filter object
static intfobject sunPdFilter1; // pure delay filter object 1
static maxfobject sunMaxFilter; // maximum filter object
static maxfobject sunMinFilter; // minimum filter object
static intfobject sunMdFilter;  // morphological derivative filter object
static intfobject sunPdFilter2; // pure delay filter object 2
static fpfobject noisBsFilter;  // band-stop filter object
static fpfobject noisLpFilter;  // low-pass filter object
static fpfobject noisPdFilter;  // pure delay filter object
static int tpkLpGainLog2;       // log2 of low-pass filter gain
static int tpkHpGainLog2;       // log2 of high-pass filter gain
static int tpkDeGainLog2;       // log2 of derivative filter gain
static int tpkMaGainLog2;       // log2 of moving-average filter gain
static int sunLpGainLog2;       // log2 of low-pass filter gain
static int sunMaGainLog2;       // log2 of moving-average filter gain
static int sunDeGainLog2;       // log2 of derivative filter gain

/*=========================================================================
 * Low-pass filter
 *=======================================================================*/
short tpklpf(short sample)
{
    return intfnewx(&tpkLpFilter, ci, sample) >> tpkLpGainLog2;
}

/*=========================================================================
 * High-pass filter
 *=======================================================================*/
short tpkhpf(short sample)
{
    return intfnewx(&tpkHpFilter, ci, sample) >> tpkHpGainLog2;
}

/*=========================================================================
 * Derivative filter
 *=======================================================================*/
short tpkdef(short sample)
{
    return intfnewx(&tpkDeFilter, ci, sample);// >> tpkDeGainLog2;
}

/*=========================================================================
 * Squaring and hardlimiting
 *=======================================================================*/
short tpksqr(short sample)
{
    int intsample = (int)sample * sample;
    return min(SHRT_MAX, intsample >> 3);
}

/*=========================================================================
 * Moving-average filter
 *=======================================================================*/
int tpkmaf(short sample)
{
    return intfnewx(&tpkMaFilter, ci, sample);// >> tpkMaGainLog2;
}

/*=========================================================================
 * Low-pass filter
 *=======================================================================*/
int sunlpf(int sample)
{
    return intfnewx(&sunLpFilter, ci, sample);// >> lpGainLog2;
}

/*=========================================================================
 * Moving-average filter
 *=======================================================================*/
int sunmaf(int sample)
{
    return intfnewx(&sunMaFilter, ci, sample);// >> maGainLog2;
}

/*=========================================================================
 * Derivative filter
 *=======================================================================*/
int sundef(int sample)
{
    return intfnewx(&sunDeFilter, ci, sample);// >> deGainLog2;
}

/*=========================================================================
 * Pure delay filter 1
 *=======================================================================*/
int sunpdf1(int sample)
{
    return intfnewx(&sunPdFilter1, ci, sample);
}

/*=========================================================================
 * Maximum filter
 *=======================================================================*/
int sunmaxf(int sample)
{
    return maxfnewx(&sunMaxFilter, ci, sample);
}

/*=========================================================================
 * Minimum filter
 *=======================================================================*/
int sunminf(int sample)
{
    return maxfnewx(&sunMinFilter, ci, sample);
}

/*=========================================================================
 * Morphological derivative filter
 *=======================================================================*/
int sunmdf(int sample)
{
    return sunmaxf(sample) + sunminf(sample) -
            (intfnewx(&sunMdFilter, ci, sample) << 1);
}

/*=========================================================================
 * Pure delay filter 2
 *=======================================================================*/
int sunpdf2(int sample)
{
    return intfnewx(&sunPdFilter2, ci, sample);
}

/*=========================================================================
 * Low-pass filter
 *=======================================================================*/
double noisbsf(double sample)
{
    return fpfnewx(&noisBsFilter, ci, sample);
}

/*=========================================================================
 * High-pass filter
 *=======================================================================*/
double noislpf(double sample)
{
    return fpfnewx(&noisLpFilter, ci, sample);
}

/*=========================================================================
 * Pure delay filter
 *=======================================================================*/
double noispdf(double sample)
{
    return fpfnewx(&noisPdFilter, ci, sample);
}

/*=========================================================================
 * Process one input sample, and return an output sample
 *=======================================================================*/
int tpkProcessSample(short sample)
{
    sample = tpklpf(sample);        // low-pass
    sample = tpkhpf(sample);        // high-pass
    sample = tpkdef(sample);        // derivative
    sample = tpksqr(sample);        // squaring and hardlimiting
    return tpkmaf(sample);          // moving-average
}

/*=========================================================================
 * Process one input sample, and return an output sample
 *=======================================================================*/
int sunProcessSample(int sample, int *mdsamp)
{
    sample = sunlpf(sample);        // low-pass
    *mdsamp = sunmdf(sample);       // morphological derivative
    *mdsamp = sunpdf2(*mdsamp);     // delay of morphological derivative
    
    sample = sunmaf(sample);        // moving-average
    sample = sundef(sample);        // derivative
    return sunpdf1(sample);         // delay of derivative
}

/*=========================================================================
 * Process one input sample, and return an output sample
 *=======================================================================*/
double noisProcessSample(double sample)
{
    sample = noisbsf(sample);       // band-stop
    sample = noislpf(sample);       // low-pass
    return noispdf(sample);         // pure delay
}

/*=========================================================================
 * Event triggered for an arriving sample
 *=======================================================================*/
void onNewSample(double sample)
{
    int mdsamp;
    
    // process sample, and save result
    filtSig1[ci] = (double)tpkProcessSample((short)sample);
    filtSig2[ci] = (double)sunProcessSample((int)sample, &mdsamp);
    filtSig4[ci] = noisProcessSample(sample);
    filtSig3[ci] = (double)mdsamp;
    
    // increment sample index
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
            "EcgToolbox:c_preprocess_filter:nrhs",
            "%d input(s) required and %d optional",
            MIN_INPUTS, MAX_INPUTS - MIN_INPUTS);
    }
    // check for proper number of output arguments
    if (nlhs < MIN_OUTPUTS || nlhs > MAX_OUTPUTS) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_preprocess_filter:nlhs",
            "%d output(s) required and %d optional",
            MIN_OUTPUTS, MAX_OUTPUTS - MIN_OUTPUTS);
    }
    // make sure all input arguments are of type double
    for (mwSize i = 0; i < nrhs; i++) {
        if (!mxIsDouble(prhs[i]) || mxIsComplex(prhs[i])) {
            mexErrMsgIdAndTxt(
                "EcgToolbox:c_preprocess_filter:notDouble",
                "Input #%d must be of type double.", i + 1);
        }
    }
    // make sure the first input argument is a vector
    if (mxGetM(prhs[0]) != 1 && mxGetN(prhs[0]) != 1) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_preprocess_filter:notVector",
            "First input must be a vector.");
    }
    // make sure the remaining input arguments are all scalars
    for (mwSize i = 1; i < nrhs; i++) {
        if (mxGetNumberOfElements(prhs[i]) != 1) {
            mexErrMsgIdAndTxt(
                "EcgToolbox:c_preprocess_filter:notScalar",
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
            "EcgToolbox:c_preprocess_filter:badSampFreq",
            "Sampling frequency must be between %d and %d.",
            MIN_SAMP_FREQ, MAX_SAMP_FREQ);
    }
    
    // get or guess the mains frequency
    if (nrhs > 2) {
        mainsFreq = mxGetScalar(prhs[2]);
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
    filtSig1 = mxGetPr(plhs[0]);
    
    // get a pointer to the first output vector
    plhs[1] = mxCreateDoubleMatrix(nrows, ncols, mxREAL);
    filtSig2 = mxGetPr(plhs[1]);
    
    // get a pointer to the second output vector
    plhs[2] = mxCreateDoubleMatrix(nrows, ncols, mxREAL);
    filtSig3 = mxGetPr(plhs[2]);
    
    // get a pointer to the output vector
    plhs[3] = mxCreateDoubleMatrix(nrows, ncols, mxREAL);
    filtSig4 = mxGetPr(plhs[3]);
}

/*=========================================================================
 * Design the filters 
 *=======================================================================*/
void designTpkFilters()
{
    // initialize filter objeects
    initintfilter(tpkLpFilter);
    initintfilter(tpkHpFilter);
    initintfilter(tpkDeFilter);
    initintfilter(tpkMaFilter);
    
    // design low-pass filter
    design_lowpass(&tpkLpFilter, sampFreq, "N,Fc,3db", TPK_LPF_ORDER,
            (double)TPK_LPF_CUTOFF);
    
    // design high-pass filter
    design_highpass(&tpkHpFilter, sampFreq, "N,Fc,3db", TPK_HPF_ORDER,
            (double)TPK_HPF_CUTOFF);
    
    // design derivative filter
    design_derivative(&tpkDeFilter, TPK_DER_ORDERN, TPK_DER_ORDERM);
    
    // design moving-average filter
    design_maverage(&tpkMaFilter, sampFreq, TPK_MAF_WIDTH);
    
    // calculate log2 of filter gains (to speedup division)
    tpkLpGainLog2 = 1 + ilogb(tpkLpFilter.gain - 1);
    tpkHpGainLog2 = 1 + ilogb(tpkHpFilter.gain - 1);
    tpkDeGainLog2 = 1 + ilogb(tpkDeFilter.gain - 1);
    tpkMaGainLog2 = 1 + ilogb(tpkMaFilter.gain - 1);
    
    // calculate overall delay
    delay = tpkLpFilter.delay +
            tpkHpFilter.delay +
            tpkDeFilter.delay +
            tpkMaFilter.delay;
}

/*=========================================================================
 * Design the filters 
 *=======================================================================*/
void designSunFilters()
{
    double delay2;
    int width = (int)round(SUN_MDF_WIDTH * sampFreq);
    
    // initialize filter objects
    initintfilter(sunLpFilter);
    initintfilter(sunMaFilter);
    initintfilter(sunDeFilter);
    initintfilter(sunPdFilter1);
    initmaxfilter(sunMaxFilter);
    initmaxfilter(sunMinFilter);
    initintfilter(sunMdFilter);
    initintfilter(sunPdFilter2);
    
    // design low-pass filter
    design_lowpass(&sunLpFilter, sampFreq, "N,Fc,3db", SUN_LPF_ORDER,
            (double)SUN_LPF_CUTOFF);
    
    // make sure the required delay is sufficiently large
    if (delay < sunLpFilter.delay + width) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_yansun_filter:badDelay",
            "Required delay must be at least equal to %.2f.",
            sunLpFilter.delay + width);
        void finalize();
        finalize();
    }
    
    // design moving-average filter
    design_maverage(&sunMaFilter, sampFreq, SUN_MAF_WIDTH);
    
    // design derivative filter
    design_derivative(&sunDeFilter, SUN_DER_ORDERN, SUN_DER_ORDERM);
    
    // create an all-pass filter for delay of the derivative signal
    delay2 = sunLpFilter.delay + sunMaFilter.delay + sunDeFilter.delay;
    design_allpass(&sunPdFilter1, 1, (int)ceil(delay - delay2));
    
    // create maximum filter
    create_maxfilter(&sunMaxFilter, width * 2 + 1);
    
    // create minimum filter
    create_minfilter(&sunMinFilter, width * 2 + 1);
    
    // create a pure delay filter for the morphological derivative
    design_allpass(&sunMdFilter, 1, width);
    
    // create an all-pass filter for delay of the morphological derivative
    delay2 = sunLpFilter.delay + width;
    design_allpass(&sunPdFilter2, 1, (int)ceil(delay - delay2));
    
    // calculate log2 of filter gains (to speedup division)
    sunLpGainLog2 = 1 + ilogb(sunLpFilter.gain - 1);
    sunMaGainLog2 = 1 + ilogb(sunMaFilter.gain - 1);
    sunDeGainLog2 = 1 + ilogb(sunDeFilter.gain - 1);
}

/*=========================================================================
 * Design the filters 
 *=======================================================================*/
void designNoisFilters()
{
    int i = (int)round(sampFreq / mainsFreq) - 3;
    
    // initialize filter objeects
    initfpfilter(noisBsFilter);
    initfpfilter(noisLpFilter);
    initfpfilter(noisPdFilter);
    
    if ((int)mainsFreq == 50) {
        // create band-stop filter
        fp_create_filter(&noisBsFilter, chebb50[i], COEFF_COUNT,
                cheba50[i], COEFF_COUNT);

        // design high-pass filter
        fp_create_filter(&noisLpFilter, buttb50[i], COEFF_COUNT,
                butta50[i], COEFF_COUNT);

        // design derivative filter
        fp_create_allpass(&noisPdFilter,
                (int)ceil(delay - chebd50[i] - buttd50[i]));
    }
    else {
        // create band-stop filter
        fp_create_filter(&noisBsFilter, chebb60[i], COEFF_COUNT,
                cheba60[i], COEFF_COUNT);

        // design high-pass filter
        fp_create_filter(&noisLpFilter, buttb60[i], COEFF_COUNT,
                butta60[i], COEFF_COUNT);

        // design derivative filter
        fp_create_allpass(&noisPdFilter,
                (int)ceil(delay - chebd60[i] - buttd60[i]));
    }
}

/*=========================================================================
 * The initialization routine 
 *=======================================================================*/
void init()
{
    // initialize sample index
    ci = 0;
    
    // design the filters
    designTpkFilters();
    designSunFilters();
    designNoisFilters();
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
    if (nlhs > 4) {
        plhs[4] = mxCreateDoubleScalar(delay);
    }
    
    // deallocate memory of the filter objects
    endintfilter(tpkLpFilter);
    endintfilter(tpkHpFilter);
    endintfilter(tpkDeFilter);
    endintfilter(tpkMaFilter);
    endintfilter(sunLpFilter);
    endintfilter(sunMaFilter);
    endintfilter(sunDeFilter);
    endintfilter(sunPdFilter1);
    endmaxfilter(sunMaxFilter);
    endmaxfilter(sunMinFilter);
    endintfilter(sunMdFilter);
    endintfilter(sunPdFilter2);
    endfpfilter(noisBsFilter);
    endfpfilter(noisLpFilter);
    endfpfilter(noisPdFilter);
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
