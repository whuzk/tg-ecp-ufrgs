/*=========================================================================
 * c_detect_fiducial_marks_double.c
 * 
 *  Title: detect fiducial marks using floating-point arithmetic
 *  Author: Diego Sogari
 *  Last modification:  27/Apr/2014
 *
 *  Outputs:
 *      1. ECG signal*
 *      2. ECG lead name*
 *      3. sampling frequency*
 *      4. mains frequency (is guessed based on sampling frequency)
 *
 *  Outputs:
 *      1. MMD filtered signal
 *      2. LAP filtered signal
 *      3. overall preprocessing delay (in samples)
 *
 *  *required
 *
 *=======================================================================*/
#include <math.h>
#include "mex.h"
#include "c_ecg_utils.h"
#include "c_time_utils.h"

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
 *  LPF     low-pass filter
 *  WMF     windowed-maximum/minimum filter
 *  MMD     multiscale morphological derivative
 *  LAP     discrete laplacian filter
 *  SIG     signal
 *  AMP     amplitude
 *  IDX     index
 *  FILT    filtered signal
 *=======================================================================*/

/*=========================================================================
 * Constants
 *=======================================================================*/
#define MIN_INPUTS      3
#define MAX_INPUTS      4
#define MIN_OUTPUTS     0
#define MAX_OUTPUTS     3

#define MIN_MAINS_FREQ  50
#define MAX_MAINS_FREQ  60
#define DEF_MAINS_FREQ  60
#define DEF_LEAD_NAME   "V4"

#define MMD_LEN_SEC     0.06

#define LPF_BUFLEN      256     // supports sampling frequencies of up to (256*mainsFreq) Hz
#define WMF_BUFLEN      2048    // 2^nextpow2(MAX_MAINS_FREQ * LPF_BUFLEN * MMD_LEN_SEC * 2 + 1)
#define MMD_BUFLEN      1024    // 2^nextpow2(MAX_MAINS_FREQ * LPF_BUFLEN * MMD_LEN_SEC)
#define LAP_BUFLEN      4       // 2^nextpow2(3)

/*=========================================================================
 * Input variables
 *=======================================================================*/
static double *inputSig;        // input signal
static char *leadName;          // ECG lead name
static int sampFreq;            // sampling frequency
static int mainsFreq;           // mains frequency

/*=========================================================================
 * Output variables
 *=======================================================================*/
static double *filtSig1;    // MMD filtered signal
static double *filtSig2;    // LAP filtered signal
static double delay;        // overall preprocessing delay

/*=========================================================================
 * Preprocessing variables
 *=======================================================================*/
static mwSize lpfWin;                       // LPF window length
static int lpfGain;                         // LPF filter gain
static mwSize wmfWin;                       // WMI window length
static mwSize mmdWin;                       // MMD window length
static double lpfBuf[LPF_BUFLEN] = {0};     // LPF buffer
static double wmaBuf[WMF_BUFLEN] = {0};     // WMaxF buffer
static unsigned int wmaAux[WMF_BUFLEN];     // WMaxF aux buffer
static double wmiBuf[WMF_BUFLEN] = {0};     // WMinF buffer
static unsigned int wmiAux[WMF_BUFLEN];     // WMinF aux buffer
static double mmdBuf[MMD_BUFLEN] = {0};     // MMD buffer
static double lapBuf[LAP_BUFLEN] = {0};     // LAP buffer
static unsigned int ci;                     // current sample index

/*=========================================================================
 * Fast lookup for filter buffers
 *=======================================================================*/
#define lpfb(I) (lpfBuf[(ci+(I))&(LPF_BUFLEN-1)])
#define wmab(I) (wmaBuf[(I)&(WMF_BUFLEN-1)])
#define wmaa(I) (wmaAux[(I)&(WMF_BUFLEN-1)])
#define wmib(I) (wmiBuf[(I)&(WMF_BUFLEN-1)])
#define wmia(I) (wmiAux[(I)&(WMF_BUFLEN-1)])
#define mmdb(I) (mmdBuf[(ci+(I))&(MMD_BUFLEN-1)])
#define lapb(I) (lapBuf[(ci+(I))&(LAP_BUFLEN-1)])

/*=========================================================================
 * Low-pass filter
 *=======================================================================*/
double lpf(double sample)
{
    static double y1 = 0.0;
    static double y2 = 0.0;
    
    // update filter output: y[n] = x[n] - 2*x[n-M/2] + x[n-M+1] + 2*y[n-1] - y[n-2]
    double y0 = sample - 2 * lpfb(-lpfWin/2) + lpfb(-lpfWin + 1) + 2 * y1 - y2;
    
    // update filter memory
    lpfb(0) = sample;
    y2 = y1;
    y1 = y0;
    
    // compute result
    return y0 / lpfGain;
}

/*=========================================================================
 * Windowed-maximum
 *=======================================================================*/
double wma(double sample)
{
    static mwSize first = 0;
    static mwSize count = 0;
    mwSize j,k;
    
    // search for the first element greater than the current sample
    j = count;
    k = first + count - 1;
    while (j > 0 && sample >= wmab(k)) {
        j--;
        k--;
    }
    // put the sample next to the element found and adjust the length
    wmab(k + 1) = sample;
    wmaa(k + 1) = ci;
    count = j + 1;
    // check if the first in line has gone out of the window length
    if (count > wmfWin || wmaa(first) == ci - wmfWin) {
        first++;
        count--;
    }
    // return the max in the window
    return wmab(first);
}

/*=========================================================================
 * Windowed-minimum
 *=======================================================================*/
double wmi(double sample)
{
    static mwSize first = 0;
    static mwSize count = 0;
    mwSize j,k;
    
    // search for the first element smaller than the current sample
    j = count;
    k = first + count - 1;
    while (j > 0 && sample <= wmib(k)) {
        j--;
        k--;
    }
    // put the sample next to the element found and adjust the length
    wmib(k + 1) = sample;
    wmia(k + 1) = ci;
    count = j + 1;
    // check if the first in line has gone out of the window length
    if (count > wmfWin || wmia(first) == ci - wmfWin) {
        first++;
        count--;
    }
    // return the max in the window
    return wmib(first);
}

/*=========================================================================
 * Multiscale morphological derivative
 *=======================================================================*/
double mmd(double sample)
{
    // update filter output: y[n] = max{x[n],..,x[n-2*s]} + min{x[n],..,x[n-2*s]} - 2*x[n-s]
    double y0 = wma(sample) + wmi(sample) - 2 * mmdb(-mmdWin);
    
    // update filter memory
    mmdb(0) = sample;
    
    // compute result
    return y0 / mmdWin;
}

/*=========================================================================
 * Discrete laplacian 
 *=======================================================================*/
double lap(double sample)
{
    // update filter output: y[n] = x[n] - 2*x[n-1] + x[n-2]
    double y0 = sample - 2 * lapb(-1) + lapb(-2);
    
    // update filter memory
    lapb(0) = sample;
    
    // compute result
    return y0 / 4;
}

/*=========================================================================
 * Preprocess one input sample 
 *=======================================================================*/
double preprocessSample(double sample, double *mmdval)
{
    static double y1 = 0.0;
    sample = lpf(sample);
    *mmdval = y1;
    y1 = mmd(sample);
    return lap(y1);
}

/*=========================================================================
 * Event triggered for an arriving sample
 *=======================================================================*/
void onNewSample(double sample)
{
    // preprocess sample
    filtSig2[ci] = preprocessSample(sample, &filtSig1[ci]);
    
    // increment global index
    ci++;
}

/*=========================================================================
 * Design the preprocessing filters 
 *=======================================================================*/
void designPreprocessingFilters()
{
    int n;
    
    // low-pass filter
    n = (int)round(sampFreq / mainsFreq);
    lpfWin = (n << 1) + 1;
    lpfGain = n * n;
    
    // windowed-maximum/minimum
    wmfWin = (int)(sampFreq * (double)MMD_LEN_SEC) * 2 + 1;
    
    // multiscale morphlogical derivative
    mmdWin = (int)(sampFreq * (double)MMD_LEN_SEC);
    
    // calculate overall preprocessing delay (LPF + WMF + MMD + LAP)
    delay = n-1 + 0 + mmdWin + 1;
}

/*=========================================================================
 * Check correct number of arguments and their types 
 *=======================================================================*/
void checkArgs( int nlhs, mxArray *plhs[],
                int nrhs, const mxArray *prhs[])
{
    /* check for proper number of input arguments */
    if (nrhs < MIN_INPUTS || nrhs > MAX_INPUTS) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_detect_fiducial_marks_double:nrhs",
                "%d input(s) required.", MIN_INPUTS);
    }
    /* check for proper number of output arguments */
    if (nlhs < MIN_OUTPUTS || nlhs > MAX_OUTPUTS) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_detect_fiducial_marks_double:nlhs",
                "%d output(s) required.", MIN_OUTPUTS);
    }
    /* make sure the first input argument is a vector */
    if (mxGetM(prhs[0]) != 1 && mxGetN(prhs[0]) != 1) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_detect_fiducial_marks_double:notVector",
                "First input must be a Vector.");
    }
    /* make sure the first input argument is of type double */
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_detect_fiducial_marks_double:notDouble",
                "First input must be of type Double.");
    }
    /* make sure the second input argument is a string */
    if (nrhs > 1 && !mxIsChar(prhs[1])) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_detect_fiducial_marks_double:notString",
                "Second input must be a string.");
    }
    /* make sure the remaining arguments are all scalars */
    for (mwSize i = 2; i < nrhs; i++) {
        if (mxGetNumberOfElements(prhs[i]) != 1) {
            mexErrMsgIdAndTxt(
                    "EcgToolbox:c_detect_fiducial_marks_double:notScalar",
                    "Input #%d must be a scalar.", i+1);
        }
    }
}

/*=========================================================================
 * Hanlde input arguments 
 *=======================================================================*/
void handleInputs( int nrhs, const mxArray *prhs[],
                   mwSize *nrows, mwSize *ncols)
{
    size_t charbuflen;
    
    /* get a pointer to the data in the input vector  */
    inputSig = mxGetPr(prhs[0]);
    
    /* get the dimensions of the input vector */
    *nrows = (mwSize)mxGetM(prhs[0]);
    *ncols = (mwSize)mxGetN(prhs[0]);
    
    /* get the ECG lead name  */
    charbuflen = mxGetN(prhs[1]) * sizeof(mxChar) + 1;
    leadName = mxMalloc((mwSize)charbuflen);
    if (mxGetString(prhs[1], leadName, (mwSize)charbuflen) != 0) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_detect_fiducial_marks_double:stringError",
                "Could not read string from input parameter.");
    }
    
    /* get the sampling frequency  */
    sampFreq = (int)mxGetScalar(prhs[2]);
    
    /* get the mains frequency  */
    if (nrhs > 3) {
        mainsFreq = (int)mxGetScalar(prhs[3]);
    }
    else if (sampFreq % 50 == 0) {
        mainsFreq = 50; // guess based on the multiplicity of Fs by 50 Hz
    }
    else if (sampFreq % 60 == 0) {
        mainsFreq = 60; // guess based on the multiplicity of Fs by 60 Hz
    }
    else {
        mainsFreq = DEF_MAINS_FREQ;     // default mains frequency
    }
    
    /* make sure the ECG lead name is one of the known ones */
    if (strcmp(leadName,"V4") != 0 && strcmp(leadName,"MLIII") != 0) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_detect_fiducial_marks_double:badLeadName",
                "ECG lead name '%s' not known.", leadName);
    }
    /* make sure the mains frequency is within pre-defined limits */
    if (mainsFreq < MIN_MAINS_FREQ || mainsFreq > MAX_MAINS_FREQ) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_detect_fiducial_marks_double:badMainsFreq",
                "Mains frequency must be between %d and %d.",
                MIN_MAINS_FREQ, MAX_MAINS_FREQ);
    }
    /* make sure the sampling frequency is within pre-defined limits */
    if (sampFreq < mainsFreq || sampFreq > mainsFreq * LPF_BUFLEN) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_detect_fiducial_marks_double:badSampFreq",
                "Sampling frequency must be between %d and %d.",
                mainsFreq, mainsFreq * LPF_BUFLEN);
    }
}

/*=========================================================================
 * Hanlde output arguments 
 *=======================================================================*/
void handleOutputs( int nlhs, mxArray *plhs[],
                    mwSize nrows, mwSize ncols)
{
    double *outVectors[MAX_OUTPUTS-1] = {NULL};
    
    /* create the output vectors */
    for (mwSize i = 0; i < min(nlhs,MAX_OUTPUTS-1); i++) {
        plhs[i] = mxCreateDoubleMatrix(nrows,ncols,mxREAL);
        outVectors[i] = mxGetPr(plhs[i]);
    }
    
    /* get pointers to the output vectors */
    filtSig1 = outVectors[0];
    filtSig2 = outVectors[1];
}

/*=========================================================================
 * Finalization routine 
 *=======================================================================*/
void finalize(int nlhs, mxArray *plhs[])
{
    /* create the last output (preprocessing delay) */
    if (nlhs > MAX_OUTPUTS-1) {
        plhs[MAX_OUTPUTS-1] = mxCreateDoubleScalar(delay);
    }
    
    // deallocate memory
    mxFree(leadName);
}

/*=========================================================================
 * The gateway function 
 *=======================================================================*/
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    mwSize nrows;
    mwSize ncols;
    mwSize inLen;
    double time;
    
    /* check argument correctness */
    checkArgs(nlhs, plhs, nrhs, prhs);
    
    /* handle input arguments */
    handleInputs(nrhs, prhs, &nrows, &ncols);
    
    /* handle output arguments */
    handleOutputs(nlhs, plhs, nrows, ncols);
    
    /* calculate the length of the signal */
    inLen = max(nrows,ncols);
    
    /* design the preprocessing filters */
    designPreprocessingFilters();
    
    /* process one input sample at a time */
    tic();
    for (mwSize i = 0; i < inLen; i++) {
        onNewSample(inputSig[i]);
    }
    time = toc();
    mexPrintf("Total processing time: %.2f ms\n", 1000*time);
    mexPrintf("Average time per sample: %.2f ns\n", 1000000000*time/inLen);
    
    /* perform some final adjustments */
    finalize(nlhs, plhs);
}
