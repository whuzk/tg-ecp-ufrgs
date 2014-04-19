/*=========================================================================
 * c_detect_qrs.c
 *=======================================================================*/
#include <math.h>
#include "mex.h"

#define DEFAULT_M_FACTOR    3
#define INTEGRATION_LEN_MS  100

static double sigThreshold;     // signal threshold
static double signalLevel;      // signal level
static double noiseLevel;       // noise level
static mwSize qrsCount;         // count of QRS complex in main QRS buffer
static mwSize qrsCount2;        // current position in second QRS buffer
static mwSize searchBackIdx;    // index of searchback starting point
static mwSize lastPeakIdx;      // index of the last peak detected
static double lastPeakAmp;      // amplitude of the last peak detected
static mwSize lastQrsIdx;       // index of last detected QRS complex
static double qrsAmpMean;       // running average qrs amplitude
static double rrIntMean;        // running average of RR intervals
static mwSize rrIntMiss;        // interval for qrs to be assumed as missed
static double levelEstRatio;    // ratio of signal/noise level estimation
static bool isSignalRising;     // flag to indicate a rise in the signal

/* Iterative estimator */
double estimate(double x, double r, double v)
{
    return (1-r)*x + r*v;
}

/* Limit a value by lower and upper bounds */
double limit(double x, double a, double b)
{
    return fmin(fmax(x,a),b);
}

/* Find the position of the maximum value */
mwSize findmax(double *x, mwSize nx)
{
    mwSize pos = 0;
    double y = x[0];
    
    for (mwSize i = 1; i < nx; i++) {
        if (x[i] > y) {
            y = x[i];
            pos = i;
        }
    }
    return pos;
}

/* Get the maximum slope of the signal */
double maxdiff(double *x, mwSize nx)
{
    double d = 0;
    for (mwSize i = 1; i < nx; i++) {
        if (x[i]-x[i-1] > d) {
            d = x[i]-x[i-1];
        }
    }
    return d;
}

/* Check for T wave */
bool istwave(double *x, mwSize candQrs, mwSize lastQrs, int Fs)
{
    // check if the canditate QRS occurs near the previous one
    if (candQrs-lastQrs < (int)round(0.36*Fs)) {
        // length of QRS waves
        mwSize len = (int)round(0.05*Fs);

        // max slope of waveforms
        double slope1 = maxdiff(x+candQrs-len+1, len);
        double slope2 = maxdiff(x+lastQrs-len+1, len);
        
        // check condition for T wave
        return (slope1 < 0.5*slope2);
    }
    else return false;
}

/* Initialize global variables */
void initVariables(int Fs)
{
    sigThreshold = 0.0;
    signalLevel = 0.0;
    noiseLevel = 0.0;
    qrsCount = 0;
    qrsCount2 = 0;
    searchBackIdx = 0;
    lastPeakIdx = 0;
    lastPeakAmp = 0.0;
    lastQrsIdx = 0;
    qrsAmpMean = 0.0;
    rrIntMean = Fs;
    rrIntMiss = (int)round(1.8*Fs);
    levelEstRatio = 0.125;
    isSignalRising = false;
}

/* Peak detection */
bool detectPeak(double *x, mwSize i, mwSize *peakIdx, double *peakAmp)
{
    // check if the new amplitude is greater than that of the last peak
    if (x[i] > lastPeakAmp) {
        // signalize beginning or continuation of positive slope
        isSignalRising = true;
        // update peak info
        lastPeakAmp = x[i];
        lastPeakIdx = i;
        return false;
    }
    else if (!isSignalRising) {
        // update current amplitude
        lastPeakAmp = x[i];
        return false;
    }
    else if (x[i] < 0.5*lastPeakAmp) {
        // report the new peak amplitude and location
        *peakAmp = lastPeakAmp;
        *peakIdx = lastPeakIdx;
        // reset state of positive slope
        isSignalRising = false;
        lastPeakAmp = x[i];
        // signalize detection
        return true;
    }
    else return false;
}

/* QRS detection */
bool detectQrs(double *x, mwSize i, mwSize peakIdx, double peakAmp, int Fs)
{
    if (lastQrsIdx == 0 || peakIdx-lastQrsIdx > (int)round(0.2*Fs)) {
        // decrease estimation ratio for times beyond the training period
        if (i == 2*Fs) {
            levelEstRatio = levelEstRatio*0.5;
        }
        
        // adjust levels
        if (peakAmp >= sigThreshold) {
            signalLevel = estimate(signalLevel,levelEstRatio,peakAmp);
        }
        else {
            noiseLevel = estimate(noiseLevel,levelEstRatio,peakAmp);
        }
        
        // assert qrs
        if (i > 2*Fs && peakAmp >= sigThreshold && !istwave(x, peakIdx, lastQrsIdx, Fs)) {
            return true;
        }
        else return false;
    }
    else return false;
}

/* Search Back procedure */
bool searchBack(double *x, mwSize i, mwSize *peakIdx, double *peakAmp, int Fs)
{
    if (((!isSignalRising) || x[i] < sigThreshold) && searchBackIdx > 0 && i-searchBackIdx >= rrIntMiss) {
        // search back and locate the max in this interval
        mwSize len = (int)round(rrIntMean);
        mwSize begin = searchBackIdx + len/2;
        mwSize pos = findmax(x+begin, len);
        *peakIdx = begin + pos;
        *peakAmp = x[*peakIdx];
        
        // check if candidate peak is from qrs
        if (*peakAmp >= 0.5*sigThreshold && !istwave(x, *peakIdx, lastQrsIdx, Fs)) {
            // adjust signal level
            signalLevel = estimate(signalLevel,0.25,*peakAmp);
            // signalize qrs detection
            return true;
        }
        else {
            // reduce levels by half
            signalLevel = 0.5*signalLevel;
            noiseLevel = 0.5*noiseLevel;
            // postpone searchback
            searchBackIdx += len;
            return false;
        }
    }
    else return false;
}

/* Update QRS and RR interval info */
void updateInfo(mwSize peakIdx, double peakAmp, int Fs)
{
    // update QRS amplitude estimation
    qrsAmpMean = estimate(qrsAmpMean,0.2,peakAmp);

    // update RR interval
    if (lastQrsIdx > 0) {
        mwSize newRR = peakIdx - lastQrsIdx;
        if (0.6*rrIntMean < newRR && newRR < 1.4*rrIntMean) {
            rrIntMean = estimate(rrIntMean,0.2,newRR);
            rrIntMean = limit(rrIntMean,0.2*Fs,2*Fs);
        }
    }

    // calculate RR missed limit
    rrIntMiss = (int)round(1.8*rrIntMean);

    // update indices
    lastQrsIdx = peakIdx;
    searchBackIdx = peakIdx;
}   

/* Main algorithm */
loop( double *x, mwSize i, int Fs, double *qrs, double *qrs2,
      double *thrshist, double *thrnhist, double *rrinthist)
{
    mwSize peakIdx;         // current detected peak index
    double peakAmp;         // current detected peak amplitude
    bool wasPeakDetected;   // flag to indicate that a peak was detected
    bool wasQrsDetected;    // flag to indicate that a qrs was detected
    
    // detect peak
    wasPeakDetected = detectPeak(x, i, &peakIdx, &peakAmp);
    
    // detect qrs
    wasQrsDetected = wasPeakDetected && detectQrs(x, i, peakIdx, peakAmp, Fs);
    
    // search back
    //wasQrsDetected = wasQrsDetected || searchBack(x, i, &peakIdx, &peakAmp, Fs);
    if (!wasQrsDetected) {
        wasQrsDetected = searchBack(x, i, &peakIdx, &peakAmp, Fs);
        if (wasQrsDetected) {
            qrs2[qrsCount2++] = peakIdx+1;
        }
    }
    
    // update info
    if (wasQrsDetected) {
        qrs[qrsCount++] = peakIdx+1;
        updateInfo(peakIdx, peakAmp, Fs);
    }
    
    // update threshold
    sigThreshold = noiseLevel + 0.25*fabs(signalLevel - noiseLevel);
    
    // update history
    thrshist[i] = sigThreshold;
    thrnhist[i] = 0.5*sigThreshold;
    rrinthist[i] = rrIntMean;
}

/* Process input and output arguments */
void processArgs( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[],
                  mwSize *nrows, mwSize *ncols, int *sampFreq)
{
    /* check for proper number of arguments */
    if (nrhs != 2) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_detect_qrs:nrhs",
                "Two inputs required.");
    }
    if (nlhs < 1 || nlhs > 5) {
        mexErrMsgIdAndTxt(
                "EcgToolbox:c_detect_qrs:nlhs",
                "One output required and four optional.");
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
    double *inVector;
    double *outVector1;
    double *outVector2;
    double *outVector3;
    double *outVector4;
    double *outVector5;
    int sampFreq;
    
    mwSize nrows;
    mwSize ncols;
    mwSize inLen;
    
    /* process arguments */
    processArgs(nlhs, plhs, nrhs, prhs, &nrows, &ncols, &sampFreq);
    
    /* calculate the length of the vector */
    inLen = max(nrows,ncols);
    
    /* get a pointer to the data in the input vector  */
    inVector = mxGetPr(prhs[0]);
    
    /* create the output vectors */
    plhs[0] = mxCreateDoubleMatrix(nrows,ncols,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(nrows,ncols,mxREAL);
    plhs[2] = mxCreateDoubleMatrix(nrows,ncols,mxREAL);
    plhs[3] = mxCreateDoubleMatrix(nrows,ncols,mxREAL);
    plhs[4] = mxCreateDoubleMatrix(nrows,ncols,mxREAL);
    
    /* get a pointer to the data in the output vectors */
    outVector1 = mxGetPr(plhs[0]);
    outVector2 = mxGetPr(plhs[1]);
    outVector3 = mxGetPr(plhs[2]);
    outVector4 = mxGetPr(plhs[3]);
    outVector5 = mxGetPr(plhs[4]);
    
    /* initialize global variables */
    initVariables(sampFreq);
    
    /* main algorithm */
    for (mwSize i = 0; i < inLen; i++) {
        loop(inVector, i, sampFreq, outVector1, outVector2,
                outVector3, outVector4, outVector5);
    }
    
    /* adjust size of outputs */
    if (ncols == 1) {
        mxSetM(plhs[0], qrsCount);
        mxSetM(plhs[1], qrsCount2);
    }
    else {
        mxSetN(plhs[0], qrsCount);
        mxSetN(plhs[1], qrsCount2);
    }
    mxSetPr(plhs[0], mxRealloc(outVector1, qrsCount*sizeof(double)));
    mxSetPr(plhs[1], mxRealloc(outVector2, qrsCount2*sizeof(double)));
}
