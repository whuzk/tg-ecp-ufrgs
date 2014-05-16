/*=========================================================================
 * c_template_and_artifacts.c
 * 
 *  Title: construction of a template beat and assertion of artifacts
 *  Author:     Diego Sogari
 *  Modified:   11/May/2014
 *
 *  Intputs:
 *      1. List of beats (in matrix form)*
 *      2. Number of beats for the adaptive template construction
 *
 *  Outputs:
 *      1. List of template versions (in matrix form)*
 *      2. List of flags indicating detected artifacts*
 *      3. Threshold history of the algorithm (in matrix form)
 *
 *  *required
 *
 *=======================================================================*/
#include <math.h>
#include "mex.h"
#include "c_mathutils.h"
#include "c_mexutils.h"
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
 *  IDX     index
 *  FDP     fiducial points
 *=======================================================================*/

/*=========================================================================
 * Constants
 *=======================================================================*/
#define MIN_INPUTS      1       // minimum number of input arguments
#define MAX_INPUTS      2       // maximum number of input arguments
#define MIN_OUTPUTS     2       // minimum number of output arguments
#define MAX_OUTPUTS     3       // maximum number of output arguments

#define FDP_NUM_COL     7       // number of columns in the FDP matrix
#define DEF_NUM_BEATS   30      // default number of beats for the adaptive
                                //      template construction

/*=========================================================================
 * Input and output variables
 *=======================================================================*/
static double *beatList;        // list of beats
static mwSize qrsLen;           // number of beats
static mwSize numBeats;         // number of beats for adaptation
static double *tempList;        // list of template versions
static mxLogical *artFlags;     // list of flags indicating artifacts
static double *thrHist;         // history of thresholds

/*=========================================================================
 * Buffer variables
 *=======================================================================*/
static mwSize bi;               // current beat index
static double *tempBuf;         // buffer for the current template
static double *auxBuf1;         // auxiliary buffer for math operations
static double *auxBuf2;         // second auxiliary buffer
static mwSize frameSize;        // size of one beat
static double thresh[2];        // pair of thresolds
static double ratio;            // ratio of adaptation
static double rmsd;             // value of RMSD between beat and template

/*=========================================================================
 * Fast lookup for buffers
 *=======================================================================*/
#define beat(I,J)   (beatList[(I)*frameSize+(J)])
#define temp(I,J)   (tempList[(I)*frameSize+(J)])
#define thr(I,J)    (thrHist[(I)+qrsLen*(J)])

/*=========================================================================
 * Update outputs
 *=======================================================================*/
void updateOutputs()
{
    if (bi < qrsLen) {
        // save template
        memcpy(&temp(bi, 0), tempBuf, frameSize * sizeof(double));
        // save thresholds
        if (thrHist != NULL) {
            thr(bi,0) = rmsd;
            thr(bi,1) = thresh[0];
            thr(bi,2) = thresh[1];
        }
        // increment beat index
        bi++;
    }
}

/*=========================================================================
 * Event triggered for a new beat
 *=======================================================================*/
void onNewBeat()
{
    // save difference in auxBuf1 and calculate rmsd
    double_subtract(&beat(bi, 0), tempBuf, auxBuf1, frameSize);
    double_array_multiply(auxBuf1, auxBuf1, auxBuf2, frameSize);
    rmsd = sqrt(double_sum(auxBuf2, frameSize) / frameSize);
    
    if (bi == 0) {
        memcpy(tempBuf, &beat(bi, 0), frameSize * sizeof(double));
    }
    else if (bi < numBeats || rmsd < thresh[0]) {
        // very good beat
        double_scalar_multiply(auxBuf1, ratio, auxBuf1, frameSize);
        double_add(tempBuf, auxBuf1, tempBuf, frameSize);
        thresh[0] += ratio * (2 * rmsd - thresh[0]);
        thresh[1] += ratio * (4 * rmsd - thresh[1]);
    }
    else if (rmsd > thresh[1]) {
        // artifact
        artFlags[bi] = true;
    }
    
    // update the outputs
    updateOutputs();
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
            "EcgToolbox:c_template_and_artifacts:nrhs",
            "%d input(s) required and %d optional",
            MIN_INPUTS, MAX_INPUTS - MIN_INPUTS);
    }
    // check for proper number of output arguments
    if (nlhs < MIN_OUTPUTS || nlhs > MAX_OUTPUTS) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_template_and_artifacts:nlhs",
            "%d output(s) required and %d optional",
            MIN_OUTPUTS, MAX_OUTPUTS - MIN_OUTPUTS);
    }
    // make sure all input arguments are of type double
    for (mwSize i = 0; i < nrhs; i++) {
        if (!mxIsDouble(prhs[i]) || mxIsComplex(prhs[i])) {
            mexErrMsgIdAndTxt(
                "EcgToolbox:c_template_and_artifacts:notDouble",
                "Input #%d must be of type double.", i + 1);
        }
    }
    // make sure the first input argument is a matrix
    if (mxGetM(prhs[0]) == 1 || mxGetN(prhs[0]) == 1) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_template_and_artifacts:notMatrix",
            "Input #1 must be a matrix.");
    }
    // make sure the remaining input arguments are all scalars
    for (mwSize i = 1; i < nrhs; i++) {
        if (mxGetNumberOfElements(prhs[i]) != 1) {
            mexErrMsgIdAndTxt(
                "EcgToolbox:c_template_and_artifacts:notScalar",
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
    // get pointers to the data in the input vectors
    beatList = mxGetPr(prhs[0]);
    
    // get the dimensions of the first input
    *nrows = (mwSize)mxGetM(prhs[0]);
    *ncols = (mwSize)mxGetN(prhs[0]);
    
    // get the number of beats for adaptation
    if (nrhs > 1) {
        numBeats = (int)mxGetScalar(prhs[1]);
    }
    else numBeats = DEF_NUM_BEATS;
}

/*=========================================================================
 * Hanlde output arguments 
 *=======================================================================*/
void handleOutputs( int nlhs, mxArray *plhs[])
{
    // get a pointer to the first output
    plhs[0] = mxCreateDoubleMatrix(frameSize, qrsLen, mxREAL);
    tempList = mxGetPr(plhs[0]);
    
    // get a pointer to the second output
    plhs[1] = mxCreateLogicalMatrix(1, qrsLen);
    artFlags = mxGetLogicals(plhs[1]);
    
    // get a pointer to the third output
    if (nlhs > 2) {
        plhs[2] = mxCreateDoubleMatrix(qrsLen, 3, mxREAL);
        thrHist = mxGetPr(plhs[2]);
    }
    else thrHist = NULL;
}

/*=========================================================================
 * The initializarion routine
 *=======================================================================*/
void init()
{
    // initialize beat index
    bi = 0;
    
    // calculate the ratio for adaptation
    ratio = 1.0 / numBeats;
    
    // create the buffers
    tempBuf = (double *)mxCalloc(frameSize, sizeof(double));
    auxBuf1 = (double *)mxMalloc(frameSize * sizeof(double));
    auxBuf2 = (double *)mxMalloc(frameSize * sizeof(double));
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
    for (mwSize i = 0; i < qrsLen; i++) {
        onNewBeat();
    }
    
    // stop time counter
    time = toc();
    
    // display time statistics
    mexPrintf("Total processing time: %.2f ms\n", 1000 * time);
    mexPrintf("Average time per beat: %.2f ns\n", 1000000000 * time / qrsLen);
}

/*=========================================================================
 * The finalization routine 
 *=======================================================================*/
void finalize()
{
    // deallocate memory
    mxFree(tempBuf);
    mxFree(auxBuf1);
    mxFree(auxBuf2);
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
    
    // get the frame size
    frameSize = nrows;
    
    // get the number of beats in the list
    qrsLen = ncols;
    
    // handle output arguments
    handleOutputs(nlhs, plhs);
    
    // make some initializations
    init();
    
    // do the actual processing
    doTheJob();
    
    // perform final adjustments
    finalize();
}
