/*=========================================================================
 * msfcn_rminmax.c
 * 
 *  Title: S-Function block implementation of running-max/min filters.
 *  Author:     Diego Sogari
 *  Modified:   11/June/2014
 *
 *=======================================================================*/
#define S_FUNCTION_NAME  msfcn_rminmax
#define S_FUNCTION_LEVEL 2

#include "simstruc.h"
#include "minmaxfilter.h"

#define NUM_INPUTS  1
#define NUM_OUTPUTS 1

#define OBJECT  ((MinMaxFilter<int> *)ssGetPWorkValue(S, 0))
#define INPUT   ((const int_T *)ssGetInputPortSignal(S, 0))[0]
#define OUTPUT  ((int_T *)ssGetOutputPortSignal(S, 0))[0]

static void mdlInitializeSizes(SimStruct *S)
{
    int i;
    
    // number of parameters
    ssSetNumSFcnParams(S, 2);
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
        return;
    }
    
    // number of states
    ssSetNumContStates(S, 0);
    ssSetNumDiscStates(S, 0);
    
    // number of ports
    if (!ssSetNumInputPorts(S, NUM_INPUTS)) return;
    if (!ssSetNumOutputPorts(S, NUM_OUTPUTS)) return;
    
    // input port properties
    for (i = 0; i < NUM_INPUTS; i++) {
        ssSetInputPortWidth(S, i, 1);
        ssSetInputPortDirectFeedThrough(S, i, 1);
        ssSetInputPortSampleTime(S, i, INHERITED_SAMPLE_TIME);
        ssSetInputPortRequiredContiguous(S, i, 1);
    }
    ssSetInputPortDataType(S, 0, SS_INT32);
    
    // output port properties
    for (i = 0; i < NUM_OUTPUTS; i++) {
        ssSetOutputPortWidth(S, i, 1);
        ssSetOutputPortSampleTime(S, i, INHERITED_SAMPLE_TIME);
    }
    ssSetOutputPortDataType(S, 0, SS_INT32);
    
    // number of sample times
    ssSetNumSampleTimes(S, 1);
    
    // size of work and mode vectors
    ssSetNumRWork(S, 0);
    ssSetNumIWork(S, 0);
    ssSetNumPWork(S, 1);
    ssSetNumModes(S, 0);
    
    // disable zero-crossing detection
    ssSetNumNonsampledZCs(S, 0);
    
    // exception-free code
    ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE);
    
    // SimState compliance
    ssSetSimStateCompliance(S, USE_DEFAULT_SIM_STATE);
}

static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, INHERITED_SAMPLE_TIME);
    ssSetOffsetTime(S, 0, 0.0);
    ssSetModelReferenceSampleTimeDefaultInheritance(S);  
}

#define MDL_CHECK_PARAMETERS
static void mdlCheckParameters(SimStruct *S)
{
    const mxArray *m;
    
    m = ssGetSFcnParam(S, 0);
    if (mxGetNumberOfElements(m) != 1 || !mxIsNumeric(m) || mxIsComplex(m)) {
        ssSetErrorStatus(S, "First parameter must be real-valued.");
        return;
    }
    else if ((int_T)mxGetPr(m)[0] <= 0) {
        ssSetErrorStatus(S, "The window size must be greater than zero.");
        return;
    }
    
    m = ssGetSFcnParam(S, 1);
    if (mxGetNumberOfElements(m) != 1 || !mxIsLogical(m)) {
        ssSetErrorStatus(S, "Second parameter must be logical.");
        return;
    }
}

#define MDL_START
static void mdlStart(SimStruct *S)
{
    mwSize wsize = (int)mxGetPr(ssGetSFcnParam(S, 0))[0];
    bool ismax = mxGetLogicals(ssGetSFcnParam(S, 1))[0];
    ssSetPWorkValue(S, 0, new MinMaxFilter<int>(wsize, ismax));
}

static void mdlOutputs(SimStruct *S, int_T tid)
{
    MinMaxFilter<int> *obj = OBJECT;
    obj->newx(INPUT);
    OUTPUT = obj->output();
}

static void mdlTerminate(SimStruct *S)
{
    delete OBJECT;
}

#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif
