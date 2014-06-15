/*=========================================================================
 * msfcn_intrtfilter.cpp
 * 
 *  Title: S-Function block implementation of filtering.
 *  Author:     Diego Sogari
 *  Modified:   11/June/2014
 *
 *=======================================================================*/
#define S_FUNCTION_NAME  msfcn_intrtfilter
#define S_FUNCTION_LEVEL 2

#include "simstruc.h"
#include "rtfilter.h"

#define NUM_INPUTS  1
#define NUM_OUTPUTS 1
#define NUM_PARAMS  3

#define OBJECT  ((RtFilter<int> *)ssGetPWorkValue(S, 0))
#define PARAM1  ((int)mxGetPr(ssGetSFcnParam(S, 0))[0])
#define PARAM2a ((int *)mxGetPr(ssGetSFcnParam(S, 1)))
#define PARAM2b ((int)mxGetNumberOfElements(ssGetSFcnParam(S, 1)))
#define PARAM3a ((int *)mxGetPr(ssGetSFcnParam(S, 2)))
#define PARAM3b ((int)mxGetNumberOfElements(ssGetSFcnParam(S, 2)))
#define INPUT1  ((const int *)ssGetInputPortSignal(S, 0))[0]
#define OUTPUT1 ((int *)ssGetOutputPortSignal(S, 0))[0]

static void mdlInitializeSizes(SimStruct *S)
{
    // number of parameters
    ssSetNumSFcnParams(S, NUM_PARAMS);
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
    ssSetInputPortWidth(S, 0, 1);
    ssSetInputPortDirectFeedThrough(S, 0, 1);
    ssSetInputPortRequiredContiguous(S, 0, 1);
    ssSetInputPortSampleTime(S, 0, INHERITED_SAMPLE_TIME);
    ssSetInputPortDataType(S, 0, SS_INT32);
    
    // output port properties
    ssSetOutputPortWidth(S, 0, 1);
    ssSetOutputPortSampleTime(S, 0, INHERITED_SAMPLE_TIME);
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
    ssSetSampleTime(S, 0, 1.0/PARAM1);
    ssSetOffsetTime(S, 0, 0.0);
    ssSetModelReferenceSampleTimeDefaultInheritance(S);  
}

#define MDL_START
static void mdlStart(SimStruct *S)
{
    ssSetPWorkValue(S, 0, new RtFilter<int>(PARAM2a, PARAM2b, PARAM3a,
            PARAM3b));
}

static void mdlOutputs(SimStruct *S, int_T tid)
{
    OUTPUT1 = OBJECT->newx(INPUT1);
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
