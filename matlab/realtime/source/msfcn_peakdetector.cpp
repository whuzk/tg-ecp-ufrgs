/*=========================================================================
 * msfcn_peakdetector.cpp
 * 
 *  Title: S-Function block implementation of peak detection.
 *  Author:     Diego Sogari
 *  Modified:   11/June/2014
 *
 *=======================================================================*/
#define S_FUNCTION_NAME  msfcn_peakdetector
#define S_FUNCTION_LEVEL 2

#include "simstruc.h"
#include "peakdetector.h"

#define NUM_INPUTS  1
#define NUM_OUTPUTS 3

#define OBJECT  ((PeakDetector<int> *)ssGetPWorkValue(S, 0))
#define PARAM1  ((double)mxGetPr(ssGetSFcnParam(S, 0))[0])
#define INPUT   ((const int_T *)ssGetInputPortSignal(S, 0))[0]
#define OUTPUT1 ((int_T *)ssGetOutputPortSignal(S, 0))[0]
#define OUTPUT2 ((bool *)ssGetOutputPortSignal(S, 1))[0]
#define OUTPUT3 ((bool *)ssGetOutputPortSignal(S, 2))[0]

static void mdlInitializeSizes(SimStruct *S)
{
    int i;
    
    // number of parameters
    ssSetNumSFcnParams(S, 1);
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
    ssSetOutputPortDataType(S, 1, SS_BOOLEAN);
    ssSetOutputPortDataType(S, 2, SS_BOOLEAN);
    
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
    ssSetPWorkValue(S, 0, new PeakDetector<int>());
}

static void mdlOutputs(SimStruct *S, int_T tid)
{
    PeakDetector<int> *obj = OBJECT;
    obj->newx(INPUT);
    OUTPUT1 = obj->outputPeakIndex();
    OUTPUT2 = obj->outputPeakDetected();
    OUTPUT3 = obj->outputSignalRising();
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
