/*=========================================================================
 * msfcn_searchpeak.cpp
 * 
 *  Title: S-Function block implementation of peak searching.
 *  Author:     Diego Sogari
 *  Modified:   11/June/2014
 *
 *=======================================================================*/
#define S_FUNCTION_NAME  msfcn_searchpeak
#define S_FUNCTION_LEVEL 2

#include "simstruc.h"
#include "searchpeak.h"

#define NUM_INPUTS  5
#define NUM_OUTPUTS 1
#define NUM_PARAMS  1

#define OBJECT  ((SearchPeak<int> *)ssGetPWorkValue(S, 0))
#define PARAM1  ((double)mxGetPr(ssGetSFcnParam(S, 0))[0])
#define INPUT1  ((const int_T *)ssGetInputPortSignal(S, 0))
#define INPUT2  ((const int_T *)ssGetInputPortSignal(S, 1))
#define INPUT3  ((const int_T *)ssGetInputPortSignal(S, 2))[0]
#define INPUT4  ((const int_T *)ssGetInputPortSignal(S, 3))[0]
#define INPUT5  ((const int_T *)ssGetInputPortSignal(S, 4))[0]
#define OUTPUT  ((int_T *)ssGetOutputPortSignal(S, 0))[0]

static void mdlInitializeSizes(SimStruct *S)
{
    int i;
    
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
    for (i = 0; i < NUM_INPUTS; i++) {
        if (i > 1) {
            ssSetInputPortWidth(S, i, 1);
        }
        ssSetInputPortDirectFeedThrough(S, i, 1);
        ssSetInputPortSampleTime(S, i, INHERITED_SAMPLE_TIME);
        ssSetInputPortRequiredContiguous(S, i, 1);
    }
    ssSetInputPortMatrixDimensions(S, 0, DYNAMICALLY_SIZED, DYNAMICALLY_SIZED);
    ssSetInputPortMatrixDimensions(S, 1, DYNAMICALLY_SIZED, DYNAMICALLY_SIZED);
    ssSetInputPortDataType(S, 0, SS_INT32);
    ssSetInputPortDataType(S, 1, SS_INT32);
    ssSetInputPortDataType(S, 2, SS_INT32);
    ssSetInputPortDataType(S, 3, SS_INT32);
    ssSetInputPortDataType(S, 4, SS_INT32);
    
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
    ssSetSampleTime(S, 0, 1.0/PARAM1);
    ssSetOffsetTime(S, 0, 0.0);
    ssSetModelReferenceSampleTimeDefaultInheritance(S);  
}

#define MDL_START
static void mdlStart(SimStruct *S)
{
    int length = ssGetInputPortDimensions(S, 0)[0];
    ssSetPWorkValue(S, 0, new SearchPeak<int>(length));
}

static void mdlOutputs(SimStruct *S, int_T tid)
{
    SearchPeak<int> *obj = OBJECT;
    obj->newx(INPUT1, INPUT2, INPUT3, INPUT4, INPUT5);
    OUTPUT = obj->outputPeakIdx();
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
