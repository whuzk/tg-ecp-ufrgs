/*=========================================================================
 * msfcn_smoothbeat.cpp
 * 
 *  Title: S-Function block implementation of smoothing.
 *  Author:     Diego Sogari
 *  Modified:   11/June/2014
 *
 *=======================================================================*/
#define S_FUNCTION_NAME  msfcn_smoothbeat
#define S_FUNCTION_LEVEL 2

#include "simstruc.h"
#include "smoothbeat.h"

#define NUM_INPUTS  1
#define NUM_OUTPUTS 1

#define OBJECT  ((SmoothBeat<double> *)ssGetPWorkValue(S, 0))
#define PARAM1  ((double)mxGetPr(ssGetSFcnParam(S, 0))[0])
#define PARAM2  ((int)mxGetPr(ssGetSFcnParam(S, 1))[0])
#define INPUT1  ((const real_T *)ssGetInputPortSignal(S, 0))
#define OUTPUT1 ((real_T *)ssGetOutputPortSignal(S, 0))

static void mdlInitializeSizes(SimStruct *S)
{
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
    ssSetInputPortDirectFeedThrough(S, 0, 1);
    ssSetInputPortRequiredContiguous(S, 0, 1);
    ssSetInputPortSampleTime(S, 0, INHERITED_SAMPLE_TIME);
    ssSetInputPortFrameData(S, 0, FRAME_YES);
    ssSetInputPortMatrixDimensions(S, 0, DYNAMICALLY_SIZED, 1);
    ssSetInputPortDataType(S, 0, SS_DOUBLE);
    
    // output port properties
    ssSetOutputPortSampleTime(S, 0, INHERITED_SAMPLE_TIME);
    ssSetOutputPortFrameData(S, 0, FRAME_YES);
    ssSetOutputPortMatrixDimensions(S, 0, DYNAMICALLY_SIZED, 1);
    ssSetOutputPortDataType(S, 0, SS_DOUBLE);
    
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

#if defined(MATLAB_MEX_FILE)
#define MDL_SET_INPUT_PORT_DIMENSION_INFO
static void mdlSetInputPortDimensionInfo(
        SimStruct *S, int_T port, const DimsInfo_T *dimsInfo)
{
    if(!ssSetInputPortDimensionInfo(S, port, dimsInfo)) return;
    if (port == 0) {
        ssSetOutputPortDimensionInfo(S, port, dimsInfo);
    }
}

# define MDL_SET_OUTPUT_PORT_DIMENSION_INFO
static void mdlSetOutputPortDimensionInfo(
        SimStruct *S, int_T port, const DimsInfo_T *dimsInfo)
{
    if(!ssSetOutputPortDimensionInfo(S, port, dimsInfo)) return;
}
#endif

#define MDL_START
static void mdlStart(SimStruct *S)
{
    int bufflen = ssGetInputPortDimensions(S, 0)[0];
    ssSetPWorkValue(S, 0, new SmoothBeat<double>(bufflen, PARAM2));
}

static void mdlOutputs(SimStruct *S, int_T tid)
{
    OBJECT->newx(INPUT1, OUTPUT1);
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
