/*=========================================================================
 * msfcn_edgedetector.cpp
 * 
 *  Title: S-Function block implementation of block filtering w/o delay
 *  Author:     Diego Sogari
 *  Modified:   11/June/2014
 *
 *=======================================================================*/
#define S_FUNCTION_NAME  msfcn_edgedetector
#define S_FUNCTION_LEVEL 2

#include "simstruc.h"
#include "edgedetector.h"

#define NUM_INPUTS  4
#define NUM_OUTPUTS 1
#define NUM_PARAMS  3

#define OBJECT  ((EdgeDetector<double> *)ssGetPWorkValue(S, 0))
#define PARAM1  ((double)mxGetPr(ssGetSFcnParam(S, 0))[0])
#define PARAM2  ((int)mxGetPr(ssGetSFcnParam(S, 1))[0])
#define PARAM3  ((double)mxGetPr(ssGetSFcnParam(S, 2))[0])
#define INPUT1  ((const real_T *)ssGetInputPortSignal(S, 0))
#define INPUT2  ((const int *)ssGetInputPortSignal(S, 1))[0]
#define INPUT3  ((const int *)ssGetInputPortSignal(S, 2))[0]
#define INPUT4  ((const int *)ssGetInputPortSignal(S, 3))[0]
#define OUTPUT1 ((int *)ssGetOutputPortSignal(S, 0))[0]

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
        if (i > 0) {
            ssSetInputPortWidth(S, i, 1);
        }
        ssSetInputPortDirectFeedThrough(S, i, 1);
        ssSetInputPortSampleTime(S, i, INHERITED_SAMPLE_TIME);
        ssSetInputPortRequiredContiguous(S, i, 1);
    }
    ssSetInputPortMatrixDimensions(S, 0, DYNAMICALLY_SIZED, 1);
    ssSetInputPortDataType(S, 0, SS_DOUBLE);
    ssSetInputPortDataType(S, 1, SS_INT32);
    ssSetInputPortDataType(S, 2, SS_INT32);
    ssSetInputPortDataType(S, 3, SS_INT32);
    
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
    int bufflen = ssGetInputPortDimensions(S, 0)[0];
    ssSetPWorkValue(S, 0, new EdgeDetector<double>(bufflen, PARAM2, PARAM3));
}

static void mdlOutputs(SimStruct *S, int_T tid)
{
    EdgeDetector<double> *obj = OBJECT;
    obj->newx(INPUT1, INPUT2, INPUT3, INPUT4);
    OUTPUT1 = obj->outputEdgeIdx();
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
