/*=========================================================================
 * msfcn_peakevaluator.cpp
 * 
 *  Title: S-Function block implementation of peak evaluation.
 *  Author:     Diego Sogari
 *  Modified:   11/June/2014
 *
 *=======================================================================*/
#define S_FUNCTION_NAME  msfcn_peakevaluator
#define S_FUNCTION_LEVEL 2

#include "simstruc.h"
#include "peakevaluator.h"

#define NUM_INPUTS  6
#define NUM_OUTPUTS 3

#define OBJECT  ((PeakEvaluator<int> *)ssGetPWorkValue(S, 0))
#define INPUT1  ((const int_T *)ssGetInputPortSignal(S, 0))
#define INPUT2  ((const int_T *)ssGetInputPortSignal(S, 1))[0]
#define INPUT3  ((const int_T *)ssGetInputPortSignal(S, 2))[0]
#define INPUT4  ((const bool *)ssGetInputPortSignal(S, 3))[0]
#define INPUT5  ((const bool *)ssGetInputPortSignal(S, 4))[0]
#define INPUT6  ((const bool *)ssGetInputPortSignal(S, 5))[0]
#define OUTPUT1 ((int_T *)ssGetOutputPortSignal(S, 0))[0]
#define OUTPUT2 ((int_T *)ssGetOutputPortSignal(S, 1))[0]
#define OUTPUT3 ((bool *)ssGetOutputPortSignal(S, 2))[0]

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
        if (i != 0) {
            ssSetInputPortWidth(S, i, 1);
        }
        ssSetInputPortDirectFeedThrough(S, i, 1);
        ssSetInputPortSampleTime(S, i, INHERITED_SAMPLE_TIME);
        ssSetInputPortRequiredContiguous(S, i, 1);
    }
    ssSetInputPortFrameData(S, 0, FRAME_YES);
    ssSetInputPortMatrixDimensions(S, 0, DYNAMICALLY_SIZED, DYNAMICALLY_SIZED);
    ssSetInputPortDataType(S, 0, SS_INT32);
    ssSetInputPortDataType(S, 1, SS_INT32);
    ssSetInputPortDataType(S, 2, SS_INT32);
    ssSetInputPortDataType(S, 3, SS_BOOLEAN);
    ssSetInputPortDataType(S, 4, SS_BOOLEAN);
    ssSetInputPortDataType(S, 5, SS_BOOLEAN);
    
    // output port properties
    for (i = 0; i < NUM_OUTPUTS; i++) {
        ssSetOutputPortWidth(S, i, 1);
        ssSetOutputPortSampleTime(S, i, INHERITED_SAMPLE_TIME);
    }
    ssSetOutputPortDataType(S, 0, SS_INT32);
    ssSetOutputPortDataType(S, 1, SS_INT32);
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
    double Fs = mxGetPr(ssGetSFcnParam(S, 0))[0];
    ssSetSampleTime(S, 0, 1.0/Fs);
    ssSetOffsetTime(S, 0, 0.0);
    ssSetModelReferenceSampleTimeDefaultInheritance(S);  
}

#define MDL_CHECK_PARAMETERS
static void mdlCheckParameters(SimStruct *S)
{
    const mxArray *m;
    
    m = ssGetSFcnParam(S, 1);
    if (mxGetNumberOfElements(m) != 1 || !mxIsNumeric(m) || mxIsComplex(m)) {
        ssSetErrorStatus(S, "First parameter must be real-valued.");
        return;
    }
    else if ((int_T)mxGetPr(m)[0] <= 0) {
        ssSetErrorStatus(S, "The adaptation rate must be greater than zero.");
        return;
    }
}

#define MDL_START
static void mdlStart(SimStruct *S)
{
    int length = ssGetInputPortDimensions(S, 0)[0];
    int factor = (int)mxGetPr(ssGetSFcnParam(S, 1))[0];
    ssSetPWorkValue(S, 0, new PeakEvaluator<int>(length, factor));
}

static void mdlOutputs(SimStruct *S, int_T tid)
{
    PeakEvaluator<int> *obj = OBJECT;
    obj->newx(INPUT1, INPUT2, INPUT3, INPUT4, INPUT5, INPUT6);
    OUTPUT1 = obj->outputSignalLevel();
    OUTPUT2 = obj->outputNoiseLevel();
    OUTPUT3 = obj->outputQrsDetected();
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
