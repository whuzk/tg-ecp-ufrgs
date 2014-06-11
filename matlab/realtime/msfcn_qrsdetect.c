/*=========================================================================
 * msfcn_qrsdetect.c
 * 
 *  Title: S-Function block implementation of QRS detection.
 *  Author:     Diego Sogari
 *  Modified:   11/June/2014
 *
 *=======================================================================*/
#define S_FUNCTION_NAME  msfcn_qrsdetect
#define S_FUNCTION_LEVEL 2

#include "simstruc.h"
#include "qrsdetector.h"

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
    if (!ssSetNumInputPorts(S, 1)) return;
    if (!ssSetNumOutputPorts(S, 2)) return;
    
    // input port properties
    ssSetInputPortWidth(S, 0, 1);
    ssSetInputPortDataType(S, 0, SS_INT32);
    ssSetInputPortDirectFeedThrough(S, 0, 1);
    ssSetInputPortSampleTime(S, 0, INHERITED_SAMPLE_TIME);
    ssSetInputPortRequiredContiguous(S, 0, 1);
    
    // output port properties
    ssSetOutputPortWidth(S, 0, 1);
    ssSetOutputPortDataType(S, 0, SS_INT32);
    ssSetOutputPortSampleTime(S, 0, INHERITED_SAMPLE_TIME);
    ssSetOutputPortWidth(S, 1, 1);
    ssSetOutputPortDataType(S, 1, SS_INT32);
    ssSetOutputPortSampleTime(S, 1, INHERITED_SAMPLE_TIME);
    
    // number of sample times
    ssSetNumSampleTimes(S, 1);
    
    // size of work and mode vectors
    ssSetNumRWork(S, 0);
    ssSetNumIWork(S, 2);
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
    const mxArray *m = ssGetSFcnParam(S, 0);
    if (mxGetNumberOfElements(m) != 1 || !mxIsNumeric(m) || mxIsComplex(m)) {
        ssSetErrorStatus(S, "First parameter must be real-valued.");
        return;
    }
    else if ((int_T)mxGetPr(m)[0] <= 0) {
        ssSetErrorStatus(S, "The sampling frequency must be greater than zero.");
        return;
    }
}

#define MDL_START
static void mdlStart(SimStruct *S)
{
    qrsdetobject *qrsdet;
    const mxArray *m;
    double Fs;
    
    m = ssGetSFcnParam(S, 0);
    Fs = mxGetPr(m)[0];
    
    qrsdet = malloc(sizeof(qrsdetobject));
    initqrsdetector(*qrsdet);
    create_qrsdet(qrsdet, Fs);
    ssSetPWorkValue(S, 0, qrsdet);
}

static void mdlOutputs(SimStruct *S, int_T tid)
{
    int_T *output = ssGetOutputPortSignal(S, 0);
    
    output[0] = ssGetIWorkValue(S, 0);
}

#define MDL_UPDATE
static void mdlUpdate(SimStruct *S, int_T tid)
{
    const int_T *x = ssGetInputPortSignal(S, 0);
    qrsdetobject *qrsdet = ssGetPWorkValue(S, 0);
    mwSize rpeak;
    
    if (qrsdetnewx(qrsdet, *x, &rpeak)) {
        ssSetIWorkValue(S, 0, rpeak);
    }
}

static void mdlTerminate(SimStruct *S)
{
    qrsdetobject *qrsdet = ssGetPWorkValue(S, 0);
    
    endqrsdetector(*qrsdet);
    free(ssGetPWorkValue(S, 0));
}

#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif
