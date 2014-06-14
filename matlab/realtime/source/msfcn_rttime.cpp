#define S_FUNCTION_LEVEL 2
#define S_FUNCTION_NAME  msfcn_rttime

#include <simstruc.h>
#include "rttime.h"

#define OBJECT  ((RtTime *)ssGetPWorkValue(S, 0))
#define PARAM   ((double)mxGetPr(ssGetSFcnParam(S, 0))[0])
#define OUTPUT  ((real_T *)ssGetOutputPortSignal(S, 0))[0]

static void mdlInitializeSizes(SimStruct *S)
{
   ssSetNumSFcnParams(S, 1);  /* Number of expected parameters */
   if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) return;
   ssSetNumContStates(S, 0);
   ssSetNumDiscStates(S, 1);
   if (!ssSetNumInputPorts(S, 0)) return;
   if (!ssSetNumOutputPorts(S, 1)) return;
   ssSetOutputPortWidth(S, 0, 1);
   ssSetNumSampleTimes(S, 1);
   ssSetNumRWork(S, 0);
   ssSetNumIWork(S, 0);
   ssSetNumPWork(S, 1);
   ssSetNumModes(S, 0);
   ssSetNumNonsampledZCs(S, 0);
   ssSetOptions(S, 0);
}

static void mdlInitializeSampleTimes(SimStruct *S)
{
   ssSetSampleTime(S, 0, CONTINUOUS_SAMPLE_TIME);
   ssSetOffsetTime(S, 0, 0.0);
}

#define MDL_START
static void mdlStart(SimStruct *S)
{
   ssSetPWorkValue(S, 0, new RtTime(PARAM, ssGetTStart(S)));
}

static void mdlOutputs(SimStruct *S, int_T tid)
{
   RtTime *obj = OBJECT;
   obj->newx(ssGetDiscStates(S), ssGetT(S));
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
