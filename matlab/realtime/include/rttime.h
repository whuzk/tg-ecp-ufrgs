#ifndef RTTIME
#define RTTIME

#if defined(_WIN32)
    /* Include the windows SDK header for handling time functions. */
    #include <windows.h>
    #include <math.h>
    /* Function of the high performance counter (in seconds). */
    __inline double hightimer()
    {
        HANDLE hCurrentProcess = GetCurrentProcess();
        DWORD dwProcessAffinity, dwSystemAffinity;    
        LARGE_INTEGER frequency, counter;

        /* force thread on first cpu */
        GetProcessAffinityMask(hCurrentProcess,
                (PDWORD_PTR)&dwProcessAffinity,
                (PDWORD_PTR)&dwSystemAffinity);
        SetProcessAffinityMask(hCurrentProcess, 1);
        /* retrieves the frequency of the performance counter */
        QueryPerformanceFrequency(&frequency);
        /* retrieves the current value of the performance counter */
        QueryPerformanceCounter(&counter);
         /* reset thread */
        SetProcessAffinityMask(hCurrentProcess, dwProcessAffinity);
        /* time in seconds */
        return (double)counter.QuadPart / (double)frequency.QuadPart;
    }
#else
    /* Include the standard ANSI C header for handling time functions. */
    #include <time.h>
    /* Function of the high performance counter (in seconds). */
    __inline double hightimer()
    {    
        return (double)clock() / CLOCKS_PER_SEC;
    }
#endif

/*=========================================================================
 * Type definitions
 *=======================================================================*/
class RtTime {
protected:
    double scaleFactor;
    time_T prevSimTime;
    double previousState;
public:
    RtTime(double scale, time_T initial);
    ~RtTime();
    double update(time_T simTime);
};

/*=========================================================================
 * Constructor
 *=======================================================================*/
RtTime::RtTime(double scale, time_T initial)
{
    this->scaleFactor = scale;
    this->prevSimTime = initial;
    this->previousState = hightimer();
}

/*=========================================================================
 * Destructor
 *=======================================================================*/
RtTime::~RtTime()
{
}

/*=========================================================================
 * Update filter memory with an incoming sample
 *=======================================================================*/
double RtTime::update(time_T simTime)
{
   double diff = 0.0;
   double dt, t0;
   double current;
   double previous;
   double elapsed;
   double execution;
   
   /* Desired Delta time */
   dt = (simTime - prevSimTime) * scaleFactor;
   
   /* Get clock time at the beginning of this step */   
   previous = hightimer();
   t0 = previous;
   
   /* Wait to reach the desired time */
   execution = t0 - previousState;
   while (diff < (dt - min(dt, execution))) {
       current = hightimer();
       /* Look for wrapup */
       if (current < previous){
           elapsed = previous - t0;
           t0 = hightimer() - elapsed;
       }
       diff = current - t0;
       previous = current;
   }
   
   /* Store current time to be used in next time step */
   prevSimTime = simTime;
   previousState = previous;
   
   /* Return the output */
   return dt - execution;
}

#endif