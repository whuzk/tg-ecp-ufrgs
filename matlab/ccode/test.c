/*=========================================================================
 * c_test.c
 *=======================================================================*/
#include "limits.h"
#include "mex.h"

#define MIN_INPUTS  1
#define MAX_INPUTS  1
#define MIN_OUTPUTS 1
#define MAX_OUTPUTS 1

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    unsigned int i = UINT_MAX;
    printf("%u %u %u\n",i,i+1,i+2);
    
    i = 0;
    for (mwSize k = 0; k < 250; k++) {
        printf("%u ",i-k);
    }
    printf("\n");
}
