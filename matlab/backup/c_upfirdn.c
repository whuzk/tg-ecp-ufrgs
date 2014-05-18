/*=========================================================================
 *
 *  c_upfirdn.c   a MEX-file to perform multirate filtering
 *
 *  The calling syntax is:
 *      y = upfirdn(x, h, p, q)
 *
 *  Jim McClellan    21-Jan-95
 *
 *  $Revision: 1.9 $
 *  Copyright (c) 1988-98 by The MathWorks, Inc.
 *
 *  On PC platforms, must be compiled with /Alfw and /Gs option to generate
 *  proper code.
 *
 *=======================================================================*/
#include <math.h>   
#include "mex.h"   

/*=========================================================================
 * upfirdn
 *=======================================================================*/
static void upfirdn(double y[],  unsigned int Ly, 
                    double x[],  unsigned int Lx,
                    double h[],  unsigned int Lh,
                    int p, int q)   
{   
    int r, rpq_offset, Lg;
    int iv, ig, igv, iw;
    double  *pw;
    double  *pv, *pvend;
    double  *pvhi, *pvlo, *pvt;
    double  *pg, *pgend;
    double  *pghi, *pglo, *pgt;
    
    iv  = q;   
    ig  = iw = p;   
    igv = p*q;   
   
    pvend = x + Lx;   
    pgend = h + Lh;   
   
    for (r=0; r<p; r++) {   
        pw = y + r;   
        pg = h + ( (r*q)%p );   
        Lg = (int)(pgend - pg);   
        Lg = (Lg%p) ? Lg/p+1 : Lg/p ;   
        rpq_offset = (r*q)/p;   
        pv = x + rpq_offset;   
   
        /*  
         * PSEUDO-CODE for CONVOLUTION with GENERAL INCREMENTS:  
         *  
         *   w[n] = v[n] * g[n]  
         *  
         * Given:  
         *   pointers:   pg, pv, and pw  
         *   or arrays:  g[ ], v[ ], and w[ ]  
         *   increments: ig, iv, and iw  
         *   end points: h+Lh, x+Lx  
         */   
   
        /*  
         * Region #1 (running onto the data):  
         */   
        pglo = pg;   
        pghi = pg + p*rpq_offset;   
        pvlo = x;   
        pvhi = pv;   
        while ((pvhi<pvend) && (pghi<pgend)) {   
            double acc = 0.0;   
            pvt = pvhi;   
            pgt = pglo;   
            while (pgt <= pghi) {   
                acc += (*pgt) * (*pvt--);   
                pgt += ig;   
            }   
            *pw  += acc;   
            pw   += iw;   
            pvhi += iv;   
            pghi += igv;   
        }   
   
        /*  
         * Do we need to drain rest of signal?  
         */   
        if (pvhi < pvend)  {   
            /*  
             * Region #2 (complete overlap):  
             */   
            while (pghi >= pgend) {   
                pghi -= ig;   
            }   
            while (pvhi < pvend)  {   
                double acc = 0.0;   
                pvt = pvhi;   
                pgt = pglo;   
                while (pgt <= pghi) {   
                    acc += (*pgt) * (*pvt--);   
                    pgt += ig;   
                }   
                *pw  += acc;   
                pw   += iw;   
                pvhi += iv;   
            }   
        }
        else if (pghi < pgend)  {   
            /*  
             * Region #2a (drain out the filter):  
             */   
            while (pghi < pgend)  {   
                double acc = 0.0;   
                pvt = pvlo;     /* pvlo is still equal to x */   
                pgt = pghi;   
                while( pvt < pvend ) {   
                    acc += (*pgt) * (*pvt++);   
                    pgt -= ig;   
                }   
                *pw += acc;   
                pw += iw;   
                pghi += igv;   
                pvhi += iv;   
            }   
        }   
   
        while (pghi >= pgend) {   
            pghi -= ig;   
        }   
        pvlo = pvhi - Lg + 1;   
 
        while (pvlo < pvend)  {   
            /*  
             *  Region #3 (running off the data):  
             */   
            double acc = 0.0;   
            pvt = pvlo;   
            pgt = pghi;   
            while (pvt < pvend) {   
                acc += (*pgt) * (*pvt++);   
                pgt -= ig;   
            }   
            *pw += acc;   
            pw += iw;   
            pvlo += iv;   
        }   
    } 
}

/*=========================================================================
 * The gateway function
 *=======================================================================*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *yy, *xx, *hh;
    int pp, qq;
    size_t Lx, Lh, Ly, mx, mh, my, nx, nh, ny;
   
    // Check the dimensions of X
    mx = mxGetM(prhs[0]);   
    nx = mxGetN(prhs[0]); 
    Lx = max(nx,mx);
    
    // Check the dimensions of H
    mh = mxGetM(prhs[1]);
    nh = mxGetN(prhs[1]);
    Lh = max(nh, mh);
   
    // Get pointers to input and filter
    xx = mxGetPr(prhs[0]);
    hh = mxGetPr(prhs[1]);
    
    // Supply defaults for P and Q
    pp = (int)mxGetScalar(prhs[2]);
    qq = (int)mxGetScalar(prhs[3]);
    if (pp < 1 || qq < 1) {   
        mexErrMsgTxt("P and Q must be positive integers");   
    } 
    
    // Create return array
    Ly = pp * (Lx - 1) + Lh;
    Ly = (Ly % qq) ? Ly / qq + 1 : Ly / qq;
    plhs[0] = mxCreateDoubleMatrix((mwSize)Ly, 1, mxREAL);
    yy = mxGetPr(plhs[0]);
    
    // Input and filters are both purely real
    upfirdn(yy, Ly, xx, Lx, hh, Lh, pp, qq);
}