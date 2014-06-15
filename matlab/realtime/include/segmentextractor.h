/*=========================================================================
 * segmentextractor.h
 * 
 *  Title: extraction of segments
 *  Author:     Diego Sogari
 *  Modified:   12/June/2014
 *
 *=======================================================================*/
#ifndef SEGMENTEXTRACTOR
#define SEGMENTEXTRACTOR

#include <string.h>
#include <math.h>
#include "mex.h"
#include "upfirdn.h"

#define PI          3.1415926535897932384626433832795
#define MAX(a,b)    ((a) > (b) ? (a) : (b))

/*=========================================================================
 * Type definitions
 *=======================================================================*/
class SegmentExtractor {
private:
    void rational(double x, int *p, int *q, double tolerance);
    void impresp(int p, int q);
    void resamp(double *y, mwSize ny, const double *x, mwSize nx);
protected:
    mwSize bufferLen;
    mwSize segmentSize;
    mwSize filterLen;
    mwSize resampLen;
    double *filterBuffer;
    double *resampBuffer;
public:
    SegmentExtractor(mwSize bufferLen, mwSize segsize, mwSize filtlen);
    ~SegmentExtractor();
    void newx(const double *buffer, double *outbuffer, int start, int end);
};

/*=========================================================================
 * Constructor
 *=======================================================================*/
SegmentExtractor::SegmentExtractor(mwSize buflen, mwSize segsize,
        mwSize filtlen)
{
    this->bufferLen = buflen;
    this->segmentSize = segsize;
    this->filterLen = filtlen;
    this->resampLen = segsize + filtlen - 1;
    this->filterBuffer = new double[filtlen];
    this->resampBuffer = new double[resampLen];
}

/*=========================================================================
 * Destructor
 *=======================================================================*/
SegmentExtractor::~SegmentExtractor()
{
    delete[] filterBuffer;
    delete[] resampBuffer;
}

/*=========================================================================
 * Rational approximation of x by p/q, with the specified tolerance
 *=======================================================================*/
void SegmentExtractor::rational(double x, int *p, int *q,
        double tolerance)
{
    int nh[2] = {1, 0};
    int dh[2] = {0, 1};
    int d, saven, saved;
    double savex = x;
    double diff = HUGE_VAL;
    int k = 0;
    
    while (k < 100 && x != 0 && diff > tolerance) {
        k++;
        d = (int)x;
        x = x - d;
        
        saven = nh[0];
        saved = dh[0];
        nh[0] = nh[0] * d + nh[1];
        dh[0] = dh[0] * d + dh[1];
        nh[1] = saven;
        dh[1] = saved;
        
        diff = fabs(nh[0] / (double)dh[0] - savex);
        x = 1 / x;
    }
    
    if (dh[0] < 0) {
        *p = -nh[0];
    }
    else *p = nh[0];
    *q = abs(dh[0]);
}

/*=========================================================================
 * Compute filter impulse response
 *=======================================================================*/
void SegmentExtractor::impresp(int p, int q)
{
    mwSize i;
    double c, win;
    double wc = PI / (double)MAX(p, q);
    double M = filterLen - 1;
    double alpha = M / 2.0;
    
    for (i = 0; i < filterLen; i++) {
        // hamming window
        win = 0.54 - 0.46 * cos(2 * PI * i / M);
        if (i != (int)alpha) {
            // non-center taps
            c = (double)i - alpha;
            filterBuffer[i] = win * p * sin(wc * c) / (PI * c);
        }
        else {
            // center tap
            filterBuffer[i] = win * p * wc / PI;
        }
    }
}

/*=========================================================================
 * Resample
 *=======================================================================*/
void SegmentExtractor::resamp(double *y, mwSize ny, const double *x,
        mwSize nx)
{
    mwSize nr;
    int p, q;
    int half = (filterLen - 1) >> 1;
    
    // compute p and q
    rational(ny / (double)nx, &p, &q, 1.0e-5);
    
    // compute impulse response
    impresp(p, q);
    
    // calculate size of resampling
    nr = (int)ceil(((nx - 1) * p + filterLen) / (double)q);
    if (nr > resampLen) {
        printf("Warning: resampling buffer is not big enough\n");
        return;
    }
    
    // perform the resampling
    memset(resampBuffer, 0, nr * sizeof(double));
    upfirdn(resampBuffer, nr, x, nx, filterBuffer, filterLen, p, q);
    
    // remove delay and write to output
    memcpy(y, &resampBuffer[half / q], ny * sizeof(double));
}

/*=========================================================================
 * Update filter memory with an incoming sample
 *=======================================================================*/
void SegmentExtractor::newx(const double *buffer, double *outbuffer,
        int start, int end)
{
    mwSize len = end - start + 1;
    
    if (len <= 0) {
        // do nothing
        memset(outbuffer, 0, segmentSize * sizeof(double));
    }
    else if (len == segmentSize) {
        // copy input
        memcpy(outbuffer, &buffer[start], segmentSize * sizeof(double));
    }
    else {
        // resample
        resamp(outbuffer, segmentSize, &buffer[start], len);
    }
}

#endif