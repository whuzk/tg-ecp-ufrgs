/*=========================================================================
 * templatebuilder.h
 * 
 *  Title: real-time construction of template and artifact removal
 *  Author:     Diego Sogari
 *  Modified:   12/June/2014
 *
 *=======================================================================*/
#ifndef TEMPLATEBUILDER
#define TEMPLATEBUILDER

#include <string.h>
#include <math.h>
#include "mex.h"
#include "blas.h"

/*=========================================================================
 * Type definitions
 *=======================================================================*/
class TemplateBuilder {
protected:
    mwSize bufferLen;
    int tempCount;
    double *auxBuffer;
    double tempRatio;
    double meanRmsd;
    bool initialized;
    bool isArtifact;
public:
    TemplateBuilder(mwSize bufferLen, int count);
    ~TemplateBuilder();
    void newx(const double *buffer, double *outbuffer);
    bool outputIsArtifact();
    double outputMeanRmsd();
};

/*=========================================================================
 * Constructor
 *=======================================================================*/
TemplateBuilder::TemplateBuilder(mwSize buflen, int count)
{
    this->bufferLen = buflen;
    this->tempCount = count;
    this->tempRatio = 1/(double)count;
    this->meanRmsd = 0.0;
    this->initialized = false;
    this->auxBuffer = new double[buflen];
    memset(auxBuffer, 0, bufferLen * sizeof(double));
}

/*=========================================================================
 * Destructor
 *=======================================================================*/
TemplateBuilder::~TemplateBuilder()
{
    delete[] auxBuffer;
}

/*=========================================================================
 * Update filter memory with an incoming sample
 *=======================================================================*/
void TemplateBuilder::newx(const double *buffer, double *outbuffer)
{
    ptrdiff_t len = bufferLen;
    ptrdiff_t inc = 1;
    double minusone = -1.0;
    double one = 1.0;
    double rmsd;
    
    // calculate rmsd
    memcpy(auxBuffer, buffer, bufferLen * sizeof(double));
    daxpy(&len, &minusone, outbuffer, &inc, auxBuffer, &inc);
    rmsd = dnrm2(&len, auxBuffer, &inc) / sqrt((double)len);
    
    isArtifact = false;
    if (!initialized) {
        // first beat
        memcpy(outbuffer, buffer, bufferLen * sizeof(double));
        initialized = true;
    }
    else if (tempCount > 0 || rmsd < 3 * meanRmsd) {
        // very good beat
        dscal(&len, &tempRatio, auxBuffer, &inc);
        daxpy(&len, &one, auxBuffer, &inc, outbuffer, &inc);
        meanRmsd += tempRatio * (rmsd - meanRmsd);
    }
    else isArtifact = (rmsd > 4 * meanRmsd);
    
    tempCount--;
}

/*=========================================================================
 * Return the output
 *=======================================================================*/
bool TemplateBuilder::outputIsArtifact()
{
    return this->isArtifact;
}

/*=========================================================================
 * Return the output
 *=======================================================================*/
double TemplateBuilder::outputMeanRmsd()
{
    return this->meanRmsd;
}

#endif