/*=========================================================================
 * c_intfilter.h
 * 
 *  Title: utilities for the design and use of filters with integer
 *         coefficients.
 *  Author:     Diego Sogari
 *  Modified:   07/May/2014
 *
 *  Notes: supported filter types are: all-pass (pure delay), low-pass,
 *         high-pass, band-pass, derivative and moving-average. The design
 *         for a given specification is not always realizable.
 *
 *=======================================================================*/
#ifndef C_INTFILTER
#define C_INTFILTER

#include <stdarg.h>
#include <string.h>
#include <math.h>
#include "mex.h"
#include "c_mathutils.h"

/*=========================================================================
 * Constants
 *=======================================================================*/
#define M_PI    3.1415926535897932384626433832795

/*=========================================================================
 * Type definitions
 *=======================================================================*/
typedef struct {
    mwSize size;        // size of the coefficient array
    int *val;           // array of coefficient values
    mwSize *idx;        // indices of the corresponding coefficients
                        //      in the filter's recurrence equation
} intfcoeffs;

typedef struct {
    mwSize size;        // size of the buffer array
    int *val;           // array of buffer values
} intfbuffer;

typedef struct {
    int gain;           // filter gain (approximated by an integer)
    double delay;       // filter delay (in samples)
    intfcoeffs b;       // numerator coefficients
    intfcoeffs a;       // denominator coefficients
    intfbuffer x;       // X buffer of the filter
    intfbuffer y;       // Y buffer of the filter
} intfobject;

typedef struct {
    double N;           // filter order (except for moving-average)
    double Fc;          // cutoff or center frequency (in Hz)
    double Bw;          // bandwith of bandpass filter (in Hz)
} intfspecs;

/*=========================================================================
 * Macros
 *=======================================================================*/
#define xval(f,k)    (f).x.val[(k)&((f).x.size-1)]
#define yval(f,k)    (f).y.val[(k)&((f).y.size-1)]

#define initintfilter(f)\
    (f).b.size = 0;     \
    (f).b.val = NULL;   \
    (f).b.idx = NULL;   \
    (f).a.size = 0;     \
    (f).a.val = NULL;   \
    (f).a.idx = NULL;   \
    (f).x.size = 0;     \
    (f).x.val = NULL;   \
    (f).y.size = 0;     \
    (f).y.val = NULL;
    
#define endintfilter(f) \
    mxFree((f).b.val);  \
    mxFree((f).b.idx);  \
    mxFree((f).a.val);  \
    mxFree((f).a.idx);  \
    mxFree((f).x.val);  \
    mxFree((f).y.val);

/*=========================================================================
 * Global variables
 *=======================================================================*/
static const char *lp_allowed_specs[] = {
    "N,Fc,3db",         // order and 3db cutoff frequency (in Hz)
    "N,Fc,6db",         // order and 6db cutoff frequency (in Hz)
    "N,Fc,24db",        // order and 24db cutoff frequency (in Hz)
    "N,Fc,nom",         // order and nominal cutoff frequency (in Hz)
    NULL
};
static const char **hp_allowed_specs = lp_allowed_specs;
static const char *bp_allowed_specs[] = {
    "N,Fc,Bw,3db",      // order, center frequency (in Hz) and 3db bandwidth (in Hz)
    "N,Fc,Bw,6db",      // order, center frequency (in Hz) and 6db bandwidth (in Hz)
    "N,Fc,Bw,24db",     // order, center frequency (in Hz) and 24db bandwidth (in Hz)
    "N,Fc1,Fc2,3db",    // order and 3db cutoff frequencies (in Hz)
    "N,Fc1,Fc2,6db",    // order and 6db cutoff frequencies (in Hz)
    "N,Fc1,Fc2,24db",   // order and 24db cutoff frequencies (in Hz)
    NULL
};

/*=========================================================================
 * Update filter memory with an incoming sample, and return the output.
 *=======================================================================*/
int intfnewx(intfobject *filter, unsigned int ci, int sample)
{
    int result = 0;
    
    // update filter X memory
    xval(*filter, ci) = sample;
    
    // compute the first part of the recurrence equation
    for (mwSize i = 0; i < filter->b.size; i++) {
        result += filter->b.val[i] * xval(*filter, ci - filter->b.idx[i]);
    }
    
    // compute the second part of the recurrence equation
    for (mwSize i = 1; i < filter->a.size; i++) {
        result -= filter->a.val[i] * yval(*filter, ci - filter->a.idx[i]);
    }
    
    // update filter Y memory
    yval(*filter, ci) = result;
    
    return result;
}

/*=========================================================================
 * Fill the coefficient object with the non-zero portion of an array. If
 * the number of non-zero elements in the array is different from the size
 * of the coefficient object, this size is adjusted accordingly. The index
 * of the last non-zero coefficient is returned, or zero if no such
 * coefficient exists.
 *=======================================================================*/
mwSize int_construct_coeff(intfcoeffs *coeffObj, const int *a, mwSize na)
{
    mwSize count = na - count_zeros(a, na);
    
    // adjust the size of the coefficient vectors
    if (count != coeffObj->size) {
        coeffObj->val = (int *)mxRealloc(coeffObj->val, count * sizeof(int));
        coeffObj->idx = (mwSize *)mxRealloc(coeffObj->idx, count * sizeof(mwSize));
        coeffObj->size = count;
    }
    
    // fill in the non-zero portion of the array
    while (count-- > 0) {
        while (a[--na] == 0);
        coeffObj->val[count] = a[na];
        coeffObj->idx[count] = na;
    }
    
    return (coeffObj->size > 0) ? coeffObj->idx[coeffObj->size - 1] : 0;
}

/*=========================================================================
 * Reconstruct a coefficient array from a coefficient object. This is the
 * inverse of the function int_construct_coeff. It is assumed that the value
 * of the last element in coeffObj->idx is less than the length of the
 * resulting array.
 *=======================================================================*/
void int_reconstruct_coeff(intfcoeffs *coeffObj, int *result)
{
    for (mwSize i = 0; i < coeffObj->size; i++) {
        result[coeffObj->idx[i]] = coeffObj->val[i];
    }
}

/*=========================================================================
 * Initialize a filter object and allocate memory for its properties
 *=======================================================================*/
void int_create_filter(intfobject *filter, double gain, double delay,
        const int *b, mwSize nb, const int *a, mwSize na)
{
    mwSize xsize, ysize;
    
    // fill in numeric info
    filter->gain = (int)ceil(gain);
    filter->delay = delay;
    
    // initialize the coefficient objects
    xsize = int_construct_coeff(&filter->b, b, nb) + 1;
    ysize = int_construct_coeff(&filter->a, a, na) + 1;
    
    // initialize the X buffer
    if (xsize > 2) {
        xsize = 1 << (1 + ilogb(xsize - 1));
    }
    filter->x.val = (int *)mxRealloc(filter->x.val, xsize * sizeof(int));
    memset(filter->x.val, 0, xsize * sizeof(int));
    filter->x.size = xsize;
    
    // initialize the Y buffer
    if (ysize > 2) {
        ysize = 1 << (1 + ilogb(ysize - 1));
    }
    filter->y.val = (int *)mxRealloc(filter->y.val, ysize * sizeof(int));
    memset(filter->y.val, 0, ysize * sizeof(int));
    filter->y.size = ysize;
}

/*=========================================================================
 * Get the specifications from a specification string
 *=======================================================================*/
int get_specs(char *strSpec, char **strSpecArray, mwSize maxcount)
{
    int count = 0;
    char *chpos;
    
    chpos = strtok(strSpec, ",");
    while (chpos != NULL && count < maxcount) {
        strSpecArray[count++] = chpos;
        chpos = strtok(NULL, ",");
    }
    return count;
}

/*=========================================================================
 * Checks if a string is a member of a string array
 *=======================================================================*/
bool is_member(const char *str, const char **strArray)
{
    mwSize i = 0;
    while (strArray[i] != NULL) {
        if (strcmp(strArray[i++], str) == 0) {
            return true;
        }
    }
    return false;
}

/*=========================================================================
 * Process the specifications
 *=======================================================================*/
void process_specs(va_list vaList, char **strSpecArray,
        mwSize specCount, intfspecs *specObj)
{
    // process each specification
    for (mwSize i = 0; i < specCount; i++) {
        char *strSpec = strSpecArray[i];
        
        if (strcmp(strSpec, "N") == 0) {
            specObj->N = va_arg(vaList, int);
        }
        else if (strcmp(strSpec, "Fc") == 0) {
            specObj->Fc = va_arg(vaList, double);
        }
        else if (strcmp(strSpec, "Fc1") == 0) {
            specObj->Fc = va_arg(vaList, double);
        }
        else if (strcmp(strSpec, "Fc2") == 0) {
            double Fc = va_arg(vaList, double);
            specObj->Bw = (Fc - specObj->Fc);
            specObj->Fc = (Fc + specObj->Fc) / 2.0;
        }
        else if (strcmp(strSpec, "Bw") == 0) {
            specObj->Bw = va_arg(vaList, double);
        }
    }
}

/*=========================================================================
 * Check if a specification string conforms to the expected values
 *=======================================================================*/
void check_spec(const char *specStr, const char **allowedSpecStrings)
{
    if (!is_member(specStr, allowedSpecStrings)) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_intfilter:checkspec:badSpecString",
            "specification string '%s' unknown for this filter type",
            specStr);
    }
}

/*=========================================================================
 * Design a low-pass filter with the specifications provided
 *=======================================================================*/
void get_reciprocal_sense(char *senseStr)
{
    if (strcmp(senseStr, "3db") == 0) {
        strcpy(senseStr, "24db");
    }
    else if (strcmp(senseStr, "24db") == 0) {
        strcpy(senseStr, "3db");
    }
}

/*=========================================================================
 * Get an estimation of the parameter m, according to the sense
 *=======================================================================*/
double get_m(int N, double Wc, char *senseStr)
{
    if (strcmp(senseStr, "3db") == 0) {
        return (0.91823 / sqrt(N + 0.072177)) / Wc;
    }
    else if (strcmp(senseStr, "6db") == 0) {
        return (1.2992 / sqrt(N + 0.15023)) / Wc;
    }
    else if (strcmp(senseStr, "24db") == 0) {
        return (1.732 / sqrt(N + 0.28356)) / Wc;
    }
    else if (strcmp(senseStr, "nom") == 0) {
        return 2.0 / Wc;
    }
}

/*=========================================================================
 * Subtract the filter transfer function from that of an all-pass with the
 * same gain and delay
 *=======================================================================*/
void subtract_allpass(intfobject *filter)
{
    mwSize bsize, asize, tempsize;
    int *btemp, *atemp, *temp1, *temp2;
    
    // recover the coefficient arrays of the filter
    bsize = filter->b.idx[filter->b.size - 1] + 1;
    asize = filter->a.idx[filter->a.size - 1] + 1;
    btemp = (int *)mxCalloc(bsize, sizeof(int));
    atemp = (int *)mxCalloc(asize, sizeof(int));
    int_reconstruct_coeff(&filter->b, btemp);
    int_reconstruct_coeff(&filter->a, atemp);
    
    // create temporary arrays
    tempsize = bsize - asize + 1;
    temp1 = (int *)mxCalloc(tempsize, sizeof(int));
    temp2 = (int *)mxMalloc(bsize * sizeof(int));
    
    // calculate the new coefficients: b = (g * z^(-[d])) * a - b
    filter->delay = ceil(filter->delay);
    temp1[(int)filter->delay] = filter->gain;
    conv(temp1, tempsize, atemp, asize, temp2);
    subtract(temp2, btemp, btemp, bsize);
    int_construct_coeff(&filter->b, btemp, bsize);
    
    // deallocate memory
    mxFree(btemp);
    mxFree(atemp);
    mxFree(temp1);
    mxFree(temp2);
}

/*=========================================================================
 * Design a basic low-pass filter
 *=======================================================================*/
void design_basic_lp(intfobject *filter, int N, int m)
{
    const int abase[2] = {1, -1};
    int *bbase, *btemp, *atemp;
    double gain, delay;
    mwSize bsize, asize;
    
    // check argument correctness
    if (N < 0 || m < 0) return;
    
    // create temporary arrays
    bsize = 1 + N * m;
    asize = 1 + N;
    bbase = (int *)mxCalloc(m + 1, sizeof(int));
    btemp = (int *)mxMalloc(bsize * sizeof(int));
    atemp = (int *)mxMalloc(asize * sizeof(int));
    
    // design the numerator and denominator of the transfer function
    bbase[0] = 1;
    bbase[m] = -1;
    nconv(bbase, m + 1, N, btemp);
    nconv(abase, 2, N, atemp);
    
    // calculate gain and delay
    gain = pow(m, N);
    delay = N * (m - 1) / 2.0;
    
    // initialize the filter object
    int_create_filter(filter, gain, delay, btemp, bsize, atemp, asize);
    
    // deallocate memory
    mxFree(bbase);
    mxFree(btemp);
    mxFree(atemp);
}

/*=========================================================================
 * Design a basic high-pass filter
 *=======================================================================*/
void design_basic_hp(intfobject *filter, int N, int m)
{
    const int abase[2] = {1, 1};
    int *bbase, *btemp, *atemp;
    double gain, delay;
    mwSize bsize, asize;
    
    // check argument correctness
    if (N < 0 || m < 0) return;
    
    // create temporary arrays
    bsize = 1 + N * m;
    asize = 1 + N;
    bbase = (int *)mxCalloc(m + 1, sizeof(int));
    btemp = (int *)mxMalloc(bsize * sizeof(int));
    atemp = (int *)mxMalloc(asize * sizeof(int));
    
    // design the numerator and denominator of the transfer function
    bbase[0] = 1;
    bbase[m] = 1;
    if (m % 2 == 0) {
        bbase[m] = -1;
    }
    nconv(bbase, m + 1, N, btemp);
    nconv(abase, 2, N, atemp);
    
    // calculate gain and delay
    gain = pow(m, N);
    delay = N * (m - 1) / 2.0;
    
    // initialize the filter object
    int_create_filter(filter, gain, delay, btemp, bsize, atemp, asize);
    
    // deallocate memory
    mxFree(bbase);
    mxFree(btemp);
    mxFree(atemp);
}

/*=========================================================================
 * Design a basic band-pass filter
 *=======================================================================*/
void design_basic_bp(intfobject *filter, int N, int m, double theta)
{
    int abase[3] = {1, 0, 1};
    int *bbase, *btemp, *atemp;
    double gain, delay;
    mwSize bsize, asize;
    
    // check argument correctness
    if (N < 0 || m < 0) return;
    
    // create temporary arrays
    bsize = 1 + N * m;
    asize = 1 + N * 2;
    bbase = (int *)mxCalloc(m + 1, sizeof(int));
    btemp = (int *)mxMalloc(bsize * sizeof(int));
    atemp = (int *)mxMalloc(asize * sizeof(int));
    
    // design the numerator and denominator of the transfer function
    bbase[0] = 1;
    bbase[m] = -1;
    nconv(bbase, m + 1, N, btemp);
    abase[1] = (int)round(-2 * cos(theta));
    nconv(abase, 3, N, atemp);
    
    // calculate gain and delay
    gain = m / 2.0 * fabs(cos(m / 2.0 * theta)) / sin(theta);
    delay = N * (m / 2.0 - 1);
    
    // initialize the filter object
    int_create_filter(filter, gain, delay, btemp, bsize, atemp, asize);
    
    // deallocate memory
    mxFree(bbase);
    mxFree(btemp);
    mxFree(atemp);
}

/*=========================================================================
 * Design a basic derivative filter
 *=======================================================================*/
void design_basic_de(intfobject *filter, int N, int M)
{
    const int nbase[2] = {1, -1};
    const int mbase[2] = {1, 1};
    int *ntemp, *mtemp, *btemp;
    double gain, delay;
    mwSize nsize, msize, bsize;
    
    // check argument correctness
    if (N <= 0 || M < 0) return;
    
    // create temporary arrays
    nsize = 1 + N;
    msize = 1 + M;
    bsize = 1 + N + M;
    ntemp = (int *)mxMalloc(nsize * sizeof(int));
    mtemp = (int *)mxMalloc(msize * sizeof(int));
    btemp = (int *)mxMalloc(bsize * sizeof(int));
    
    // design the numerator of the transfer function
    nconv(nbase, 2, N, ntemp);
    nconv(mbase, 2, M, mtemp);
    conv(ntemp, nsize, mtemp, msize, btemp);
    
    // calculate gain and delay
    gain = 1 << max(N, M);
    delay = (N + M) / 2.0;
    
    // initialize the filter object
    int_create_filter(filter, gain, delay, btemp, bsize, NULL, 0);
    
    // deallocate memory
    mxFree(ntemp);
    mxFree(mtemp);
    mxFree(btemp);
}

/*=========================================================================
 * Design an all-pass filter with the specifications provided
 *=======================================================================*/
bool design_allpass(intfobject *filter, int gain, int delay)
{
    int *b;
    
    // check argument correctness
    if (delay < 0) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_intfilter:design_allpass:noSolution",
            "could not design an all-pass for this specfication");
        return false;
    }
    
    // create the numerator of the transfer function
    b = (int *)mxCalloc(delay + 1, sizeof(int));
    b[delay] = gain;
    
    // initialize the filter object
    int_create_filter(filter, gain, delay, b, delay + 1, NULL, 0);
    
    // deallocate memory
    mxFree(b);
    
    return true;
}

/*=========================================================================
 * Design a low-pass filter with the specifications provided
 *=======================================================================*/
bool design_lowpass(intfobject *filter, double Fs, const char *specStr, ...)
{
    va_list vaList;
    char *specArray[3];
    char specCopy[32];
    int specCount;
    intfspecs specObj;
    char *sense;
    double Wc;
    int N, m;
    
    // create a copy of the specification string
    strcpy(specCopy, specStr);
    
    // check if specification is ok
    check_spec(specCopy, lp_allowed_specs);
    
    // get the specification labels
    specCount = get_specs(specCopy, specArray, 3);
    
    // get the string argument
    sense = specArray[--specCount];
    
    // get the numeric arguments
    va_start(vaList, specStr);
    process_specs(vaList, specArray, specCount, &specObj);
    va_end(vaList);
    
    // normalize the cutoff frequency to the interval [0..1]
    Wc = specObj.Fc * 2.0 / Fs;
    
    // get the filter order
    N = (int)specObj.N;
    
    // design the filter
    if (N < 0 || Wc <= 0 || Wc > 1) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_intfilter:design_lowpass:noSolution",
            "could not design a low-pass for this specfication");
        return false;
    }
    else if (Wc == 1) {
        // design the simplest low-pass
        design_basic_lp(filter, N, 2);
    }
    else if (Wc <= 0.5 || strcmp(sense, "nom") == 0) {
        // design low-pass directly
        m = (int)ceil(get_m(N, Wc, sense));
        design_basic_lp(filter, N, m);
    }
    else {
        // subtract high-pass from all-pass
        get_reciprocal_sense(sense);
        m = (int)get_m(N, 1 - Wc, sense);
        if ((m / 2) % 2 != 0) {
            mexErrMsgIdAndTxt(
                "EcgToolbox:c_intfilter:design_lowpass:noSolution",
                "could not design a low-pass for this specfication");
            return false;
        }
        design_basic_hp(filter, N, m);
        subtract_allpass(filter);
    }
    return true;
}

/*=========================================================================
 * Design a high-pass filter with the specifications provided
 *=======================================================================*/
bool design_highpass(intfobject *filter, double Fs, const char *specStr, ...)
{
    va_list vaList;
    char *specArray[3];
    char specCopy[32];
    int specCount;
    intfspecs specObj;
    char *sense;
    double Wc;
    int N, m;
    
    // create a copy of the specification string
    strcpy(specCopy, specStr);
    
    // check if specification is ok
    check_spec(specCopy, hp_allowed_specs);
    
    // get the specification labels
    specCount = get_specs(specCopy, specArray, 3);
    
    // get the string argument
    sense = specArray[--specCount];
    
    // get the numeric arguments
    va_start(vaList, specStr);
    process_specs(vaList, specArray, specCount, &specObj);
    va_end(vaList);
    
    // normalize the cutoff frequency to the interval [0..1]
    Wc = specObj.Fc * 2.0 / Fs;
    
    // get the filter order
    N = (int)specObj.N;
    
    // design the filter
    if (N < 0 || Wc < 0 || Wc >= 1) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_intfilter:design_highpass:noSolution",
            "could not design a high-pass for this specfication");
        return false;
    }
    else if (Wc == 0) {
        // design the simplest high-pass
        design_basic_hp(filter, N, 2);
    }
    else if (Wc >= 0.5 || strcmp(sense, "nom") == 0) {
        // design high-pass directly
        m = (int)ceil(get_m(N, 1 - Wc, sense));
        design_basic_hp(filter, N, m);
    }
    else {
        // subtract low-pass from all-pass
        get_reciprocal_sense(sense);
        m = (int)get_m(N, Wc, sense);
        design_basic_lp(filter, N, m);
        subtract_allpass(filter);
    }
    return true;
}

/*=========================================================================
 * Design a band-pass filter with the specifications provided
 *=======================================================================*/
bool design_bandpass(intfobject *filter, double Fs, const char *specStr, ...)
{
    va_list vaList;
    char *specArray[4];
    char specCopy[32];
    int specCount;
    intfspecs specObj;
    char *sense;
    double Wc, Bw, temp;
    int N, m;
    
    // create a copy of the specification string
    strcpy(specCopy, specStr);
    
    // check if specification is ok
    check_spec(specCopy, bp_allowed_specs);
    
    // get the specification labels
    specCount = get_specs(specCopy, specArray, 4);
    
    // get the string argument
    sense = specArray[--specCount];
    
    // get the numeric arguments
    va_start(vaList, specStr);
    process_specs(vaList, specArray, specCount, &specObj);
    va_end(vaList);
    
    // normalize the frequencies to the interval [0..1]
    Wc = specObj.Fc * 2.0 / Fs;
    Bw = specObj.Bw * 2.0 / Fs;
    
    // get the filter order
    N = (int)specObj.N;
    
    // check if it is possible to design the filter
    if (N < 0 || Bw <= 0 || Bw >= 2) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_intfilter:design_bandpass:noSolution",
            "could not design a band-pass for this specfication");
        return false;
    }
    else if (fabs(Wc - 1.0 / 3.0) <= 0.01) {
        Wc = 1.0 / 3.0;
    }
    else if (fabs(Wc - 1.0 / 2.0) <= 0.01) {
        Wc = 1.0 / 2.0;
    }
    else if (fabs(Wc - 2.0 / 3.0) <= 0.01) {
        Wc = 2.0 / 3.0;
    }
    else {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_intfilter:design_bandpass:noSolution",
            "could not design a band-pass for this specfication");
        return false;
    }
    
    // design the filter
    temp = get_m(N, Bw / 2.0, sense);
    m = (int)round(ceil(temp * Wc / 2.0) / Wc * 2.0);
    design_basic_bp(filter, N, m, M_PI * Wc);
    
    return true;
}

/*=========================================================================
 * Design a derivative filter with the specifications provided
 *=======================================================================*/
bool design_derivative(intfobject *filter, int N, int M)
{
    // design the filter
    if (N <= 0 || M < 0) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_intfilter:design_derivative:noSolution",
            "N must be positive for design of a derivative filter");
        return false;
    }
    else {
        design_basic_de(filter, N, M);
        return true;
    }
}

/*=========================================================================
 * Design a moving-average filter with the specifications provided
 *=======================================================================*/
bool design_maverage(intfobject *filter, double Fs, double width)
{
    // get the m-factor
    int m = (int)round(width * Fs);
    
    // design the filter
    if (m <= 0) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_intfilter:design_maverage:noSolution",
            "could not design a moving-average for this specfication");
        return false;
    }
    else {
        design_basic_lp(filter, 1, m);
        return true;
    }
}

#endif