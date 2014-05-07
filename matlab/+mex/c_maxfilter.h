/*=========================================================================
 * c_maxfilter.h
 * 
 *  Title: utilities for the use of running-maximum/minimum filters.
 *  Author:     Diego Sogari
 *  Modified:   07/May/2014
 *
 *=======================================================================*/
#ifndef C_MAXFILTER
#define C_MAXFILTER

/*=========================================================================
 * Type definitions
 *=======================================================================*/
typedef struct {
    mwSize size;        // size of the buffer array
    int *val;           // array of buffer values
    unsigned int *idx;  // indices of the corresponding values
} maxfbuffer;

typedef struct {
    double delay;       // filter delay (in samples)
    mwSize width;       // width of the filter window
    mwSize count;       // number of elements in the filter buffer
    mwSize first;       // index of the first element in the buffer
    maxfbuffer b;       // filter buffer
    bool ismax;         // flag to indicate whether the filter is a max
                        //      filter, rather than a min filter
} maxfobject;

/*=========================================================================
 * Macros
 *=======================================================================*/
#define val(f,k)    (f).b.val[((f).first+(k))&((f).b.size-1)]
#define idx(f,k)    (f).b.idx[((f).first+(k))&((f).b.size-1)]

#define initmaxfilter(f)\
    (f).b.val = NULL;   \
    (f).b.idx = NULL;
    
#define endmaxfilter(f) \
    mxFree((f).b.val);  \
    mxFree((f).b.idx);

/*=========================================================================
 * Update filter memory with an incoming sample, and return the output.
 *=======================================================================*/
int maxfnewx(maxfobject *filter, unsigned int ci, int sample)
{
    mwSize j;
    
    if (filter->ismax) {
        // search for the first element greater than the current sample
        for (j = filter->count; j > 0 && sample >= val(*filter, j - 1); j--);
    }
    else {
        // search for the first element smaller than the current sample
        for (j = filter->count; j > 0 && sample <= val(*filter, j - 1); j--);
    }
    
    // put the sample next to the element found and adjust the count
    val(*filter, j) = sample;
    idx(*filter, j) = ci;
    filter->count = j + 1;
    
    // check if the first in line has gone out of the window length
    if (filter->count > filter->width ||
            idx(*filter, 0) == ci - filter->width) {
        filter->first++;
        filter->count--;
    }
    
    // return the max in the window
    return val(*filter, 0);
}

/*=========================================================================
 * Initialize a maxfilter object and allocate memory for its properties
 *=======================================================================*/
void create_minmax(maxfobject *filter, mwSize width, bool ismax)
{
    mwSize size = width;
    if (size > 2) {
        size = 1 << (1 + ilogb(size - 1));
    }
    filter->delay = width / 2.0;
    filter->width = width;
    filter->count = 0;
    filter->first = 0;
    filter->b.size = size;
    filter->b.val = (int *)mxRealloc(filter->b.val, size * sizeof(int));
    filter->b.idx = (unsigned int *)mxRealloc(filter->b.idx, size * sizeof(unsigned int));
    filter->ismax = ismax;
}

/*=========================================================================
 * Create a maximum filter with the specified width (in samples)
 *=======================================================================*/
bool create_maxfilter(maxfobject *filter, int width)
{
    // check argument correctness
    if (width <= 0) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_maxfilter:create_maxfilter:noSolution",
            "could not create a max filter for this specfication");
        return false;
    }
    
    // initialize the filter object
    create_minmax(filter, width, true);
    
    return true;
}

/*=========================================================================
 * Create a minimum filter with the specified width (in samples)
 *=======================================================================*/
bool create_minfilter(maxfobject *filter, int width)
{
    // check argument correctness
    if (width <= 0) {
        mexErrMsgIdAndTxt(
            "EcgToolbox:c_maxfilter:create_minfilter:noSolution",
            "could not create a min filter for this specfication");
        return false;
    }
    
    // initialize the filter object
    create_minmax(filter, width, false);
    
    return true;
}

#endif