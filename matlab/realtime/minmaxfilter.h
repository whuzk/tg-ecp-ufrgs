/*=========================================================================
 * minmaxfilter.h
 * 
 *  Title: utilities for the use of running-maximum/minimum filters.
 *  Author:     Diego Sogari
 *  Modified:   11/June/2014
 *
 *=======================================================================*/
#ifndef MINMAXFILTER
#define MINMAXFILTER

#include <stdlib.h>
#include <math.h>

/*=========================================================================
 * Type definitions
 *=======================================================================*/
typedef struct {
    mwSize size;        // size of the buffer array
    int *val;           // array of buffer values
    unsigned int *idx;  // indices of the corresponding values
} maxfbuffer;

typedef struct {
    unsigned int ci;    // current sample index
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
#define ILOG2(x)    (int)(log10((double)x) / log10(2.0))
#define VAL(f,k)    (f).b.val[((f).first+(k))&((f).b.size-1)]
#define IDX(f,k)    (f).b.idx[((f).first+(k))&((f).b.size-1)]

#define initmaxfilter(f)\
    (f).b.val = NULL;   \
    (f).b.idx = NULL;
    
#define endmaxfilter(f) \
    free((f).b.val);  \
    free((f).b.idx);

/*=========================================================================
 * Update filter memory with an incoming sample, and return the output.
 *=======================================================================*/
static int maxfnewx(maxfobject *object, int sample)
{
    mwSize j;
    
    if (object->ismax) {
        // search for the first element greater than the current sample
        for (j = object->count; j > 0 && sample >= VAL(*object, j - 1); j--);
    }
    else {
        // search for the first element smaller than the current sample
        for (j = object->count; j > 0 && sample <= VAL(*object, j - 1); j--);
    }
    
    // put the sample next to the element found and adjust the count
    VAL(*object, j) = sample;
    IDX(*object, j) = object->ci;
    object->count = j + 1;
    
    // check if the first in line has gone out of the window length
    if (object->count > object->width ||
            IDX(*object, 0) == object->ci - object->width) {
        object->first++;
        object->count--;
    }
    
    // increment sample index
    object->ci++;
    
    // return the max in the window
    return VAL(*object, 0);
}

/*=========================================================================
 * Initialize a maxfilter object and allocate memory for its properties
 *=======================================================================*/
static void create_minmax(maxfobject *filter, mwSize width, bool ismax)
{
    mwSize size = width;
    if (size > 2) {
        size = 1 << (1 + ILOG2(size - 1));
    }
    filter->delay = width / 2.0;
    filter->width = width;
    filter->count = 0;
    filter->first = 0;
    filter->b.size = size;
    filter->b.val = (int *)realloc(filter->b.val, size * sizeof(int));
    filter->b.idx = (unsigned int *)realloc(filter->b.idx, size * sizeof(unsigned int));
    filter->ismax = ismax;
}

#endif