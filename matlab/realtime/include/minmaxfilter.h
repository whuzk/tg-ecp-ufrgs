/*=========================================================================
 * minmaxfilter.h
 * 
 *  Title: utilities for the use of running-maximum/minimum filters.
 *  Author:     Diego Sogari
 *  Modified:   12/June/2014
 *
 *=======================================================================*/
#ifndef MINMAXFILTER
#define MINMAXFILTER

#include <math.h>
#include "mex.h"

/*=========================================================================
 * Macros
 *=======================================================================*/
#define ILOG2(x)    (int)(log10((double)x) / log10(2.0))
#define VAL(k)      b.val[(first+(k))&(b.size-1)]
#define IDX(k)      b.idx[(first+(k))&(b.size-1)]

/*=========================================================================
 * Types
 *=======================================================================*/
template <class type>
class MinMaxBuff {
public:
    mwSize size;        // size of the buffer array
    type *val;          // array of buffer values
    unsigned int *idx;  // indices of the corresponding values
public:
    MinMaxBuff();
    ~MinMaxBuff();
    void init(mwSize width);
};

template <class type>
class MinMaxFilter {
private:
    unsigned int ci;    // current sample index
protected:
    double delay;       // filter delay (in samples)
    mwSize width;       // width of the filter window
    mwSize count;       // number of elements in the filter buffer
    mwSize first;       // index of the first element in the buffer
    MinMaxBuff<type> b; // filter buffer
    bool ismax;         // flag to indicate if the filter is a max filter
public:
    MinMaxFilter(mwSize width, bool ismax);
    ~MinMaxFilter();
    void newx(type x);
    type output();
};

/*=========================================================================
 * Constructor
 *=======================================================================*/
template <class type>
MinMaxBuff<type>::MinMaxBuff()
{
    this->size = 0;
    this->val = NULL;
    this->idx = NULL;
}

/*=========================================================================
 * Destructor
 *=======================================================================*/
template <class type>
MinMaxBuff<type>::~MinMaxBuff()
{
    delete[] this->val;
    delete[] this->idx;
}

/*=========================================================================
 * Buffer initialization
 *=======================================================================*/
template <class type>
void MinMaxBuff<type>::init(mwSize width)
{
    this->size = (width <= 2) ? width : 1 << (1 + ILOG2(width - 1));
    delete[] this->val;
    delete[] this->idx;
    this->val = new type[size]{0};
    this->idx = new unsigned int[size]{0};
}

/*=========================================================================
 * Constructor
 *=======================================================================*/
template <class type>
MinMaxFilter<type>::MinMaxFilter(mwSize width, bool ismax)
{
    this->delay = width / 2.0;
    this->width = width;
    this->count = 0;
    this->first = 0;
    this->b.init(width);
    this->ismax = ismax;
}

/*=========================================================================
 * Destructor
 *=======================================================================*/
template <class type>
MinMaxFilter<type>::~MinMaxFilter()
{
}

/*=========================================================================
 * Update filter memory with an incoming sample
 *=======================================================================*/
template <class type>
void MinMaxFilter<type>::newx(type x)
{
    mwSize j;
    
    if (ismax) {
        // search for the first element greater than the current sample
        for (j = count; j > 0 && x >= VAL(j - 1); j--);
    }
    else {
        // search for the first element smaller than the current sample
        for (j = count; j > 0 && x <= VAL(j - 1); j--);
    }
    
    // put the sample next to the element found and adjust the count
    VAL(j) = x;
    IDX(j) = ci;
    count = j + 1;
    
    // check if the first in line has gone out of the window length
    if (count > width || IDX(0) == ci - width) {
        first++;
        count--;
    }
    
    // increment sample index
    ci++;
}

/*=========================================================================
 * Return the output
 *=======================================================================*/
template <class type>
type MinMaxFilter<type>::output()
{
    return VAL(0);
}

#endif