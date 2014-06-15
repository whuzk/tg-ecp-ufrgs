/*=========================================================================
 * rtfilter.h
 * 
 *  Title: real-time implementation of digital IIR filters
 *  Author:     Diego Sogari
 *  Modified:   13/June/2014
 *
 *=======================================================================*/
#ifndef RTFILTER
#define RTFILTER

#include <string.h>
#include <math.h>
#include "mex.h"

/*=========================================================================
 * Macros
 *=======================================================================*/
#define ILOG2(x)    (int)(log10((double)x) / log10(2.0))
#define LENLOG2(x)  (((x) <= 2) ? (x) : 1 << (1 + ILOG2((x) - 1)))
#define XVAL(k)     (xmem.val[(ci+(k))&(xmem.size-1)])
#define YVAL(k)     (ymem.val[(ci+(k))&(ymem.size-1)])

/*=========================================================================
 * Type definitions
 *=======================================================================*/
template <class type>
class RtFilterCoeff {
private:
    mwSize count_nonzeros(const type *coeffs, mwSize len);
    void fill_coeffs(const type *coeffs, mwSize len);
public:
    mwSize size;        // size of the coefficient array
    type *val;          // array of coefficient values
    mwSize *idx;        // indices of the corresponding values
public:
    RtFilterCoeff();
    ~RtFilterCoeff();
    void init(const type *coeffs, mwSize len);
    void print();
};

template <class type>
class RtFilterBuff {
public:
    mwSize size;        // size of the buffer array
    type *val;          // array of buffer values
public:
    RtFilterBuff();
    ~RtFilterBuff();
    void init(mwSize width);
};

template <class type>
class RtFilter {
private:
    unsigned int ci;            // current sample index
protected:
    RtFilterCoeff<type> num;    // numerator coefficients
    RtFilterCoeff<type> den;    // denominator coefficients
    RtFilterBuff<type> xmem;    // X buffer of the filter
    RtFilterBuff<type> ymem;    // Y buffer of the filter
    bool divide;                // divide by leading coefficient
public:
    RtFilter(const type *b, mwSize nb, const type *a, mwSize na,
            bool divide = false);
    ~RtFilter();
    type newx(type x);
    void print();
};

/*=========================================================================
 * Constructor
 *=======================================================================*/
template <class type>
RtFilterCoeff<type>::RtFilterCoeff()
{
    this->size = 0;
    this->val = NULL;
    this->idx = NULL;
}

/*=========================================================================
 * Destructor
 *=======================================================================*/
template <class type>
RtFilterCoeff<type>::~RtFilterCoeff()
{
    delete[] this->val;
    delete[] this->idx;
}

/*=========================================================================
 * Count non zero coefficients
 *=======================================================================*/
template <class type>
mwSize RtFilterCoeff<type>::count_nonzeros(const type *coeffs, mwSize len)
{
    mwSize i, count = 0;
    for (i = 0; i < len; i++) {
        if (abs((double)coeffs[i]) >= (double)1E-10) {
            count++;
        }
    }
    return count;
}

/*=========================================================================
 * Fill in the non-zero portion of the array
 *=======================================================================*/
template <class type>
void RtFilterCoeff<type>::fill_coeffs(const type *coeffs, mwSize len)
{
    mwSize i, count = 0;
    for (i = 0; i < len; i++) {
        if (abs((double)coeffs[i]) >= (double)1E-10) {
            val[count] = coeffs[i];
            idx[count] = i;
            count++;
        }
    }
}

/*=========================================================================
 * Coefficient initialization
 *=======================================================================*/
template <class type>
void RtFilterCoeff<type>::init(const type *coeffs, mwSize len)
{
    this->size = count_nonzeros(coeffs, len);
    this->val = new type[size];
    this->idx = new mwSize[size];
    fill_coeffs(coeffs, len);
}

/*=========================================================================
 * Print the coefficients
 *=======================================================================*/
template <class type>
void RtFilterCoeff<type>::print()
{
    mwSize i;
    printf("[ ");
    for (i = 0; i < size; i++) {
        printf("%.2f ", (double)val[i]);
    }
    printf("] at [ ");
    for (i = 0; i < size; i++) {
        printf("%d ", idx[i]);
    }
    printf("]\n");
}

/*=========================================================================
 * Constructor
 *=======================================================================*/
template <class type>
RtFilterBuff<type>::RtFilterBuff()
{
    this->size = 0;
    this->val = NULL;
}

/*=========================================================================
 * Destructor
 *=======================================================================*/
template <class type>
RtFilterBuff<type>::~RtFilterBuff()
{
    delete[] this->val;
}

/*=========================================================================
 * Buffer initialization
 *=======================================================================*/
template <class type>
void RtFilterBuff<type>::init(mwSize width)
{
    this->size = LENLOG2(width);
    this->val = new type[size];
    memset(val, 0, size * sizeof(type));
}

/*=========================================================================
 * Constructor
 *=======================================================================*/
template <class type>
RtFilter<type>::RtFilter(const type *b, mwSize nb, const type *a,
        mwSize na, bool divide)
{
    this->num.init(b, nb);
    this->den.init(a, na);
    this->xmem.init(nb);
    this->ymem.init(na);
    this->divide = divide;
}

/*=========================================================================
 * Destructor
 *=======================================================================*/
template <class type>
RtFilter<type>::~RtFilter()
{
}

/*=========================================================================
 * Update filter memory with an incoming sample, and return the output.
 *=======================================================================*/
template <class type>
type RtFilter<type>::newx(type x)
{
    type acc = (type)0;
    mwSize i;
    
    // update filter X memory
    XVAL(0) = x;
    
    // compute the first part of the recurrence equation
    for (i = 0; i < num.size; i++) {
        acc += num.val[i] * XVAL(-num.idx[i]);
    }
    
    // compute the second part of the recurrence equation
    for (i = 1; i < den.size; i++) {
        acc -= den.val[i] * YVAL(-den.idx[i]);
    }
    
    // update filter Y memory
    YVAL(0) = acc;
    
    // increment sample index
    ci++;
    
    // report result
    if (divide) {
        // divide by leading coefficient
        return acc / den.val[0];
    }
    else return acc;
}

/*=========================================================================
 * Print the filter coefficients
 *=======================================================================*/
template <class type>
void RtFilter<type>::print()
{
    printf("Numerator:   ");
    num.print();
    printf("Denominator: ");
    den.print();
}

#endif