#ifndef yeti_sortimpl_h
#define yeti_sortimpl_h

#include "sort.h"
#include "permutation.h"

namespace yeti {

template <class T>
static void
qs(T* item, int* index, int left, int right)
{
    register int i,j;
    T x;
    int y;

    i=left; j=right;
    x=item[index[(left+right)/2]];

    do {
        while(item[index[i]]<x && i<right) i++;
        while(x<item[index[j]] && j>left) j--;

        if (i<=j) {
            if (item[index[i]] != item[index[j]]) {
                y=index[i];
                index[i]=index[j];
                index[j]=y;
            }
            i++; j--;
        }
    } while(i<=j);

    if (left<j) qs(item,index,left,j);
    if (i<right) qs(item,index,i,right);
}

template <class T>
void
quicksort(T* item, int* index, uli n)
{
    int i;
    if (n<=0) return;
    for (i=0; i<n; ++i) {
        index[i] = i;
    }
    qs<T>(item,index,0,n-1);
}

template <class T, class U>
void
Sort::accumulate(const T* src, U* dst, U scale) const
{
    if (scale == 0) //use the permutation sign
    {
        scale = p_->sign();
    }
    //sending tmp variable is necessary since value is passed by reference
    U* targetptr = dst;
    stride_accumulate<T,U>(0, src, targetptr, scale);
}

template <class T, class U>
void
Sort::stride_accumulate(usi sortlevel, const T* src, U*& target, U scale) const
{
      //not yet finished
    size_t src_stride = lengths_[sortlevel];
    size_t nstrides = nstrides_[sortlevel];

    ++sortlevel;
    if (sortlevel == nindex_) //finish off
    {
        const T* srcptr = src;
        U* dstptr = target;
        if (src_stride == 1)
        {
            for (size_t i=0; i < nstrides; ++srcptr, ++dstptr, ++i)
            {
                (*dstptr) += (*srcptr) * scale;
            }
        }
        else
        {
            for (size_t i=0; i < nstrides; srcptr += src_stride, ++dstptr, ++i)
            {
                (*dstptr) += (*srcptr) * scale;
            }
        }

        //increment the target pointer and return
        target += nstrides;
        return;
    }
    else
    {
        for (size_t i=0; i < nstrides; ++i, src += src_stride)
            stride_accumulate(sortlevel, src, target, scale);
    }
}

template <class T>
T*
Sort::sort(const T* src) const
{
    T* dst = new T[ntot_];
    sort<T>(src, dst);
    return dst;
}

template <class T, class U>
void
Sort::sort_noscale(const T* src, U* dst) const
{
    //sending tmp variable is necessary since value is passed by reference
    U* targetptr = dst;
    stride_sort(0, src, targetptr);
}

template <class T, class U>
void
Sort::sort(const T* src, U* dst) const
{
    sort_noscale<T,U>(src, dst);

    if (p_->sign() == -1) //scale all of the elements
    {
        U* dptr = dst;
        for (size_t i=0; i < ntot_; ++i, ++dptr)
            (*dptr) *= -1;
    }
}

template <class T, class U>
void
Sort::stride_sort(usi sortlevel, const T* src, U*& target) const
{
      //not yet finished
    size_t src_stride = lengths_[sortlevel];
    size_t nstrides = nstrides_[sortlevel];

    ++sortlevel;
    if (sortlevel == nindex_) //finish off
    {
        if (src_stride == 1)
        {
            ::memcpy(target, src, nstrides * sizeof(T));
        }
        else
        {
            const T* srcptr = src;
            U* dstptr = target;
            for (size_t i=0; i < nstrides; srcptr += src_stride, ++dstptr, ++i)
            {
                T entry = *srcptr;
                (*dstptr) = entry;
            }
        }

        //increment the target pointer and return
        target += nstrides;
        return;
    }
    else
    {
        for (size_t i=0; i < nstrides; ++i, src += src_stride)
            stride_sort(sortlevel, src, target);
    }
}

}

#endif

