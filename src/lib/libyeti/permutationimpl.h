#ifndef yeti_permutationimpl_h
#define yeti_permutationimpl_h

#include "permutation.h"

namespace yeti {

template <class T>
T*
Permutation::permute(const T* indices) const
{
    T* image = new T[nindex_];
    permute<T>(indices, image);
    return image;
}

template <class T>
void
Permutation::permute(const T* src, T* dest) const
{
    const usi* pptr = pmap_;
    T* imgptr = dest;
    for (usi i = 0; i < nindex_; ++i, ++imgptr, ++pptr)
    {
       (*imgptr) = src[*pptr];
    }
}

}

#endif

