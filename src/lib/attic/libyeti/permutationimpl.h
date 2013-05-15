/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

#ifndef yeti_permutationimpl_h
#define yeti_permutationimpl_h

#include "permutation.h"

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

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

template <typename data_t>
bool
Permutation::fixes_arr(const data_t *indices) const
{
    for (usi i=0; i < nindex_; ++i)
    {
        if (indices[i] != indices[pmap_[i]])
        {
            return false;
        }
    }
    return true;
}

template <typename data_t>
PermutationGroupPtr
PermutationGroup::isotropy_grp(const data_t* indices)
{
    if (order() == 1) //this is just the identity group... send me back
        return this;

    usi nidx = nindex();
    PermutationGroupPtr grp = new PermutationGroup(nidx);
    PermutationGroup::iterator it(begin());
    PermutationGroup::iterator stop(end());
    for ( ; it != stop; ++it)
    {
        Permutation* p(*it);
        if (!p->is_identity() && p->fixes_arr<data_t>(indices))
            grp->add(p);
    }
    //automatically closed
    grp->closed_ = true;
    return grp;
}

}

#ifdef redefine_size_t
#undef size_t
#endif

#endif

