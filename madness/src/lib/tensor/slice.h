/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680


  $Id: slice.h 1602 2009-12-27 19:53:06Z rjharrison $
*/

/// \file slice.h
/// \brief Declares and implements Slice
/// \ingroup tensor

#ifndef MADNESS_TENSOR_SLICE_H__INCLUDED
#define MADNESS_TENSOR_SLICE_H__INCLUDED

#include <tensor/vector_factory.h>

namespace madness {

    /*!
      \brief A slice defines a sub-range or patch of a dimension.
      \ingroup tensor

      \par Introduction
     
      The slice \c Slice(start,end,step) defines the subvector
      \code
      [start, start+step, start+2*step, ..., end]
      \endcode
      with indices as if generated from these loops
      \code
      if (step > 0) 
          for (i=start; i<=end; i+=step) i;
      else if (step < 0) 
          for (i=start; i>=end; i+=step) i;
      else if (step == 0  && start == end) 
          defines a scalar
      else 
          runtime error detected when slice is used
      \endcode
      
      Note that \c start and \c end are \em inclusive, unlike the Python
      convention of specifying end+1 (note that we \em cannot do
      this easily in C/C++ unless we also define a special value
      to indicate the end of a dimension of unknown size).
      
      \par Special meanings and conventions

      Negative values for \c start or \c end (similar to Python)
      are relative to the end of the (possibly unknown) dimension.
      E.g.,
      - \c end=-1  is equivalent to    \c end=dim-1
      - \c start=-4 is equivalent to \c start=dim-4
      -  \c step=0 and \c start==end implies dimension will be eliminated
      when the slice is used to index a tensor
      -  The default constructor specifies the entire dimension
      -  If the input length is not an exact multiple of step, end is
      rounded towards start to recover the behaviour of the
      \c <= or \c >= bounds in the loops specified above.
      
      Special slices have been defined as constants
      - \c _ (1 underscore) = \c Slice(0,-1,1) = full slice in current order
      - \c ___ (3 underscores) = Array of Slices with value \c _ so that \c t(___)
      will generate an assignable view of the entire tensor \c t .
      - _reverse = \c Slice(-1,0,-1) = full dimension reversed

      \par Examples

      - \c Slice() --- full slice in current order
      - \c Slice(0,-1,1) --- full slice in current order
      - \c Slice(3,3,0) --- eliminate this dimension setting index=3 (step=0)
      - \c Slice(3,3,1) --- reduce a dimension to length 1 using index=3 (step=1)
      - \c Slice(-1,0,-1) --- reverse a dimension
      - \c Slice(0,-1,2) --- use all even numbered elements
      - \c Slice(1,-1,2) --- use all odd numbered elements
    */
    class Slice {
        
    public:
        long start;
        long end;
        long step;
        
        inline Slice() : start(0), end(-1), step(1) {};
        inline Slice(long s, long e, long stp=1) : start(s), end(e), step(stp) {};
        inline Slice& operator=(const Slice& s) {
            start=s.start;
            end=s.end;
            step=s.step;
            return *this;
        };
    };

    std::ostream& operator<<(std::ostream& stream, const Slice& s);

    static const Slice _(0,-1,1);	/// Entire dimension
    static const std::vector<Slice> ___ = vector_factory(_,_,_,_,_,_); /// Entire tensor
    static const Slice _reverse(-1,0,-1); /// Reversed dimension

}
#endif // MADNESS_TENSOR_SLICE_H__INCLUDED
