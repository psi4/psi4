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

  
  $Id: tensorexcept.h 1602 2009-12-27 19:53:06Z rjharrison $
*/

  
#ifndef MADNESS_TENSOR_TENSOREXCPT_H__INCLUDED
#define MADNESS_TENSOR_TENSOREXCPT_H__INCLUDED

/// \file tensorexcept.h
/// \brief Declares and implements TensorException

/// \ingroup tensor

#include <iosfwd>
#include <exception>

namespace madness {
/// Tensor is intended to throw only TensorExceptions
    class TensorException : public std::exception {
    public:
        const char* msg;
        const char* assertion;
        const int value;
        const BaseTensor* t;
        const int line;
        const char *function;
        const char *filename;

        // Capturing the line/function/filename info is best done with the macros below
        TensorException(const char* s, const char *a, int err, const BaseTensor* tp,
                        int lin, const char *func, const char *file)
                : msg(s)
                , assertion(a)
                , value(err)
                , t(tp)
                , line(lin)
                , function(func)
        , filename(file) {}

        virtual const char* what() const throw() {
            return msg;
        }
    };

// implemented in tensor.cc
    std::ostream& operator <<(std::ostream& out, const TensorException& e);

#define TENSOR_STRINGIZE(X) #X
#define TENSOR_EXCEPTION_AT(F, L) TENSOR_STRINGIZE(F) "(" TENSOR_STRINGIZE(L) ")"

#define TENSOR_EXCEPTION(msg,value,t) \
    throw ::madness::TensorException("TENSOR EXCPETION: " TENSOR_EXCEPTION_AT( __FILE__, __LINE__ ) ": " msg , \
    0,value,t,__LINE__,__FUNCTION__,__FILE__)

#define TENSOR_ASSERT(condition,msg,value,t) \
do {if (!(condition)) \
        throw ::madness::TensorException("TENSOR ASSERTION FAILED: " TENSOR_EXCEPTION_AT( __FILE__, __LINE__ ) ": " msg , \
        #condition,value,t,__LINE__,__FUNCTION__,__FILE__); \
   } while (0)

}

#endif // MADNESS_TENSOR_TENSOREXCPT_H__INCLUDED
