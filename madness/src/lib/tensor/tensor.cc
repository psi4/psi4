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


  $Id: tensor.cc 2173 2011-02-23 21:40:46Z justus.c79@gmail.com $
*/


#define TENSOR_CC



#ifdef STATIC
#  undef STATIC
#endif

#if HAVE_UNQUALIFIED_STATIC_DECL
#  define STATIC static
#else
// Cray X1 compiler won't instantiate static function templates (mxm*)
#  define STATIC
#endif

#include <madness_config.h>
#include <misc/ran.h>
#include <tensor/tensor.h>
#include <tensor/mtxmq.h>
#include <tensor/aligned.h>

#include <stdlib.h>
#include <algorithm>
#include <complex>
#include <cmath>
#include <cstring>
#include <iostream>

namespace madness {
#include <tensor/mxm.h>
}

/// \file tensor.cc
/// \brief Completes the implementation of Tensor and instantiates all specializations for fast compiles.
/// \ingroup tensor
/// @{

namespace madness {

    template <class T> class Tensor;

    template <class T> class SliceTensor;

    /// Print a TensorException to the stream (for human consumption)
    std::ostream& operator <<(std::ostream& out, const TensorException& e) {
        out << "TensorException: msg='";
        if (e.msg) out << e.msg;
        out << "'\n";
        if (e.assertion) out << "                 failed assertion='" <<
            e.assertion << "'\n";
        out << "                 value=" << e.value << "\n";
        if (e.line) out << "                 line=" << e.line << "\n";
        if (e.function) out << "                 function='" <<
            e.function << "'\n";
        if (e.filename) out << "                 filename='" <<
            e.filename << "'\n";
        if (e.t) {
            out << "                 tensor=Tensor<";
            if (e.t->id()>=0 && e.t->id()<=TENSOR_MAX_TYPE_ID) {
                out << tensor_type_names[e.t->id()] << ">(";
            }
            else {
                out << "invalid_type_id>(";
            }
            if (e.t->ndim()>=0 && e.t->ndim()<TENSOR_MAXDIM) {
                for (int i=0; i<e.t->ndim(); ++i) {
                    out << e.t->dim(i);
                    if (i != (e.t->ndim()-1)) out << ",";
                }
                out << ")";
            }
            else {
                out << "invalid_dimensions)";
            }
            out << " at 0x" << (void *)(e.t) << "\n";
        }

        return out;
    }


    std::ostream& operator<<(std::ostream& stream, const Slice& s) {
        stream << "Slice(" << s.start << "," << s.end << "," << s.step << ")";
        return stream;
    }

    template<> float_complex Tensor<float_complex>::min(long* ind) const {
        TENSOR_EXCEPTION("cannot perform min on complex types",0,this);
        return 0;
    }
    template<> double_complex Tensor<double_complex>::min(long* ind) const {
        TENSOR_EXCEPTION("cannot perform min on complex types",0,this);
        return 0;
    }
    template<> float_complex Tensor<float_complex>::max(long* ind) const {
        TENSOR_EXCEPTION("cannot perform max on complex types",0,this);
        return 0;
    }
    template<> double_complex Tensor<double_complex>::max(long* ind) const {
        TENSOR_EXCEPTION("cannot perform max on complex types",0,this);
        return 0;
    }

    //#include "tensor_spec.h"
}
/// @}
