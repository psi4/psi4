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


  $Id: transform3d.cc 2173 2011-02-23 21:40:46Z justus.c79@gmail.com $
*/


// #include <iostream>
// #include <tensor/tensor.h>

// namespace {
//   template <typename T>
//   inline T** create_2d_array(unsigned int d1, unsigned int d2);
//   template <typename T>
//   inline void delete_2d_array(T** A);
// }
//
//namespace madness {

/// optomized 3d transform
/// c must be 2d tensor(matrix) and must square
/// t must be 3d tensor and must be a cube
/// t and c dims must match
/// t and c must be contiguous
///
/// \code
/// result(i,j,k) <-- sum(i',j',k') A(i',j',k') C(i',i) C(j',j) C(k',k)
/// \endcode
///
template <class T>
Tensor<T> transform3d(const Tensor<T>& t, const Tensor<T>& c) {
    TENSOR_ASSERT(c.ndim == 2,"c must be 2d tensor(matrix)",c.ndim,&c);
    if ((t.ndim != 3) ||
            (t.dim[0] != t.dim[1]) ||
            (t.dim[0] != t.dim[2]) ||
            (c.dim[0] != t.dim[0]) ||
            (c.dim[0] != c.dim[1]) ||
            (!t.iscontiguous()) ||
            (!c.iscontiguous()))
        return transform(t,c);

    long d0 = c.dim[0];
    long d0_squared = d0*d0;
    long d0_cubed = d0_squared*d0;

    // Tensor constructor returns a zero filled tensor
    Tensor<T> result = Tensor<T>(d0,d0,d0);
#ifdef IBMXLC
    TENSOR_ASSERT(c.dim[0] < 36,"hard dimension failure",c.dim[0],&c);
#endif
    T* tmp = new T[d0_cubed];

    T* restrict r_p = result.ptr();
    T* restrict t_p = t.ptr();
    T* restrict c_p = c.ptr();
    T* restrict tmp_p = tmp;

    for (long i=0; i<d0_cubed; ++i)
        tmp[i] = (T)0;

    // both result and tmp are zero filled at this point
    // Transform along 1st dimension
    // result gets "result"
    mTxm(d0_squared, d0, d0, r_p, t_p, c_p);

    // Transform along 2nd dimension
    // tmp gets "result"
    mTxm(d0_squared, d0, d0, tmp_p, r_p, c_p);

    // Transform along 3rd dimension
    result.fill(0);
    // result gets "result"
    mTxm(d0_squared, d0, d0, r_p, tmp_p, c_p);

    delete[] tmp;
    return result;
}


Tensor<double_complex> transform3d(const Tensor<double_complex>& t, const Tensor<double>& c) {
    TENSOR_ASSERT(c.ndim == 2,"c must be 2d tensor(matrix)",c.ndim,&c);
    if ((t.ndim != 3) ||
            (t.dim[0] != t.dim[1]) ||
            (t.dim[0] != t.dim[2]) ||
            (c.dim[0] != t.dim[0]) ||
            (c.dim[0] != c.dim[1]) ||
            (!t.iscontiguous()) ||
            (!c.iscontiguous()))
        TENSOR_EXCEPTION("transform3d:double_complex*double",0,&t);

    Tensor<double_complex> result = Tensor<double_complex>(c.dim[1],c.dim[1],c.dim[1]);
#if 0
    double_complex v_jpkpi[c.dim[0]][c.dim[0]];
    double_complex v_kpij[c.dim[0]];
#endif
    //double_complex** v_jpkpi = ::create_2d_array<double_complex>(c.dim[0],c.dim[0]);
    TENSOR_ASSERT(c.dim[0] < 36,"hard dimension failure",c.dim[0],&c);
    double_complex v_kpjpi[36][36];
    double_complex* v_kpij = new double_complex[c.dim[0]];
    double_complex v_ijk;
    long c_d0 = c.dim[0];
    long c_d1 = c.dim[1];

    // Transform along 1st dimension
    for (long i = 0; i < c_d1; ++i) {
        for (long kp = 0; kp < c_d0; ++kp) {
            for (long jp = 0; jp < c_d0; ++jp) {
                double_complex sum = 0.0;
                for (long ip = 0; ip < c_d0; ++ip)  sum += t(ip,jp,kp) * c(ip,i);
                v_kpjpi[kp][jp] = sum;
            }
        }

        // Transform along 2nd dimension
        for (long j = 0; j < c_d1; ++j) {
            for (long kp = 0; kp < c_d0; ++kp) {
                v_kpij[kp] = (double_complex)0;
                for (long jp = 0; jp < c_d0; ++jp)
                    v_kpij[kp] +=  v_kpjpi[kp][jp] * c(jp,j);
            }

            // Transform along 3rd dimension
            for (long k = 0; k < c_d1; ++k) {
                v_ijk = (double_complex)0;
                for (long kp = 0; kp < c_d0; ++kp)
                    v_ijk +=  v_kpij[kp] * c(kp,k);

                result(i,j,k) = v_ijk;
            }
        }
    }

    delete[] v_kpij;
    return result;
}

#if 0

// this seems to perform best on my G4

/// optomized 3d transform
/// c must be 2d tensor(matrix) and must square
/// t must be 3d tensor and must be a cube
/// t and c dims must match
/// t and c must be contiguous
///
/// \code
/// result(i,j,k) <-- sum(i',j',k') A(i',j',k') C(i',i) C(j',j) C(k',k)
/// \endcode
///
template <class T>
Tensor<T> transform3d(const Tensor<T>& t, const Tensor<T>& c) {
    TENSOR_ASSERT(c.ndim == 2,"c must be 2d tensor(matrix)",c.ndim,&c);
    if ((t.ndim != 3) ||
            (t.dim[0] != t.dim[1]) ||
            (t.dim[0] != t.dim[2]) ||
            (c.dim[0] != t.dim[0]) ||
            (c.dim[0] != c.dim[1]) ||
            (!t.iscontiguous()) ||
            (!c.iscontiguous()))
        return transform(t,c);

    Tensor<T> result = Tensor<T>(c.dim[1],c.dim[1],c.dim[1]);
    T v_jpkpi[c.dim[0]][c.dim[0]];
    T v_kpij[c.dim[0]];
    T v_ijk;
    long c_d0 = c.dim[0];
    long c_d1 = c.dim[1];

    // Transform along 1st dimension
    for (long i = 0; i < c_d1; ++i) {
        for (long kp = 0; kp < c_d0; ++kp) {
            for (long jp = 0; jp < c_d0; ++jp) {
                v_jpkpi[jp][kp] = (T)0;
                for (long ip = 0; ip < c_d0; ++ip)
                    v_jpkpi[jp][kp] += t(ip,jp,kp) * c(ip,i);
            }
        }

        // Transform along 2nd dimension
        for (long j = 0; j < c_d1; ++j) {
            for (long kp = 0; kp < c_d0; ++kp) {
                v_kpij[kp] = (T)0;
                for (long jp = 0; jp < c_d0; ++jp)
                    v_kpij[kp] +=  v_jpkpi[jp][kp] * c(jp,j);
            }

            // Transform along 3rd dimension
            for (long k = 0; k < c_d1; ++k) {
                v_ijk = (T)0;
                for (long kp = 0; kp < c_d0; ++kp)
                    v_ijk +=  v_kpij[kp] * c(kp,k);

                result(i,j,k) = v_ijk;
            }
        }
    }
    return result;
}
#endif

/// optomized 3d transform with three different C's
/// all three c must be 2d tensor(matrix) and must square
/// all three c dims must match
/// t must be 3d tensor and must be a cube
/// t and c dims must match
/// t and c must be contiguous
///
/// \code
/// result(i,j,k) <-- sum(i',j',k') A(i',j',k') C0(i',i) C1(j',j) C2(k',k)
/// \endcode
///
template <class T>
Tensor<T> transform3d_3c(const Tensor<T>& t,
                         const Tensor<T>& c0, const Tensor<T>& c1, const Tensor<T>& c2) {
    TENSOR_ASSERT(c0.ndim == 2,"c must be 2d tensor(matrix)",c0.ndim,&c0);
    TENSOR_ASSERT(t.ndim == 3,"t must be 3d tensor",t.ndim,&t);
    TENSOR_ASSERT(((t.dim[0] == t.dim[1]) && (t.dim[0] == t.dim[2])),"t must be a cube",t.ndim,&t);
    TENSOR_ASSERT(c0.dim[0] == t.dim[0],"t and c dims must match",t.ndim,&t);
    TENSOR_ASSERT(c0.dim[0] == c0.dim[1],"c0 must square",0,&c0);
    TENSOR_ASSERT(c1.dim[0] == c1.dim[1],"c1 must square",0,&c1);
    TENSOR_ASSERT(c2.dim[0] == c2.dim[1],"c2 must square",0,&c2);
    TENSOR_ASSERT((c0.dim[0] == c1.dim[0]) && (c0.dim[0] == c2.dim[0]),"c's dims must be equal",0,0);
    TENSOR_ASSERT((t.iscontiguous() && c0.iscontiguous() && c1.iscontiguous() && c2.iscontiguous()),"t and c must be contiguous",t.ndim,&t);

    long d0 = c0.dim[0];
    long d0_squared = d0*d0;
    long d0_cubed = d0_squared*d0;

    // Tensor constructor returns a zero filled tensor
    Tensor<T> result = Tensor<T>(d0,d0,d0);
#ifdef IBMXLC
    TENSOR_ASSERT(c0.dim[0] < 36,"hard dimension failure",c0.dim[0],&c0);
#endif
    T* tmp = new T[d0_cubed];

    T* restrict r_p = result.ptr();
    T* restrict t_p = t.ptr();
    T* restrict c0_p = c0.ptr();
    T* restrict c1_p = c1.ptr();
    T* restrict c2_p = c2.ptr();
    T* restrict tmp_p = tmp;

    for (long i=0; i<d0_cubed; ++i)
        tmp[i] = (T)0;

    // both result and tmp are zero filled at this point
    // Transform along 1st dimension with C0
    // result gets "result"
    mTxm(d0_squared, d0, d0, r_p, t_p, c0_p);

    // Transform along 2nd dimension with C1
    // tmp gets "result"
    mTxm(d0_squared, d0, d0, tmp_p, r_p, c1_p);

    // Transform along 3rd dimension with C2
    result.fill(0);
    // result gets "result"
    mTxm(d0_squared, d0, d0, r_p, tmp_p, c2_p);

    delete[] tmp;
    return result;
}


template Tensor<double> transform3d(const Tensor<double>& t, const Tensor<double>& c);

template Tensor<double> transform3d_3c(const Tensor<double>& t,
                                       const Tensor<double>& c0,
                                       const Tensor<double>& c1,
                                       const Tensor<double>& c2);


//}

/*
  namespace {

  template <typename T>
  inline T** create_2d_array(unsigned int d1, unsigned int d2) {
  if (d1 == 0 || d2 == 0) return 0;
  T** result = new T*[d1];
  result[0] = new T[d1*d2];
  for(unsigned int i=1; i<d1; ++i) result[i] = result[i-1] + d2;

  return result;
  }

  template <typename T>
  inline void delete_2d_array(T** A) {
  delete[] A[0];
  delete[] A;
  }
  }
*/
