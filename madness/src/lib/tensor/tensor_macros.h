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


  $Id: tensor_macros.h 2173 2011-02-23 21:40:46Z justus.c79@gmail.com $
*/


#ifndef MADNESS_TENSOR_TENSOR_MACROS_H__INCLUDED
#define MADNESS_TENSOR_TENSOR_MACROS_H__INCLUDED

/*!
  \file tensor_macros.h
  \brief Macros for easy and efficient iteration over tensors.

  \ingroup tensor

Several different macros have been defined to make it
easy to iterate over expressions involving tensors.  They
vary in their generality, ease of use, and efficiency.

The most general, most easy to use, but also most inefficient,
and least safe, is
\code
ITERATOR(t, expression)
\endcode
where \c t is a Tensor of any type, size or dimension that is used to
define the range of the loop indices, and expression can be nearly
anything, including multi-line expressions performing arbitrary
operations on multiple tensors.  The loop indices, going
from left to right in the dimensions, are
\code
_i, _j, _k, ...
\endcode

E.g., to add two matrices together (there are more efficient ways
to do this, such as \c a+=b )
\code
Tensor<long> a(4,2), b(4,2);
ITERATOR(a, a(_i,_j) += b(_i,_j));
\endcode

E.g., to print out the indices of all elements of a matrix
greater than 0.5;
\code
Tensor<float> m(5,5);
m.fillrandom();
ITERATOR(m, if (m(_i,_j) > 0.5) {
               cout << _i << " " << _j << endl;
            });
\endcode

To make it possible to index arbitrary dimension tensors, the macro
\c IND has been defined as the indices for the highest supported
dimension.  E.g., to elementwise divide the contents of two tensors of
unknown dimension
\code
ITERATOR(x, x(IND)/y(IND));
\endcode

Note that using \c IND employs bounds checking where as direct indexing
with \c _i , etc., does not.

The generality of these macros is offset by their inefficiency and
lack of safety.  The inefficiency is twofold.  First, the \c ITERATOR
macro generates a separate block of code for each possible dimension.
This could cause code bloat and increased compilation time.  To
solve this problem, the macros \c ITERATOR1 , \c ITERATOR2, etc., have
been defined, with the corresponding \c IND1 , \c IND2 , etc.  These
macros may be applied to tensor expressions of the appropriate
dimension.

The second inefficiency is at runtime, due to the explicit indexing of
all the tensor expressions and the inability to optimize the order
in which memory is traversed.  The lack of safety is the inability
to check that the tensors in the expression conform and that the
indices are not out of bounds.

The safety and cost of explicit indexing are addressed by the macros
\c UNARYITERATOR , \c BINARYITERATOR , and \c TERNARYITERATOR , along
with their specialization to specific numbers of dimensions (again by
appending the dimension number to the name of the macro).  These
macros are safer since you have to explicitly name the tensors you are
iterating over, so that the macro can now check that the input tensors
conform.  The cost of looping is reduced by replacing explicit
indexing with pointer arithmetic.  These macros still define the loop
indices \c _i , \c _j , etc., but also define \c _p0 , \c _p1 , etc.,
as pointers to the current elements of tensor argument 0, tensor
argument 1, etc..

E.g., set elements of a 3-d tensor, \c t , of type \c double to a
function of the indices
\code
UNARYITERATOR(double, t, *_p0 = 1.0/(_i + _j + 1.0));
\endcode

E.g., to merge two \c double tensors as real and imaginary parts
of complex tensor of any dimension
\code
TERNARYITERATOR(double_complex, c, double, r, double, i,
.               *_p0 = double_complex(*_p1, *_p2));
\endcode

However, we still have the problems that if the dimensions of a tensor
have been reordered, the loops will go through memory inefficiently,
and the dimension independent macros still generate redundant code
blocks.  Also, the innermost loop might not be very long and will
be inefficient.

The most general, efficient and code-compact macros internally use the
\c TensorIterator , which you could also use directly.  Since there is
no nest of explicit loops, the tensor indices are no longer available
as \c _i , \c _j , etc..  Furthermore, the \c TensorIterator can
reorder the loops to optimize the memory traversal, and fuse
dimensions to make the innermost loop longer for better vectorization
and reduced loop overhead.

The most efficient macros for iteration are \c UNARY_OPTIMIZED_ITERATOR ,
\c BINARY_OPTIMIZED_ITERATOR , and \c TERNARY_OPTIMIZED_ITERATOR .
As before, these define the pointers \c _p0 , \c _p1, \c _p2 , which
point to the current (and corresponding) element of each argument
tensor.  However, unlike the previous macros there is no guarantee
that the elements are looped thru in the order expected by a
simple nest of loops.  Furthermore, the indices are completely
unvailable.  In addition to using the iterators for optimal
traversal, these macros attempt to use a single loop
for optimal vector performance.

E.g., the most efficient and safe way to perform the previous example
of merging two \c double tensors as real and imaginary parts
of a complex tensor of any dimension
\code
TERNARY_OPTIMIZED_ITERATOR(double_complex, c, double, r, double, i,
           .               *_p0 = double_complex(*_p1, *_p2));
\endcode
This is precisely how most internal operations are implemented.

In some situations it is necessary to preserve the expected
order of loops and to not fuse dimensions.  The macros
\c UNARY_UNOPTIMIZED_ITERATOR , \c BINARY_UNOPTIMIZED_ITERATOR ,
and \c TERNARY_UNOPTIMIZED_ITERATOR use the \c TensorIterator
but disable loop reordering and fusing.  Once these optimizations
have been turned off, the loop indices are avaiable, if needed,
from the \c ind[] member of the iterator (which is named
\c _iter ).

E.g., the fillindex() method is implemented as follows
\code
long count = 0;
UNARY_UNOPTIMIZED_ITERATOR(T, (*this), *_p0 = (T) count++);
\endcode

\em NB: None of the above iterator macros can be nested ... use the
actual tensor iterator to do this.

Recommendation --- for both efficiency and safety, use the optimized
macros (\c UNARY_OPTMIZED_ITERATOR , etc.), unless it is necessary to
preserve loop order, in which case use the unoptimized versions.  If
you need the loop indices, use the macros \c UNARY_ITERATOR, etc.,
unless you have a very general expression that they cannot handle.  In
this last instance, or for ease of rapid implementation, use the general
\c ITERATOR macro first described.

*/

// don't change this without changing the iterator macros
#define TENSOR_MAXDIM 6

/// Macros IND1, ..., IND6, and IND are a convenience for indexing in macro iterators.

#define IND1 _i
#define IND2 _i,_j
#define IND3 _i,_j,_k
#define IND4 _i,_j,_k,_l
#define IND5 _i,_j,_k,_l,_m
#define IND6 _i,_j,_k,_l,_m,_n
#define IND  IND6


#define ITERATOR1(t,exp) do { \
        long __xd0=t.dim(0),_index=0;                                   \
        for (long _i=0; _i<__xd0; ++_i) {exp;_index++;} } while (0)

#define ITERATOR2(t,exp) do { \
        long __xd0=t.dim(0), __xd1=t.dim(1), _index=0;  \
for (long _i=0; _i<__xd0; ++_i) { \
  for (long _j=0; _j<__xd1; ++_j) {exp;_index++;} } } while (0)

#define ITERATOR3(t,exp) do { \
        long __xd0=t.dim(0), __xd1=t.dim(1), __xd2=t.dim(2), _index=0;  \
for (long _i=0; _i<__xd0; ++_i) { \
  for (long _j=0; _j<__xd1; ++_j) { \
    for (long _k=0; _k<__xd2; ++_k) {exp;_index++;} } } } while (0)

#define ITERATOR4(t,exp) do { \
        long __xd0=t.dim(0), __xd1=t.dim(1), __xd2=t.dim(2),        \
            __xd3=t.dim(3), _index=0;                               \
for (long _i=0; _i<__xd0; ++_i) { \
  for (long _j=0; _j<__xd1; ++_j) { \
    for (long _k=0; _k<__xd2; ++_k) { \
      for (long _l=0; _l<__xd3; ++_l) {exp;_index++;} } } } } while (0)

#define ITERATOR5(t,exp) do { \
        long __xd0=t.dim(0), __xd1=t.dim(1), __xd2=t.dim(2),      \
            __xd3=t.dim(3), __xd4=t.dim(4), _index=0;             \
for (long _i=0; _i<__xd0; ++_i) { \
  for (long _j=0; _j<__xd1; ++_j) { \
    for (long _k=0; _k<__xd2; ++_k) { \
      for (long _l=0; _l<__xd3; ++_l) { \
        for (long _m=0; _m<__xd4; ++_m) {exp;_index++;} } } } } } while (0)

#define ITERATOR6(t,exp) do { \
        long __xd0=t.dim(0), __xd1=t.dim(1), __xd2=t.dim(2),            \
            __xd3=t.dim(3), __xd4=t.dim(4), __xd5=t.dim(5), _index=0;;  \
for (long _i=0; _i<__xd0; ++_i) { \
  for (long _j=0; _j<__xd1; ++_j) { \
    for (long _k=0; _k<__xd2; ++_k) { \
      for (long _l=0; _l<__xd3; ++_l) { \
        for (long _m=0; _m<__xd4; ++_m) { \
          for (long _n=0; _n<__xd5; ++_n) {exp;_index++;} } } } } } } while(0)

#define ITERATOR(t,exp) do { \
  long _j=0, _k=0, _l=0, _m=0, _n=0; \
  if (t.ndim() == 1) {ITERATOR1(t,exp);} \
  else if (t.ndim() == 2) {ITERATOR2(t,exp);} \
  else if (t.ndim() == 3) {ITERATOR3(t,exp);} \
  else if (t.ndim() == 4) {ITERATOR4(t,exp);} \
  else if (t.ndim() == 5) {ITERATOR5(t,exp);} \
  else if (t.ndim() == 6) {ITERATOR6(t,exp);} \
  else {TENSOR_ASSERT(t.ndim() <= 6,"ndim confused?",t.ndim(),&t);} \
 } while(0)

// Inside iterator access pointer to current element as _p0 (pointer to
// argument number 0).  _i, _j, _k, ..., also defined
#define UNARYITERATOR1(X,x,exp) do { \
        long __xd0=x.dim(0);         \
        long __xs0=x.stride(0);      \
X* restrict _p0=x.ptr(); \
for (long _i=0; _i<__xd0; ++_i,_p0+=__xs0) { \
  exp; \
} } while(0)

#define UNARYITERATOR2(X,x,exp) do { \
        long __xd0=x.dim(0), __xd1=x.dim(1);   \
        long __xs0=x.stride(0), __xs1=x.stride(1);      \
X* restrict __xp0=x.ptr(); \
for (long _i=0; _i<__xd0; ++_i,__xp0+=__xs0) { \
  X* restrict _p0=__xp0; \
  for (long _j=0; _j<__xd1; ++_j, _p0+=__xs1) { \
    exp; \
  } } } while(0)

#define UNARYITERATOR3(X,x,exp) do { \
        long __xd0=x.dim(0), __xd1=x.dim(1), __xd2=x.dim(2);        \
        long __xs0=x.stride(0), __xs1=x.stride(1), __xs2=x.stride(2);   \
X* restrict __xp0=x.ptr(); \
for (long _i=0; _i<__xd0; ++_i,__xp0+=__xs0) { \
  X* restrict __xp1=__xp0; \
  for (long _j=0; _j<__xd1; ++_j, __xp1+=__xs1) { \
    X* restrict _p0=__xp1; \
    for (long _k=0; _k<__xd2; ++_k, _p0+=__xs2) { \
       exp; \
    } } } } while(0)

#define UNARYITERATOR4(X,x,exp) do { \
        long __xd0=x.dim(0), __xd1=x.dim(1), __xd2=x.dim(2),    \
            __xd3=x.dim(3);                                         \
        long __xs0=x.stride(0), __xs1=x.stride(1), __xs2=x.stride(2),   \
            __xs3=x.stride(3);                                          \
X* restrict __xp0=x.ptr(); \
for (long _i=0; _i<__xd0; ++_i,__xp0+=__xs0) { \
  X* restrict __xp1=__xp0; \
  for (long _j=0; _j<__xd1; ++_j, __xp1+=__xs1) { \
    X* restrict __xp2=__xp1; \
    for (long _k=0; _k<__xd2; ++_k, __xp2+=__xs2) { \
       X* restrict _p0=__xp2; \
       for (long _l=0; _l<__xd3; ++_l, _p0+=__xs3) { \
          exp; \
       } } } } } while(0)

#define UNARYITERATOR5(X,x,exp) do { \
        long __xd0=x.dim(0), __xd1=x.dim(1), __xd2=x.dim(2),    \
            __xd3=x.dim(3), __xd4=x.dim(4);                         \
        long __xs0=x.stride(0), __xs1=x.stride(1), __xs2=x.stride(2),   \
            __xs3=x.stride(3), __xs4=x.stride(4);                       \
X* restrict __xp0=x.ptr(); \
for (long _i=0; _i<__xd0; ++_i,__xp0+=__xs0) { \
  X* restrict __xp1=__xp0; \
  for (long _j=0; _j<__xd1; ++_j, __xp1+=__xs1) { \
    X* restrict __xp2=__xp1; \
    for (long _k=0; _k<__xd2; ++_k, __xp2+=__xs2) { \
       X* restrict __xp3=__xp2; \
       for (long _l=0; _l<__xd3; ++_l, __xp3+=__xs3) { \
          X* restrict _p0 =__xp3; \
          for (long _m=0; _m<__xd4; ++_m, _p0+=__xs4) { \
            exp; \
          } } } } } } while(0)

#define UNARYITERATOR6(X,x,exp) do { \
        long __xd0=x.dim(0), __xd1=x.dim(1), __xd2=x.dim(2),    \
            __xd3=x.dim(3), __xd4=x.dim(4), __xd5=x.dim(5);         \
        long __xs0=x.stride(0), __xs1=x.stride(1), __xs2=x.stride(2),   \
            __xs3=x.stride(3), __xs4=x.stride(4), __xs5=x.stride(5);    \
X* restrict __xp0=x.ptr(); \
for (long _i=0; _i<__xd0; ++_i,__xp0+=__xs0) { \
  X* restrict __xp1=__xp0; \
  for (long _j=0; _j<__xd1; ++_j, __xp1+=__xs1) { \
    X* restrict __xp2=__xp1; \
    for (long _k=0; _k<__xd2; ++_k, __xp2+=__xs2) { \
       X* restrict __xp3=__xp2; \
       for (long _l=0; _l<__xd3; ++_l, __xp3+=__xs3) { \
          X* restrict __xp4=__xp3; \
          for (long _m=0; _m<__xd4; ++_m, __xp4+=__xs4) { \
            X* restrict _p0=__xp4; \
            for (long _n=0; _n<__xd5; ++_n, _p0+=__xs5) { \
              exp; \
          } } } } } } } while(0)

#define UNARYITERATOR(X,x,exp) do { \
  long _j=0, _k=0, _l=0, _m=0, _n=0; \
  if (x.ndim() == 1) UNARYITERATOR1(X,x,exp); \
  else if (x.ndim() == 2) UNARYITERATOR2(X,x,exp); \
  else if (x.ndim() == 3) UNARYITERATOR3(X,x,exp); \
  else if (x.ndim() == 4) UNARYITERATOR4(X,x,exp); \
  else if (x.ndim() == 5) UNARYITERATOR5(X,x,exp); \
  else if (x.ndim() == 6) UNARYITERATOR6(X,x,exp); \
  else {TENSOR_ASSERT(x.ndim() <= 6,"ndim confused?",x.ndim(),&x);} } while(0)

// Inside iterator access pointers to current elements as _p0 & _p1
// _i, _j, _k, ... also defined
#define BINARYITERATOR1(X,x,Y,y,exp) do { \
TENSOR_ASSERT(x.conforms(y),"first and second tensors do not conform",0,&x); \
 long __xd0=x.dim(0);                                                   \
 long __xs0=x.stride(0);                                                \
 long __ys0=y.stride(0);                                                \
X* restrict _p0=x.ptr(); \
Y* restrict _p1=y.ptr(); \
for (long _i=0; _i<__xd0; ++_i, _p0+=__xs0, _p1+=__ys0) { \
   exp; \
} } while(0)

#define BINARYITERATOR2(X,x,Y,y,exp) do { \
TENSOR_ASSERT(x.conforms(y),"first and second tensors do not conform",0,&x); \
 long __xd0=x.dim(0), __xd1=x.dim(1);                                   \
 long __xs0=x.stride(0), __xs1=x.stride(1);                             \
 long __ys0=y.stride(0), __ys1=y.stride(1);                             \
X* restrict __xp0=x.ptr(); \
Y* restrict __yp0=y.ptr(); \
for (long _i=0; _i<__xd0; ++_i, __xp0+=__xs0, __yp0+=__ys0) { \
  X* restrict _p0=__xp0; \
  Y* restrict _p1=__yp0; \
  for (long _j=0; _j<__xd1; ++_j, _p0+=__xs1, _p1+=__ys1) { \
     exp; \
  } } } while(0)

#define BINARYITERATOR3(X,x,Y,y,exp) do { \
TENSOR_ASSERT(x.conforms(y),"first and second tensors do not conform",0,&x); \
 long __xd0=x.dim(0), __xd1=x.dim(1), __xd2=x.dim(2);                   \
 long __xs0=x.stride(0), __xs1=x.stride(1), __xs2=x.stride(2);          \
 long __ys0=y.stride(0), __ys1=y.stride(1), __ys2=y.stride(2);          \
X* restrict __xp0=x.ptr(); \
Y* restrict __yp0=y.ptr(); \
for (long _i=0; _i<__xd0; ++_i, __xp0+=__xs0, __yp0+=__ys0) { \
  X* restrict __xp1=__xp0; \
  Y* restrict __yp1=__yp0; \
  for (long _j=0; _j<__xd1; ++_j, __xp1+=__xs1, __yp1+=__ys1) { \
    X* restrict _p0=__xp1; \
    Y* restrict _p1=__yp1; \
    for (long _k=0; _k<__xd2; ++_k, _p0+=__xs2, _p1+=__ys2) { \
       exp; \
    } } } } while(0)

#define BINARYITERATOR4(X,x,Y,y,exp) do { \
TENSOR_ASSERT(x.conforms(y),"first and second tensors do not conform",0,&x); \
 long __xd0=x.dim(0), __xd1=x.dim(1), __xd2=x.dim(2),                   \
     __xd3=x.dim(3);                                                    \
 long __xs0=x.stride(0), __xs1=x.stride(1), __xs2=x.stride(2),          \
     __xs3=x.stride(3);                                                 \
 long __ys0=y.stride(0), __ys1=y.stride(1), __ys2=y.stride(2),          \
     __ys3=y.stride(3);                                                 \
X* restrict __xp0=x.ptr(); \
Y* restrict __yp0=y.ptr(); \
for (long _i=0; _i<__xd0; ++_i, __xp0+=__xs0, __yp0+=__ys0) { \
  X* restrict __xp1=__xp0; \
  Y* restrict __yp1=__yp0; \
  for (long _j=0; _j<__xd1; ++_j, __xp1+=__xs1, __yp1+=__ys1) { \
    X* restrict __xp2=__xp1; \
    Y* restrict __yp2=__yp1; \
    for (long _k=0; _k<__xd2; ++_k, __xp2+=__xs2, __yp2+=__ys2) { \
      X* restrict _p0=__xp2; \
      Y* restrict _p1=__yp2; \
      for (long _l=0; _l<__xd3; ++_l, _p0+=__xs3, _p1+=__ys3) { \
         exp; \
      } } } } } while(0)

#define BINARYITERATOR5(X,x,Y,y,exp) do { \
TENSOR_ASSERT(x.conforms(y),"first and second tensors do not conform",0,&x); \
 long __xd0=x.dim(0), __xd1=x.dim(1), __xd2=x.dim(2),                   \
     __xd3=x.dim(3), __xd4=x.dim(4);                                    \
 long __xs0=x.stride(0), __xs1=x.stride(1), __xs2=x.stride(2),          \
     __xs3=x.stride(3), __xs4=x.stride(4);                              \
 long __ys0=y.stride(0), __ys1=y.stride(1), __ys2=y.stride(2),          \
     __ys3=y.stride(3), __ys4=y.stride(4);                              \
X* restrict __xp0=x.ptr(); \
Y* restrict __yp0=y.ptr(); \
for (long _i=0; _i<__xd0; ++_i, __xp0+=__xs0, __yp0+=__ys0) { \
  X* restrict __xp1=__xp0; \
  Y* restrict __yp1=__yp0; \
  for (long _j=0; _j<__xd1; ++_j, __xp1+=__xs1, __yp1+=__ys1) { \
    X* restrict __xp2=__xp1; \
    Y* restrict __yp2=__yp1; \
    for (long _k=0; _k<__xd2; ++_k, __xp2+=__xs2, __yp2+=__ys2) { \
      X* restrict __xp3=__xp2; \
      Y* restrict __yp3=__yp2; \
      for (long _l=0; _l<__xd3; ++_l, __xp3+=__xs3, __yp3+=__ys3) { \
        X* restrict _p0=__xp3; \
        Y* restrict _p1=__yp3; \
        for (long _m=0; _m<__xd4; ++_m, _p0+=__xs4, _p1+=__ys4) { \
           exp; \
        } } } } } } while(0)

#define BINARYITERATOR6(X,x,Y,y,exp) do { \
TENSOR_ASSERT(x.conforms(y),"first and second tensors do not conform",0,&x); \
 long __xd0=x.dim(0), __xd1=x.dim(1), __xd2=x.dim(2),                   \
     __xd3=x.dim(3), __xd4=x.dim(4), __xd5=x.dim(5);                    \
 long __xs0=x.stride(0), __xs1=x.stride(1), __xs2=x.stride(2),          \
     __xs3=x.stride(3), __xs4=x.stride(4), __xs5=x.stride(5);           \
 long __ys0=y.stride(0), __ys1=y.stride(1), __ys2=y.stride(2),          \
     __ys3=y.stride(3), __ys4=y.stride(4), __ys5=y.stride(5);           \
X* restrict __xp0=x.ptr(); \
Y* restrict __yp0=y.ptr(); \
for (long _i=0; _i<__xd0; ++_i, __xp0+=__xs0, __yp0+=__ys0) { \
  X* restrict __xp1=__xp0; \
  Y* restrict __yp1=__yp0; \
  for (long _j=0; _j<__xd1; ++_j, __xp1+=__xs1, __yp1+=__ys1) { \
    X* restrict __xp2=__xp1; \
    Y* restrict __yp2=__yp1; \
    for (long _k=0; _k<__xd2; ++_k, __xp2+=__xs2, __yp2+=__ys2) { \
      X* restrict __xp3=__xp2; \
      Y* restrict __yp3=__yp2; \
      for (long _l=0; _l<__xd3; ++_l, __xp3+=__xs3, __yp3+=__ys3) { \
        X* restrict __xp4=__xp3; \
        Y* restrict __yp4=__yp3; \
        for (long _m=0; _m<__xd4; ++_m, __xp4+=__xs4, __yp4+=__ys4) { \
          X* restrict _p0=__xp4; \
          Y* restrict _p1=__yp4; \
          for (long _n=0; _n<__xd5; ++_n, _p0+=__xs5, _p1+=__ys5) { \
             exp; \
          } } } } } } } while(0)

#define BINARYITERATOR(X,x,Y,y,exp) do { \
  long _j=0, _k=0, _l=0, _m=0, _n=0; \
  if (x.ndim() == 1) BINARYITERATOR1(X,x,Y,y,exp); \
  else if (x.ndim() == 2) BINARYITERATOR2(X,x,Y,y,exp); \
  else if (x.ndim() == 3) BINARYITERATOR3(X,x,Y,y,exp); \
  else if (x.ndim() == 4) BINARYITERATOR4(X,x,Y,y,exp); \
  else if (x.ndim() == 5) BINARYITERATOR5(X,x,Y,y,exp); \
  else if (x.ndim() == 6) BINARYITERATOR6(X,x,Y,y,exp); \
  else {TENSOR_ASSERT(x.ndim() <= 6,"ndim confused?",x.ndim(),&x);} \
} while(0)

// Inside iterator access pointers to current elements as _p0, _p1, _p2
// _i, _j, _k, ... also defined
#define TERNARYITERATOR1(X,x,Y,y,Z,z,exp) do { \
TENSOR_ASSERT(x.conforms(y),"first and second tensors do not conform",0,&x); \
TENSOR_ASSERT(x.conforms(z),"first and third tensors do not conform",0,&x); \
 long __xd0=x.dim(0);                                                   \
 long __xs0=x.stride(0);                                                \
 long __ys0=y.stride(0);                                                \
 long __zs0=z.stride(0);                                                \
X* restrict _p0=x.ptr(); \
Y* restrict _p1=y.ptr(); \
Z* restrict _p2=z.ptr(); \
for (long _i=0; _i<__xd0; ++_i, _p0+=__xs0, _p1+=__ys0, _p2+=__zs0) { \
  exp; \
} } while(0)

#define TERNARYITERATOR2(X,x,Y,y,Z,z,exp) do { \
TENSOR_ASSERT(x.conforms(y),"first and second tensors do not conform",0,&x); \
TENSOR_ASSERT(x.conforms(z),"first and third tensors do not conform",0,&x); \
 long __xd0=x.dim(0), __xd1=x.dim(1);                                   \
 long __xs0=x.stride(0), __xs1=x.stride(1);                             \
 long __ys0=y.stride(0), __ys1=y.stride(1);                             \
 long __zs0=z.stride(0), __zs1=z.stride(1);                             \
X* restrict __xp0=x.ptr(); \
Y* restrict __yp0=y.ptr(); \
Z* restrict __zp0=z.ptr(); \
for (long _i=0; _i<__xd0; ++_i, __xp0+=__xs0, __yp0+=__ys0, __zp0+=__zs0) { \
  X* restrict _p0=__xp0; \
  Y* restrict _p1=__yp0; \
  Z* restrict _p2=__zp0; \
  for (long _j=0; _j<__xd1; ++_j, _p0+=__xs1, _p1+=__ys1, _p2+=__zs1) { \
    exp; \
  } } } while(0)

#define TERNARYITERATOR3(X,x,Y,y,Z,z,exp) do { \
TENSOR_ASSERT(x.conforms(y),"first and second tensors do not conform",0,&x); \
TENSOR_ASSERT(x.conforms(z),"first and third tensors do not conform",0,&x); \
 long __xd0=x.dim(0), __xd1=x.dim(1), __xd2=x.dim(2);                   \
 long __xs0=x.stride(0), __xs1=x.stride(1), __xs2=x.stride(2);          \
 long __ys0=y.stride(0), __ys1=y.stride(1), __ys2=y.stride(2);          \
 long __zs0=z.stride(0), __zs1=z.stride(1), __zs2=z.stride(2);          \
X* restrict __xp0=x.ptr(); \
Y* restrict __yp0=y.ptr(); \
Z* restrict __zp0=z.ptr(); \
for (long _i=0; _i<__xd0; ++_i, __xp0+=__xs0, __yp0+=__ys0, __zp0+=__zs0) { \
  X* restrict __xp1=__xp0; \
  Y* restrict __yp1=__yp0; \
  Z* restrict __zp1=__zp0; \
  for (long _j=0; _j<__xd1; ++_j, __xp1+=__xs1, __yp1+=__ys1, __zp1+=__zs1) { \
    X* restrict _p0=__xp1; \
    Y* restrict _p1=__yp1; \
    Z* restrict _p2=__zp1; \
    for (long _k=0; _k<__xd2; ++_k, _p0+=__xs2, _p1+=__ys2, _p2+=__zs2) { \
      exp; \
    } } } } while(0)

#define TERNARYITERATOR4(X,x,Y,y,Z,z,exp) do { \
TENSOR_ASSERT(x.conforms(y),"first and second tensors do not conform",0,&x); \
TENSOR_ASSERT(x.conforms(z),"first and third tensors do not conform",0,&x); \
 long __xd0=x.dim(0), __xd1=x.dim(1), __xd2=x.dim(2),                   \
     __xd3=x.dim(3);                                                    \
 long __xs0=x.stride(0), __xs1=x.stride(1), __xs2=x.stride(2),          \
     __xs3=x.stride(3);                                                 \
 long __ys0=y.stride(0), __ys1=y.stride(1), __ys2=y.stride(2),          \
     __ys3=y.stride(3);                                                 \
 long __zs0=z.stride(0), __zs1=z.stride(1), __zs2=z.stride(2),          \
     __zs3=z.stride(3);                                                 \
X* restrict __xp0=x.ptr(); \
Y* restrict __yp0=y.ptr(); \
Z* restrict __zp0=z.ptr(); \
for (long _i=0; _i<__xd0; ++_i, __xp0+=__xs0, __yp0+=__ys0, __zp0+=__zs0) { \
  X* restrict __xp1=__xp0; \
  Y* restrict __yp1=__yp0; \
  Z* restrict __zp1=__zp0; \
  for (long _j=0; _j<__xd1; ++_j, __xp1+=__xs1, __yp1+=__ys1, __zp1+=__zs1) { \
    X* restrict __xp2=__xp1; \
    Y* restrict __yp2=__yp1; \
    Z* restrict __zp2=__zp1; \
    for (long _k=0; _k<__xd2; ++_k, __xp2+=__xs2, __yp2+=__ys2, __zp2+=__zs2) { \
      X* restrict _p0=__xp2; \
      Y* restrict _p1=__yp2; \
      Z* restrict _p2=__zp2; \
      for (long _l=0; _l<__xd3; ++_l, _p0+=__xs3, _p1+=__ys3, _p2+=__zs3) { \
        exp; \
      } } } } } while(0)

#define TERNARYITERATOR5(X,x,Y,y,Z,z,exp) do { \
TENSOR_ASSERT(x.conforms(y),"first and second tensors do not conform",0,&x); \
TENSOR_ASSERT(x.conforms(z),"first and third tensors do not conform",0,&x); \
 long __xd0=x.dim(0), __xd1=x.dim(1), __xd2=x.dim(2),                   \
     __xd3=x.dim(3), __xd4=x.dim(4);                                    \
 long __xs0=x.stride(0), __xs1=x.stride(1), __xs2=x.stride(2),          \
     __xs3=x.stride(3), __xs4=x.stride(4);                              \
 long __ys0=y.stride(0), __ys1=y.stride(1), __ys2=y.stride(2),          \
     __ys3=y.stride(3), __ys4=y.stride(4);                              \
 long __zs0=z.stride(0), __zs1=z.stride(1), __zs2=z.stride(2),          \
     __zs3=z.stride(3), __zs4=z.stride(4);                              \
X* restrict __xp0=x.ptr(); \
Y* restrict __yp0=y.ptr(); \
Z* restrict __zp0=z.ptr(); \
for (long _i=0; _i<__xd0; ++_i, __xp0+=__xs0, __yp0+=__ys0, __zp0+=__zs0) { \
  X* restrict __xp1=__xp0; \
  Y* restrict __yp1=__yp0; \
  Z* restrict __zp1=__zp0; \
  for (long _j=0; _j<__xd1; ++_j, __xp1+=__xs1, __yp1+=__ys1, __zp1+=__zs1) { \
    X* restrict __xp2=__xp1; \
    Y* restrict __yp2=__yp1; \
    Z* restrict __zp2=__zp1; \
    for (long _k=0; _k<__xd2; ++_k, __xp2+=__xs2, __yp2+=__ys2, __zp2+=__zs2) { \
      X* restrict __xp3=__xp2; \
      Y* restrict __yp3=__yp2; \
      Z* restrict __zp3=__zp2; \
      for (long _l=0; _l<__xd3; ++_l, __xp3+=__xs3, __yp3+=__ys3, __zp3+=__zs3) { \
        X* restrict _p0=__xp3; \
        Y* restrict _p1=__yp3; \
        Z* restrict _p2=__zp3; \
        for (long _m=0; _m<__xd4; ++_m, _p0+=__xs4, _p1+=__ys4, _p2+=__zs4) { \
          exp; \
        } } } } } } while(0)

#define TERNARYITERATOR6(X,x,Y,y,Z,z,exp) do { \
TENSOR_ASSERT(x.conforms(y),"first and second tensors do not conform",0,&x); \
TENSOR_ASSERT(x.conforms(z),"first and third tensors do not conform",0,&x); \
 long __xd0=x.dim(0), __xd1=x.dim(1), __xd2=x.dim(2),                   \
     __xd3=x.dim(3), __xd4=x.dim(4), __xd5=x.dim(5);                    \
 long __xs0=x.stride(0), __xs1=x.stride(1), __xs2=x.stride(2),          \
     __xs3=x.stride(3), __xs4=x.stride(4), __xs5=x.stride(5);           \
 long __ys0=y.stride(0), __ys1=y.stride(1), __ys2=y.stride(2),          \
     __ys3=y.stride(3), __ys4=y.stride(4), __ys5=y.stride(5);           \
 long __zs0=z.stride(0), __zs1=z.stride(1), __zs2=z.stride(2),          \
     __zs3=z.stride(3), __zs4=z.stride(4), __zs5=z.stride(5);           \
X* restrict __xp0=x.ptr(); \
Y* restrict __yp0=y.ptr(); \
Z* restrict __zp0=z.ptr(); \
for (long _i=0; _i<__xd0; ++_i, __xp0+=__xs0, __yp0+=__ys0, __zp0+=__zs0) { \
  X* restrict __xp1=__xp0; \
  Y* restrict __yp1=__yp0; \
  Z* restrict __zp1=__zp0; \
  for (long _j=0; _j<__xd1; ++_j, __xp1+=__xs1, __yp1+=__ys1, __zp1+=__zs1) { \
    X* restrict __xp2=__xp1; \
    Y* restrict __yp2=__yp1; \
    Z* restrict __zp2=__zp1; \
    for (long _k=0; _k<__xd2; ++_k, __xp2+=__xs2, __yp2+=__ys2, __zp2+=__zs2) { \
      X* restrict __xp3=__xp2; \
      Y* restrict __yp3=__yp2; \
      Z* restrict __zp3=__zp2; \
      for (long _l=0; _l<__xd3; ++_l, __xp3+=__xs3, __yp3+=__ys3, __zp3+=__zs3) { \
        X* restrict __xp4=__xp3; \
        Y* restrict __yp4=__yp3; \
        Z* restrict __zp4=__zp3; \
        for (long _m=0; _m<__xd4; ++_m, __xp4+=__xs4, __yp4+=__ys4, __zp4+=__zs4) { \
          X* restrict _p0=__xp4; \
          Y* restrict _p1=__yp4; \
          Z* restrict _p2=__zp4; \
          for (long _n=0; _n<__xd5; ++_n, _p0+=__xs5, _p1+=__ys5, _p2+=__zs5) { \
            exp; \
          } } } } } } } while(0)

#define TERNARYITERATOR(X,x,Y,y,Z,z,exp) do { \
  long _j=0, _k=0, _l=0, _m=0, _n=0; \
  if (x.ndim() == 1) TERNARYITERATOR1(X,x,Y,y,Z,z,exp); \
  else if (x.ndim() == 2) TERNARYITERATOR2(X,x,Y,y,Z,z,exp); \
  else if (x.ndim() == 3) TERNARYITERATOR3(X,x,Y,y,Z,z,exp); \
  else if (x.ndim() == 4) TERNARYITERATOR4(X,x,Y,y,Z,z,exp); \
  else if (x.ndim() == 5) TERNARYITERATOR5(X,x,Y,y,Z,z,exp); \
  else if (x.ndim() == 6) TERNARYITERATOR6(X,x,Y,y,Z,z,exp); \
  else {TENSOR_ASSERT(x.ndim() <= 6,"ndim confused?",x.ndim(),&x);} \
} while(0)

#define UNARY_OPTIMIZED_ITERATOR(X,x,exp) do { \
  if (x.iscontiguous()) { \
    X* restrict _p0 = x.ptr(); \
    for (long _j=0; _j<x.size(); ++_j,++_p0) {exp;} \
  } \
  else { \
    for (TensorIterator<REMCONST(X)> iter=x.unary_iterator(1); iter._p0; ++iter) { \
      long _dimj = iter.dimj; \
      X* restrict _p0 = iter._p0; \
      long _s0 = iter._s0; \
      for (long _j=0; _j<_dimj; ++_j, _p0+=_s0) { \
        exp; \
      } \
    } \
  } \
} while(0)


// Can optimize these by moving definition of stride out of loop.

#define UNARY_UNOPTIMIZED_ITERATOR(X,x,exp) do { \
    for (TensorIterator<REMCONST(X)> iter=x.unary_iterator(1,false,false); iter._p0; ++iter) { \
    long _dimj = iter.dimj; \
    X* restrict _p0 = iter._p0; \
    long _s0 = iter._s0; \
    for (long _j=0; _j<_dimj; ++_j, _p0+=_s0) { \
      exp; \
    } \
  } } while(0)

// Use this inside another iterator macro ... pointers are _q* instead of _p*
// Iterator is iter2.

#define UNARY_UNOPTIMIZED_ITERATOR_NESTED(X,x,exp) do { \
    for (TensorIterator<REMCONST(X)> iter2=x.unary_iterator(1,false,false); iter2._p0; ++iter2) { \
    long _dimj2 = iter2.dimj; \
    X* restrict _q0 = iter2._p0; \
    long _s20 = iter2._s0; \
    for (long _j2=0; _j2<_dimj2; ++_j2, _q0+=_s20) { \
      exp; \
    } \
  } } while(0)

#define BINARY_OPTIMIZED_ITERATOR(X,x,Y,y,exp) do { \
  if (x.iscontiguous() && y.iscontiguous() && x.size()==y.size()) { \
    X* restrict _p0 = x.ptr(); \
    Y* restrict _p1 = y.ptr(); \
    for (long _j=0; _j<x.size(); ++_j,++_p0,++_p1) {exp;} \
  } \
  else { \
    for (TensorIterator<REMCONST(X),REMCONST(Y)> iter=x.binary_iterator(y,1); iter._p0; ++iter) { \
        long _dimj = iter.dimj; \
        X* restrict _p0 = iter._p0; \
        Y* restrict _p1 = iter._p1; \
        long _s0 = iter._s0; \
        long _s1 = iter._s1; \
        for (long _j=0; _j<_dimj; ++_j, _p0+=_s0, _p1+=_s1) { \
          exp; \
        } \
  } } } while(0)

#define TERNARY_OPTIMIZED_ITERATOR(X,x,Y,y,Z,z,exp) do { \
  if (x.iscontiguous() && y.iscontiguous() && z.iscontiguous() && x.size()==y.size() && x.size()==z.size()) { \
    X* restrict _p0 = x.ptr(); \
    Y* restrict _p1 = y.ptr(); \
    Z* restrict _p2 = z.ptr(); \
    for (long _j=0; _j<x.size(); ++_j,++_p0,++_p1,++_p2) {exp;} \
  } \
  else { \
    for (TensorIterator<REMCONST(X),REMCONST(Y),REMCONST(Z)> iter=x.ternary_iterator(y,z,1); iter._p0; ++iter) { \
        long _dimj = iter.dimj; \
        X* restrict _p0 = iter._p0; \
        Y* restrict _p1 = iter._p1; \
        Z* restrict _p2 = iter._p2; \
        long _s0 = iter._s0; \
        long _s1 = iter._s1; \
        long _s2 = iter._s2; \
        for (long _j=0; _j<_dimj; ++_j, _p0+=_s0, _p1+=_s1, _p2+=_s2) { \
          exp; \
        } \
  } } } while(0)

#endif // MADNESS_TENSOR_TENSOR_MACROS_H__INCLUDED
