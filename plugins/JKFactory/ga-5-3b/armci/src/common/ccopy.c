/***************************************************************************

                  COPYRIGHT

The following is a notice of limited availability of the code, and disclaimer
which must be included in the prologue of the code and in all source listings
of the code.

Copyright Notice
 + 2009 University of Chicago

Permission is hereby granted to use, reproduce, prepare derivative works, and
to redistribute to others.  This software was authored by:

Jeff R. Hammond
Leadership Computing Facility
Argonne National Laboratory
Argonne IL 60439 USA
phone: (630) 252-5381
e-mail: jhammond@anl.gov

                  GOVERNMENT LICENSE

Portions of this material resulted from work developed under a U.S.
Government Contract and are subject to the following license: the Government
is granted for itself and others acting on its behalf a paid-up, nonexclusive,
irrevocable worldwide license in this computer software to reproduce, prepare
derivative works, and perform publicly and display publicly.

                  DISCLAIMER

This computer code material was prepared, in part, as an account of work
sponsored by an agency of the United States Government.  Neither the United
States, nor the University of Chicago, nor any of their employees, makes any
warranty express or implied, or assumes any legal liability or responsibility
for the accuracy, completeness, or usefulness of any information, apparatus,
product, or process disclosed, or represents that its use would not infringe
privately owned rights.

 ***************************************************************************/
#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include "copy.h"

/* ONE-DIMENSIONAL COPY OPERATIONS */

#if 0
      subroutine dcopy1d_n(A, B, n)
      integer n,i 
      double precision A(n), B(n)
ccdir$ no_cache_alloc a,b
      do i = 1, n 
            B(i) = A(i)
      end do
      end
#endif

void c_dcopy1d_n_(const double* const restrict A,
                        double* const restrict B,
                  const int*    const restrict n)
{
    int i;
    for ( i = 0 ; i < (*n) ; i++ ){
        B[i] = A[i];
    }
    return;
}

#if 0
      subroutine dcopy1d_u(A, B, n)
      integer n,n1,i
      double precision A(n), B(n)
      double precision d1, d2, d3, d4
      n1 = iand(max0(n,0),3)
      do i = 1, n1
            B(i) = A(i)
      end do
      do i = n1+1, n, 4
         d1 = a(i)
         d2 = a(i+1)
         d3 = a(i+2)
         d4 = a(i+3)
         b(i) = d1
         b(i+1) = d2
         b(i+2) = d3
         b(i+3) = d4
      end do
      end
#endif

void c_dcopy1d_u_(const double* const restrict A,
                        double* const restrict B,
                  const int*    const restrict n)
{
    int i;
    int m = (*n) - ((*n)%4);
    for ( i = 0 ; i < m ; i+=4 ){
        B[i  ] = A[i  ];
        B[i+1] = A[i+1];
        B[i+2] = A[i+2];
        B[i+3] = A[i+3];
    }
    for ( i = m ; i < (*n) ; i++ ){
        B[i] = A[i];
    }
    return;
}

/* TWO-DIMENSIONAL COPY OPERATIONS */

#if 0
      subroutine dcopy21(rows, cols, A, ald, buf, cur)
      integer rows, cols
      integer c, r, ald, cur 
      double precision A(ald,*), buf(ald) 
      cur = 0
      do c = 1, cols
         do r = 1, rows
            cur = cur+1
            buf(cur) = A(r,c)
         end do
      end do
      end
#endif

void c_dcopy21_(const int*    const restrict rows,
                const int*    const restrict cols,
                const double* const restrict A,
                const int*    const restrict ald,
                      double* const restrict buf,
                      int*    const restrict cur)
{
    int r, c, i=0;
    for ( c = 0 ; c < (*cols) ; c++ ){
        for ( r = 0 ; r < (*rows) ; r++ ){
            buf[i++] = A[ c * (*ald) + r ];
        }
    }
    (*cur) = i;
    return;
}

#if 0
      subroutine dcopy12(rows, cols, A, ald, buf, cur)
      integer rows, cols
      integer c, r, ald, cur
      double precision A(ald,*), buf(ald)
      cur = 0
      do c = 1, cols
         do r = 1, rows
            cur = cur+1
            A(r,c) = buf(cur)
         end do
      end do
      end
#endif

void c_dcopy12_(const int*    const restrict rows,
                const int*    const restrict cols,
                      double* const restrict A,
                const int*    const restrict ald,
                const double* const restrict buf,
                      int*    const restrict cur)
{
    int r, c, i=0;
    i = 0;
    for ( c = 0 ; c < (*cols) ; c++ ){
        for ( r = 0 ; r < (*rows) ; r++ ){
            A[ c * (*ald) + r ] = buf[i++];
        }
    }
    (*cur) = i;
    return;
}

#if 0
      subroutine dcopy2d_n(rows, cols, A, ald, B, bld)
      integer rows, cols
      integer c, r, ald, bld
      double precision A(ald,*), B(bld,*)
      do c = 1, cols
         do r = 1, rows
            B(r,c) = A(r,c)
         end do
      end do
      end
#endif

void c_dcopy2d_n_(const int*    const restrict rows,
                  const int*    const restrict cols,
                  const double* const restrict A,
                  const int*    const restrict ald,
                        double* const restrict B,
                  const int*    const restrict bld)
{
    int r, c;
    for ( c = 0 ; c < (*cols) ; c++ ){
        for ( r = 0 ; r < (*rows) ; r++ ){
            B[ c * (*bld) + r ] = A[ c * (*ald) + r ];
        }
    }
    return;
}

#if 0
      subroutine dcopy2d_u(rows, cols, A, ald, B, bld)
      integer rows, cols
      integer c, r, ald, bld
      double precision A(ald,*), B(bld,*)
      integer r1
      double precision d1, d2, d3, d4
      do c = 1, cols
      r1 = iand(max0(rows,0),3)
      do r = 1, r1
c$$$         b(r,c) = a(r,c) + b(r,c) * 0
         b(r,c) = a(r,c)
      end do
      do r = r1 + 1, rows, 4
         d1 = a(r,c)
         d2 = a(r+1,c)
         d3 = a(r+2,c)
         d4 = a(r+3,c)
         b(r,c) = d1
         b(r+1,c) = d2
         b(r+2,c) = d3
         b(r+3,c) = d4
c$$$         b(r,c) = a(r,c) + b(r,c) * 0
c$$$         b(r+1,c) = a(r+1,c) + b(r+1,c) * 0
c$$$         b(r+2,c) = a(r+2,c) + b(r+2,c) * 0
c$$$         b(r+3,c) = a(r+3,c) + b(r+3,c) * 0
      enddo
      enddo
      end
#endif

void c_dcopy2d_u_(const int*    const restrict rows,
                  const int*    const restrict cols,
                  const double* const restrict A,
                  const int*    const restrict ald,
                        double* const restrict B,
                  const int*    const restrict bld)
{
    int r, c;
    for ( c = 0 ; c < (*cols) ; c++ ){
        int m = (*rows) - ((*rows)%4);
        for ( r = 0 ; r < m ; r+=4 ){
            B[ c * (*bld) + r   ] = A[ c * (*ald) + r   ];
            B[ c * (*bld) + r+1 ] = A[ c * (*ald) + r+1 ];
            B[ c * (*bld) + r+2 ] = A[ c * (*ald) + r+2 ];
            B[ c * (*bld) + r+3 ] = A[ c * (*ald) + r+3 ];
        }
        for ( r = m ; r < (*rows) ; r++ ){
            B[ c * (*bld) + r   ] = A[ c * (*ald) + r   ];
        }
    }
    return;
}

/* THREE-DIMENSIONAL COPY OPERATIONS */

#if 0
      subroutine dcopy31(rows, cols, planes, A, aldr, aldc, buf, cur)
      integer rows, cols, planes
      integer c, r, p, aldr, aldc, cur
      double precision A(aldr, aldc, *), buf(aldr)
      cur = 0
      do p = 1, planes 
         do c = 1, cols
            do r = 1, rows
               cur = cur+1
               buf(cur) = A(r,c,p)
            end do
         end do
      end do
      end
#endif

void c_dcopy31_(const int*    const restrict rows,
                const int*    const restrict cols,
                const int*    const restrict plns,
                const double* const restrict A,
                const int*    const restrict aldr,
                const int*    const restrict aldc,
                      double* const restrict buf,
                      int*    const restrict cur)
{
    int r, c, p, i=0;
    for ( p = 0 ; p < (*plns) ; p++ ){
        for ( c = 0 ; c < (*cols) ; c++ ){
            for ( r = 0 ; r < (*rows) ; r++ ){
                buf[i++] = A[ p * (*aldc) * (*aldr) + c * (*aldr) + r ];
            }
        }
    }
    (*cur) = i;
    return;
}

#if 0
      subroutine dcopy13(rows, cols, planes, A, aldr, aldc, buf, cur)
      integer rows, cols, planes
      integer c, r, p, aldr, aldc, cur
      double precision A(aldr, aldc, *), buf(aldr)
      cur = 0
      do p = 1, planes
         do c = 1, cols
            do r = 1, rows
               cur = cur+1
               A(r,c,p) = buf(cur)
            end do
         end do
      end do
      end
#endif

void c_dcopy13_(const int*    const restrict rows,
                const int*    const restrict cols,
                const int*    const restrict plns,
                      double* const restrict A,
                const int*    const restrict aldr,
                const int*    const restrict aldc,
                const double* const restrict buf,
                      int*    const restrict cur)
{
    int r, c, p, i=0;
    for ( p = 0 ; p < (*plns) ; p++ ){
        for ( c = 0 ; c < (*cols) ; c++ ){
            for ( r = 0 ; r < (*rows) ; r++ ){
                A[ p * (*aldc) * (*aldr) + c * (*aldr) + r ] = buf[i++];
            }
        }
    }
    (*cur) = i;
    return;
}
