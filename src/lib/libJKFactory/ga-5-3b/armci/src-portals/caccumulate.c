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

/***********************************************************************
 *     accumulate operation for the following datatypes:
 *            real, double precision, complex, double complex, integer
 *
 *     WARNING: This file must be compiled WITH optimization under AIX.
 *              IBM fortran compilers generate bad code with -g option. 
 *
 *     Two versions of each routine are provided: 
 *         original and unrolled loops.
 *
 ***********************************************************************/
#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include "acc.h"

#if 0
      subroutine d_accumulate_1d(alpha,  A,  B, rows)
      integer rows, r
      double precision A(*), B(*), alpha
ccdir$ no_cache_alloc a,b
         do r = 1, rows
            A(r) = A(r)+ alpha*B(r)
         enddo
      end
#endif

void c_d_accumulate_1d_(const double* const restrict alpha,
                              double* restrict A,
                        const double* const restrict B,
                        const int*    const restrict rows)
{
    int i;
    for ( i = 0 ; i < (*rows) ; i++ ){
        A[i] += (*alpha) * B[i];
    }
    return;
}


#if 0
      subroutine d_accumulate_2d(alpha, rows, cols, A, ald, B, bld)
      integer rows, cols
      integer c, r, ald, bld
      double precision A(ald,*), B(bld,*), alpha
ccdir$ no_cache_alloc a,b
      do c = 1, cols
         do r = 1, rows
            A(r,c) = A(r,c)+ alpha*B(r,c)
         enddo
      enddo
      end
#endif

void c_d_accumulate_2d_(const double* const alpha,
                        const int* const rows,
                        const int* const cols,
                        double* restrict A,
                        const int* const ald,
                        const double* const B,
                        const int* const bld)
{
    int r, c;
    for ( c = 0 ; c < (*cols) ; c++ ){
        for ( r = 0 ; r < (*rows) ; r++ ){
            A[ c * (*ald) + r ] += (*alpha) * B[ c * (*bld) + r ];
        }
    }
    return;
}

#if 0
      subroutine f_accumulate_1d(alpha,  A,  B, rows)
      integer rows, r
      real A(*), B(*), alpha
         do r = 1, rows
            A(r) = A(r)+ alpha*B(r)
         enddo
      end
#endif

void c_f_accumulate_1d_(const float* const restrict alpha,
                              float* const restrict A,
                        const float* const restrict B,
                        const int*   const restrict rows)
{
    int i;
    for ( i = 0 ; i < (*rows) ; i++ ){
        A[i] += (*alpha) * B[i];
    }
    return;
}

#if 0
      subroutine f_accumulate_2d(alpha, rows, cols, A, ald, B, bld)
      integer rows, cols
      integer c, r, ald, bld
      real A(ald,*), B(bld,*), alpha
      do c = 1, cols
         do r = 1, rows
            A(r,c) = A(r,c)+ alpha*B(r,c)
         enddo
      enddo
      end
#endif

void c_f_accumulate_2d_(const float* const alpha,
                        const int* const rows,
                        const int* const cols,
                        float* restrict A,
                        const int* const ald,
                        const float* const B,
                        const int* const bld)
{
    int r, c;
    for ( c = 0 ; c < (*cols) ; c++ ){
        for ( r = 0 ; r < (*rows) ; r++ ){
            A[ c * (*ald) + r ] += (*alpha) * B[ c * (*bld) + r ];
        }
    }
    return;
}

#if 0
      subroutine z_accumulate_1d(alpha,  A,  B, rows)
      integer rows, r
      double complex  A(*), B(*), alpha
         do r = 1, rows
            A(r) = A(r)+ alpha*B(r)
         enddo
      end
#endif

void c_c_accumulate_1d_(const complex_t* const restrict alpha,
                              complex_t* const restrict A,
                        const complex_t* const restrict B,
                        const int*       const restrict rows)
{
    int i;
    for ( i = 0 ; i < (*rows) ; i++ ){
        A[i].real += (*alpha).real * B[i].real - (*alpha).imag * B[i].imag;
        A[i].imag += (*alpha).imag * B[i].real + (*alpha).real * B[i].imag;
    }
    return;
}

#if 0
      subroutine z_accumulate_2d(alpha, rows, cols, A, ald, B, bld)
      integer rows, cols
      integer c, r, ald, bld
      double complex A(ald,*), B(bld,*), alpha
      do c = 1, cols
         do r = 1, rows
            A(r,c) = A(r,c)+ alpha*B(r,c)
         enddo
      enddo
      end
#endif

void c_c_accumulate_2d_(const complex_t* const alpha,
                        const int* const rows,
                        const int* const cols,
                        complex_t* restrict A,
                        const int* const ald,
                        const complex_t* const B,
                        const int* const bld)
{
    int r, c;
    for ( c = 0 ; c < (*cols) ; c++ ) {
        for ( r = 0 ; r < (*rows) ; r++ ) {
            A[ c * (*ald) + r ].real += (*alpha).real * B[ c * (*bld) + r ].real - (*alpha).imag * B[ c * (*bld) + r ].imag;
            A[ c * (*ald) + r ].imag += (*alpha).imag * B[ c * (*bld) + r ].real + (*alpha).real * B[ c * (*bld) + r ].imag;
        }
    }
    return;
}

#if 0
      subroutine c_accumulate_1d(alpha,  A,  B, rows)
      integer rows, r
      complex  A(*), B(*), alpha
         do r = 1, rows
            A(r) = A(r)+ alpha*B(r)
         enddo
      end
#endif

void c_z_accumulate_1d_(const dcomplex_t* const restrict alpha,
                              dcomplex_t* const restrict A,
                        const dcomplex_t* const restrict B,
                        const int*        const restrict rows)
{
    int i;
    for ( i = 0 ; i < (*rows) ; i++ ){
        A[i].real += (*alpha).real * B[i].real - (*alpha).imag * B[i].imag;
        A[i].imag += (*alpha).imag * B[i].real + (*alpha).real * B[i].imag;
    }
    return;
}

#if 0
      subroutine c_accumulate_2d(alpha, rows, cols, A, ald, B, bld)
      integer rows, cols
      integer c, r, ald, bld
      complex A(ald,*), B(bld,*), alpha
      do c = 1, cols
         do r = 1, rows
            A(r,c) = A(r,c)+ alpha*B(r,c)
         enddo
      enddo
      end
#endif

void c_z_accumulate_2d_(const dcomplex_t* const alpha,
                        const int* const rows,
                        const int* const cols,
                        dcomplex_t* restrict A,
                        const int* const ald,
                        const dcomplex_t* const B,
                        const int* const bld)
{
    int r, c;
    for ( c = 0 ; c < (*cols) ; c++ ) {
        for ( r = 0 ; r < (*rows) ; r++ ) {
            A[ c * (*ald) + r ].real += (*alpha).real * B[ c * (*bld) + r ].real - (*alpha).imag * B[ c * (*bld) + r ].imag;
            A[ c * (*ald) + r ].imag += (*alpha).imag * B[ c * (*bld) + r ].real + (*alpha).real * B[ c * (*bld) + r ].imag;
        }
    }
    return;
}

#if 0
      subroutine i_accumulate_2d(alpha, rows, cols, A, ald, B, bld)
      integer rows, cols
      integer c, r, ald, bld
      integer A(ald,*), B(bld,*), alpha
      do c = 1, cols
         do r = 1, rows
            A(r,c) = A(r,c)+ alpha*B(r,c)
         enddo
      enddo
      end
#endif

void c_i_accumulate_1d_(const int* const restrict alpha,
                              int* const restrict A,
                        const int* const restrict B,
                        const int* const restrict rows)
{
    int i;
    for ( i = 0 ; i < (*rows) ; i++ ){
        A[i] += (*alpha) * B[i];
    }
    return;
}

void c_l_accumulate_1d_(const long* const restrict alpha,
                              long* const restrict A,
                        const long* const restrict B,
                        const int*  const restrict rows)
{
    int i;
    for ( i = 0 ; i < (*rows) ; i++ ){
        A[i] += (*alpha) * B[i];
    }
    return;
}

void c_ll_accumulate_1d_(const long long* const restrict alpha,
                               long long* const restrict A,
                         const long long* const restrict B,
                         const int*       const restrict rows)
{
    int i;
    for ( i = 0 ; i < (*rows) ; i++ ){
        A[i] += (*alpha) * B[i];
    }
    return;
}

#if 0
      subroutine i_accumulate_1d(alpha,  A,  B, rows)
      integer rows, r
      integer A(*), B(*), alpha
         do r = 1, rows
            A(r) = A(r)+ alpha*B(r)
         enddo
      end
#endif

void c_i_accumulate_2d_(const int* const alpha,
                        const int* const rows,
                        const int* const cols,
                        int* restrict A,
                        const int* const ald,
                        const int* const B,
                        const int* const bld)
{
    int r, c;
    for ( c = 0 ; c < (*cols) ; c++ ){
        for ( r = 0 ; r < (*rows) ; r++ ){
            A[ c * (*ald) + r ] += (*alpha) * B[ c * (*bld) + r ];
        }
    }
    return;
}

void c_l_accumulate_2d_(const long* const alpha,
                        const int* const rows,
                        const int* const cols,
                        long* restrict A,
                        const int* const ald,
                        const long* const B,
                        const int* const bld)
{
    int r, c;
    for ( c = 0 ; c < (*cols) ; c++ ){
        for ( r = 0 ; r < (*rows) ; r++ ){
            A[ c * (*ald) + r ] += (*alpha) * B[ c * (*bld) + r ];
        }
    }
    return;
}

void c_ll_accumulate_2d_(const long long* const alpha,
                        const int* const rows,
                        const int* const cols,
                        long long* restrict A,
                        const int* const ald,
                        const long long* const B,
                        const int* const bld)
{
    int r, c;
    for ( c = 0 ; c < (*cols) ; c++ ){
        for ( r = 0 ; r < (*rows) ; r++ ){
            A[ c * (*ald) + r ] += (*alpha) * B[ c * (*bld) + r ];
        }
    }
    return;
}

#if 0
      subroutine d_accumulate_2d_u(alpha, rows, cols, A, ald, B, bld)
      integer rows, cols
      integer c, r, ald, bld
      double precision A(ald,*), B(bld,*), alpha
      integer r1
      doubleprecision d1, d2, d3, d4
      do c = 1, cols
        r1 = iand(max0(rows,0),3)
        do r = 1, r1
            a(r,c) = a(r,c) + alpha*b(r,c)
        end do
        do r = r1 + 1, rows, 4
            d1 = a(r,c) + alpha*b(r,c)
            d2 = a(r+1,c) + alpha*b(r+1,c)
            d3 = a(r+2,c) + alpha*b(r+2,c)
            d4 = a(r+3,c) + alpha*b(r+3,c)
            a(r,c) = d1
            a(r+1,c) = d2
            a(r+2,c) = d3
            a(r+3,c) = d4
        enddo
      enddo
      end
#endif

void c_d_accumulate_2d_u_(const double* const alpha,
                          const int* const rows,
                          const int* const cols,
                          double* restrict A,
                          const int* const ald,
                          const double* const B,
                          const int* const bld)
{
    int r, c;
    int m = (*rows) - ((*rows)%4);
    for ( c = 0 ; c < (*cols) ; c++ ){
        for ( r = 0 ; r < m ; r+=4 ){
            A[ c * (*ald) + r   ] += (*alpha) * B[ c * (*bld) + r   ];
            A[ c * (*ald) + r+1 ] += (*alpha) * B[ c * (*bld) + r+1 ];
            A[ c * (*ald) + r+2 ] += (*alpha) * B[ c * (*bld) + r+2 ];
            A[ c * (*ald) + r+3 ] += (*alpha) * B[ c * (*bld) + r+3 ];
        }
        for ( r = m ; r < (*rows) ; r++ ){
            A[ c * (*ald) + r ] += (*alpha) * B[ c * (*bld) + r ];
        }
    }
    return;
}

#if 0
      subroutine f_accumulate_2d_u(alpha, rows, cols, A, ald, B, bld)
      integer rows, cols
      integer c, r, ald, bld
      real A(ald,*), B(bld,*), alpha
      integer r1
      real d1, d2, d3, d4
      do c = 1, cols
      r1 = iand(max0(rows,0),3)
      do r = 1, r1
         a(r,c) = a(r,c) + alpha*b(r,c)
      end do
      do r = r1 + 1, rows, 4
         d1 = a(r,c) + alpha*b(r,c)
         d2 = a(r+1,c) + alpha*b(r+1,c)
         d3 = a(r+2,c) + alpha*b(r+2,c)
         d4 = a(r+3,c) + alpha*b(r+3,c)
         a(r,c) = d1
         a(r+1,c) = d2
         a(r+2,c) = d3
         a(r+3,c) = d4
      enddo
      enddo
      end
#endif

void c_f_accumulate_2d_u_(const float* const alpha,
                          const int* const rows,
                          const int* const cols,
                          float* restrict A,
                          const int* const ald,
                          const float* const B,
                          const int* const bld)
{
    int r, c;
    int m = (*rows) - ((*rows)%4);
    for ( c = 0 ; c < (*cols) ; c++ ){
        for ( r = 0 ; r < m ; r+=4 ){
            A[ c * (*ald) + r   ] += (*alpha) * B[ c * (*bld) + r   ];
            A[ c * (*ald) + r+1 ] += (*alpha) * B[ c * (*bld) + r+1 ];
            A[ c * (*ald) + r+2 ] += (*alpha) * B[ c * (*bld) + r+2 ];
            A[ c * (*ald) + r+3 ] += (*alpha) * B[ c * (*bld) + r+3 ];
        }
        for ( r = m ; r < (*rows) ; r++ ){
            A[ c * (*ald) + r ] += (*alpha) * B[ c * (*bld) + r ];
        }
    }
    return;
}

#if 0
      subroutine z_accumulate_2d_u(alpha, rows, cols, A, ald, B, bld)
      integer rows, cols
      integer c, r, ald, bld
      double complex A(ald,*), B(bld,*), alpha
      integer r1
      double complex x1, x2, x3, x4
      do c = 1, cols
      r1 = iand(max0(rows,0),3)
      do r = 1, r1
         a(r,c) = a(r,c) + alpha*b(r,c)
      end do
      do r = r1 + 1, rows, 4
         x1 = a(r,c) + alpha*b(r,c)
         x2 = a(r+1,c) + alpha*b(r+1,c)
         x3 = a(r+2,c) + alpha*b(r+2,c)
         x4 = a(r+3,c) + alpha*b(r+3,c)
         a(r,c) = x1
         a(r+1,c) = x2
         a(r+2,c) = x3
         a(r+3,c) = x4
      enddo
      enddo
      end
#endif

void c_c_accumulate_2d_u_(const complex_t* const alpha,
                          const int* const rows,
                          const int* const cols,
                          complex_t* restrict A,
                          const int* const ald,
                          const complex_t* const B,
                          const int* const bld)
{
    int r, c;
    int jA, jB;
    int m = (*rows) - ((*rows)%4);
    for ( c = 0 ; c < (*cols) ; c++ ){
        for ( r = 0 ; r < m ; r+=4 ){
            jA = c * (*ald) + r;
            jB = c * (*bld) + r;
            A[ jA   ].real += (*alpha).real * B[ jB   ].real - (*alpha).imag * B[ jB   ].imag;
            A[ jA   ].imag += (*alpha).imag * B[ jB   ].real + (*alpha).real * B[ jB   ].imag;
            A[ jA+1 ].real += (*alpha).real * B[ jB+1 ].real - (*alpha).imag * B[ jB+1 ].imag;
            A[ jA+1 ].imag += (*alpha).imag * B[ jB+1 ].real + (*alpha).real * B[ jB+1 ].imag;
            A[ jA+2 ].real += (*alpha).real * B[ jB+2 ].real - (*alpha).imag * B[ jB+2 ].imag;
            A[ jA+2 ].imag += (*alpha).imag * B[ jB+2 ].real + (*alpha).real * B[ jB+2 ].imag;
            A[ jA+3 ].real += (*alpha).real * B[ jB+3 ].real - (*alpha).imag * B[ jB+3 ].imag;
            A[ jA+3 ].imag += (*alpha).imag * B[ jB+3 ].real + (*alpha).real * B[ jB+3 ].imag;
        }
        for ( r = m ; r < (*rows) ; r++ ){
            A[ c * (*ald) + r ].real += (*alpha).real * B[ c * (*bld) + r ].real - (*alpha).imag * B[ c * (*bld) + r ].imag;
            A[ c * (*ald) + r ].imag += (*alpha).imag * B[ c * (*bld) + r ].real + (*alpha).real * B[ c * (*bld) + r ].imag;
        }
    }
    return;
}

#if 0
      subroutine c_accumulate_2d_u(alpha, rows, cols, A, ald, B, bld)
      integer rows, cols
      integer c, r, ald, bld
      complex A(ald,*), B(bld,*), alpha
      integer r1
      complex x1, x2, x3, x4
      do c = 1, cols
      r1 = iand(max0(rows,0),3)
      do r = 1, r1
         a(r,c) = a(r,c) + alpha*b(r,c)
      end do
      do r = r1 + 1, rows, 4
         x1 = a(r,c) + alpha*b(r,c)
         x2 = a(r+1,c) + alpha*b(r+1,c)
         x3 = a(r+2,c) + alpha*b(r+2,c)
         x4 = a(r+3,c) + alpha*b(r+3,c)
         a(r,c) = x1
         a(r+1,c) = x2
         a(r+2,c) = x3
         a(r+3,c) = x4
      enddo
      enddo
      end
#endif

void c_z_accumulate_2d_u_(const dcomplex_t* const alpha,
                          const int* const rows,
                          const int* const cols,
                          dcomplex_t* restrict A,
                          const int* const ald,
                          const dcomplex_t* const B,
                          const int* const bld)
{
    int r, c;
    int jA, jB;
    int m = (*rows) - ((*rows)%4);
    for ( c = 0 ; c < (*cols) ; c++ ){
        for ( r = 0 ; r < m ; r+=4 ){
            jA = c * (*ald) + r;
            jB = c * (*bld) + r;
            A[ jA   ].real += (*alpha).real * B[ jB   ].real - (*alpha).imag * B[ jB   ].imag;
            A[ jA   ].imag += (*alpha).imag * B[ jB   ].real + (*alpha).real * B[ jB   ].imag;
            A[ jA+1 ].real += (*alpha).real * B[ jB+1 ].real - (*alpha).imag * B[ jB+1 ].imag;
            A[ jA+1 ].imag += (*alpha).imag * B[ jB+1 ].real + (*alpha).real * B[ jB+1 ].imag;
            A[ jA+2 ].real += (*alpha).real * B[ jB+2 ].real - (*alpha).imag * B[ jB+2 ].imag;
            A[ jA+2 ].imag += (*alpha).imag * B[ jB+2 ].real + (*alpha).real * B[ jB+2 ].imag;
            A[ jA+3 ].real += (*alpha).real * B[ jB+3 ].real - (*alpha).imag * B[ jB+3 ].imag;
            A[ jA+3 ].imag += (*alpha).imag * B[ jB+3 ].real + (*alpha).real * B[ jB+3 ].imag;
        }
        for ( r = m ; r < (*rows) ; r++ ){
            A[ c * (*ald) + r ].real += (*alpha).real * B[ c * (*bld) + r ].real - (*alpha).imag * B[ c * (*bld) + r ].imag;
            A[ c * (*ald) + r ].imag += (*alpha).imag * B[ c * (*bld) + r ].real + (*alpha).real * B[ c * (*bld) + r ].imag;
        }
    }
    return;
}

#if 0
      subroutine i_accumulate_2d_u(alpha, rows, cols, A, ald, B, bld)
      integer rows, cols
      integer c, r, ald, bld
      integer A(ald,*), B(bld,*), alpha

      integer r1, j2, j3, j4, j5
      do c = 1, cols
      r1 = iand(max0(rows,0),3)
      do r = 1, r1
         a(r,c) = a(r,c) + alpha*b(r,c)
      end do
      do r = r1 + 1, rows, 4
         j2 = a(r,c) + alpha*b(r,c)
         j3 = a(r+1,c) + alpha*b(r+1,c)
         j4 = a(r+2,c) + alpha*b(r+2,c)
         j5 = a(r+3,c) + alpha*b(r+3,c)
         a(r,c) = j2
         a(r+1,c) = j3
         a(r+2,c) = j4
         a(r+3,c) = j5
      enddo
      enddo
      end
#endif

void c_i_accumulate_2d_u_(const int* const alpha,
                          const int* const rows,
                          const int* const cols,
                          int* restrict A,
                          const int* const ald,
                          const int* const B,
                          const int* const bld)
{
    int r, c;
    int m = (*rows) - ((*rows)%4);
    for ( c = 0 ; c < (*cols) ; c++ ){
        for ( r = 0 ; r < m ; r+=4 ){
            A[ c * (*ald) + r   ] += (*alpha) * B[ c * (*bld) + r   ];
            A[ c * (*ald) + r+1 ] += (*alpha) * B[ c * (*bld) + r+1 ];
            A[ c * (*ald) + r+2 ] += (*alpha) * B[ c * (*bld) + r+2 ];
            A[ c * (*ald) + r+3 ] += (*alpha) * B[ c * (*bld) + r+3 ];
        }
        for ( r = m ; r < (*rows) ; r++ ){
            A[ c * (*ald) + r ] += (*alpha) * B[ c * (*bld) + r ];
        }
    }
    return;
}

void c_l_accumulate_2d_u_(const long* const alpha,
                          const int* const rows,
                          const int* const cols,
                          long* restrict A,
                          const int* const ald,
                          const long* const B,
                          const int* const bld)
{
    int r, c;
    int m = (*rows) - ((*rows)%4);
    for ( c = 0 ; c < (*cols) ; c++ ){
        for ( r = 0 ; r < m ; r+=4 ){
            A[ c * (*ald) + r   ] += (*alpha) * B[ c * (*bld) + r   ];
            A[ c * (*ald) + r+1 ] += (*alpha) * B[ c * (*bld) + r+1 ];
            A[ c * (*ald) + r+2 ] += (*alpha) * B[ c * (*bld) + r+2 ];
            A[ c * (*ald) + r+3 ] += (*alpha) * B[ c * (*bld) + r+3 ];
        }
        for ( r = m ; r < (*rows) ; r++ ){
            A[ c * (*ald) + r ] += (*alpha) * B[ c * (*bld) + r ];
        }
    }
    return;
}

void c_ll_accumulate_2d_u_(const long long* const alpha,
                           const int* const rows,
                           const int* const cols,
                           long long* restrict A,
                           const int* const ald,
                           const long long* const B,
                           const int* const bld)
{
    int r, c;
    int m = (*rows) - ((*rows)%4);
    for ( c = 0 ; c < (*cols) ; c++ ){
        for ( r = 0 ; r < m ; r+=4 ){
            A[ c * (*ald) + r   ] += (*alpha) * B[ c * (*bld) + r   ];
            A[ c * (*ald) + r+1 ] += (*alpha) * B[ c * (*bld) + r+1 ];
            A[ c * (*ald) + r+2 ] += (*alpha) * B[ c * (*bld) + r+2 ];
            A[ c * (*ald) + r+3 ] += (*alpha) * B[ c * (*bld) + r+3 ];
        }
        for ( r = m ; r < (*rows) ; r++ ){
            A[ c * (*ald) + r ] += (*alpha) * B[ c * (*bld) + r ];
        }
    }
    return;
}

#if 0
c---------- operations used in armci gops --------------
c
      subroutine fort_dadd(n, x, work)
      integer n,i
      double precision x(n), work(n)
      do i= 1,n
         x(i) = x(i) + work(i)
      enddo
      end
#endif

void c_dadd_(const int*    const restrict n,
                   double* const restrict x,
             const double* const restrict work)
{
    int i;
    for ( i = 0 ; i < (*n) ; i++ ){
        x[i] += work[i];
    }
    return;
}

#if 0
      subroutine fort_dadd2(n, x, work, work2)
      integer n,i
      double precision x(n), work(n), work2(n)
      do i= 1,n
         x(i) = work(i) + work2(i)
      enddo
      end
#endif

void c_dadd2_(const int*    const restrict n,
                    double* const restrict x,
              const double* const restrict work,
              const double* const restrict work2)
{
    int i;
    for ( i = 0 ; i < (*n) ; i++ ){
        x[i] = work[i] + work2[i];
    }
    return;
}

#if 0
      subroutine fort_dmult(n, x, work)
      integer n,i
      double precision x(n), work(n)
      do i= 1,n
         x(i) = x(i) * work(i)
      enddo
      end
#endif

void c_dmult_(const int*    const restrict n,
                    double* const restrict x,
              const double* const restrict work)
{
    int i;
    for ( i = 0 ; i < (*n) ; i++ ){
        x[i] *= work[i];
    }
    return;
}

#if 0
      subroutine fort_dmult2(n, x, work,work2)
      integer n,i
      double precision x(n), work(n)
      do i= 1,n
         x(i) = work(i)*work2(i)
      enddo
      end
#endif

void c_dmult2_(const int*    const restrict n,
                     double* const restrict x,
               const double* const restrict work,
               const double* const restrict work2)
{
    int i;
    for ( i = 0 ; i < (*n) ; i++ ){
        x[i] = work[i] * work2[i];
    }
    return;
}


// specific to src-portals && to src-gemini
void  RA_ACCUMULATE_2D_(long* alpha, int* rows, int* cols, long* a,
                      int* lda, long* b, int* ldb)
{
int i,j;
   for(j=0;j< *cols; j++){
     long *aa = a + j* *lda;
     long *bb = b + j* *ldb;
     for(i=0;i< *rows; i++)
       aa[i] ^= bb[i];
   }
}
