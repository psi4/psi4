/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
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
 * @END LICENSE
 */

/*!
** \file
** \brief Wrapper for the mmult function
** \ingroup CIOMR
*/
 
#include <cstdio>
#include <cstdlib>
#include "libciomr.h"

namespace psi {

  static void mxmbol(double **, int, int, double **, int, int, 
    double **, int, int, int, int, int) 
  {
  abort();
  }

  
/*!
** mxmb: multiplies two rectangular matrices together (wrapper for mmult).  
** Deprecated; use C_DGEMM instead.
**
** \param a    = first matrix to multiply
** \param ia   = if 1, normal multiplication of a
** \param ja   = if 1, transpose a before multiplication
** \param b    = second matrix to multiply
** \param ib   = if 1, normal multiplication of b
** \param jb   = if 1, transpose b before multiplication
** \param c    = matrix to store the result
** \param ic   = if 1, normal multiplication into c
** \param jc   = if 1, transpose c after multiplication
** \param nrow = number of rows of a
** \param nlnk = number of columns of a and rows of b
** \param ncol = number of columns of b
** 
** Returns: none
**
** \ingroup CIOMR
*/
void mxmb(double **a, int ia, int ja, double **b, int ib, int jb, 
          double **c, int ic, int jc, int nrow, int nlnk, int ncol)
   {
      if (ic == 1) {
         if (ia == 1) {
            if (ib == 1) {
               mmult(a,0,b,0,c,0,nrow,nlnk,ncol,0);
               }
            else {
               if (jb == 1) {
                  mmult(a,0,b,1,c,0,nrow,nlnk,ncol,0);
                  }
               else {
                  mxmbol(a,ia,ja,b,ib,jb,c,ic,jc,nrow,nlnk,ncol);
                  }
               }
            }
         else {
            if (ja == 1) {
               if (ib == 1) {
                  mmult(a,1,b,0,c,0,nrow,nlnk,ncol,0);
                  }
               else {
                  if (jb == 1) {
                     mmult(a,1,b,1,c,0,nrow,nlnk,ncol,0);
                     }
                  else {
                     mxmbol(a,ia,ja,b,ib,jb,c,ic,jc,nrow,nlnk,ncol);
                     }
                  }
               }
            else {
               mxmbol(a,ia,ja,b,ib,jb,c,ic,jc,nrow,nlnk,ncol);
               }
            }
         }
      else {
         if (jc == 1) {
            if (ia == 1) {
               if (ib == 1) {
                  mmult(a,0,b,0,c,1,nrow,nlnk,ncol,0);
                  }
               else {
                  if (jb == 1) {
                     mmult(a,0,b,1,c,1,nrow,nlnk,ncol,0);
                     }
                  else {
                     mxmbol(a,ia,ja,b,ib,jb,c,ic,jc,nrow,nlnk,ncol);
                     }
                  }
               }
            else {
               if (ja == 1) {
                  if (ib == 1) {
                     mmult(a,1,b,0,c,1,nrow,nlnk,ncol,0);
                     }
                  else {
                     if (jb == 1) {
                        mmult(a,1,b,1,c,1,nrow,nlnk,ncol,0);
                        }
                     else {
                        mxmbol(a,ia,ja,b,ib,jb,c,ic,jc,nrow,nlnk,ncol);
                        }
                     }
                  }
               else {
                  mxmbol(a,ia,ja,b,ib,jb,c,ic,jc,nrow,nlnk,ncol);
                  }
               }
            }
         else {
            mxmbol(a,ia,ja,b,ib,jb,c,ic,jc,nrow,nlnk,ncol);
            }
         }
      }

}
