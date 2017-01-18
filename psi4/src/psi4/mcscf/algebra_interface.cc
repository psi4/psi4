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

#include "psi4/libmoinfo/libmoinfo.h"
#include "algebra_interface.h"

namespace psi{

 namespace mcscf{

/*
** C_DGEMM_12()
**
** This function calculates C(m,n)=alpha*A(k,m)*B(n,k)+ beta*C(m,n)
**
** nra = number of rows in A
** ncb = number of columns in B
** ncc = number of columns in C
*/
void C_DGEMM_12(int m, int n, int k, double alpha,
           double *A, int nra, double *B, int ncb, double beta, double *C,
           int ncc)
{
  //  the only strange thing we need to do is reverse everything
  //  since the stride runs differently in C vs. Fortran

  /* also, do nothing if a dimension is 0 */
  if (m == 0 || n == 0 || k == 0) return;

  F_DGEMM("t","t",&n,&m,&k,&alpha,B,&ncb,A,&nra,&beta,C,&ncc);
}

/*
** C_DGEMM_22()
**
** This function calculates C(m,n)=alpha*A(m,k)*B(n,k)+ beta*C(m,n)
**
** nra = number of columns in A
** ncb = number of columns in B
** ncc = number of columns in C
*/
void C_DGEMM_22(int m, int n, int k, double alpha,
           double *A, int nca, double *B, int ncb, double beta, double *C,
           int ncc)
{
  //  the only strange thing we need to do is reverse everything
  //  since the stride runs differently in C vs. Fortran

  /* also, do nothing if a dimension is 0 */
  if (m == 0 || n == 0 || k == 0) return;

  F_DGEMM("t","n",&n,&m,&k,&alpha,B,&ncb,A,&nca,&beta,C,&ncc);
}


}} /* End Namespaces */
