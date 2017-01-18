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
** \brief Diagonalize a symmetric matrix in packed (lower triangular) form
** \ingroup CIOMR
*/

#include "psi4/psifiles.h"
#include <cstdlib>
#include "libciomr.h"
#include "psi4/psi4-dec.h"
namespace psi {

extern void tred2(int n,double** a,double* d,double* e,int matz);
extern void tqli(int n, double *d, double **z, double *e, int matz,
  double toler);

/* translation into c of a translation into FORTRAN77 of the EISPACK */
/* matrix diagonalization routines */

/*!
** rsp(): diagonalize a symmetric matrix in packed (lower triangular) form
** in 'array'. For square symmetric matrices, see sq_rsp().
**
** \param nm     = rows of matrix
** \param n      = columns of matrix
** \param nv     = number of elements in lower triangle (n*(n+1)/2)
** \param array  = matrix to diagonalize (packed as linear array)
** \param e_vals = array to hold eigenvalues
** \param matz   = 0 (no eigenvectors, eigenvals in ascending order)
**               = 1 (eigenvectors and eigenvalues in ascending order)
**               = 2 (no eigenvectors, eigenvalues in descending order)
**               = 3 (eigenvectors and eigenvalues in descending order)
** \param e_vecs = matrix of eigenvectors (one column for each eigvector)
** \param toler  = tolerance for eigenvalues?  Often 1.0E-14.
**
** Returns: none
**
** \ingroup CIOMR
*/
void rsp(int nm, int n,int nv,double *array, double *e_vals, int matz,
         double ** e_vecs, double toler)
{
      int i, j, ij, ierr;
      int ascend_order;
      double *fv1=NULL;
      double sw;

      /* Modified by Ed - matz can have values 0 through 3 */

      if ((matz > 3) || (matz < 0)) {
        matz = 0;
        ascend_order = 1;
        }
      else
        if (matz < 2)
          ascend_order = 1;	/* Eigenvalues in ascending order */
        else {
          matz -= 2;
          ascend_order = 0;	/* Eigenvalues in descending order */
          }

      fv1 = (double *) init_array(n);
      /*temp = (double **) init_matrix(n,n);*/

      if (n > nm) {
         ierr = 10*n;
         outfile->Printf("n = %d is greater than nm = %d in rsp\n",n,nm);
         exit(PSI_RETURN_FAILURE);
         }

      if (nv < (n*(n+1)/2)) {
         int num = n*(n+1)/2;
         ierr = 20*n;
         outfile->Printf("nv = %d is less than n*(n+1)/2 = %d in rsp\n",nv,num);
         exit(PSI_RETURN_FAILURE);
         }

      for (i=0,ij=0; i < n; i++) {
         for (j=0; j <= i; j++,ij++) {
            e_vecs[i][j] = array[ij];
            e_vecs[j][i] = array[ij];
            }
          }

      tred2(n,e_vecs,e_vals,fv1,matz);

      for (i=0; i < n; i++)
         for (j=0; j < i; j++){
            sw = e_vecs[i][j];
            e_vecs[i][j] = e_vecs[j][i];
            e_vecs[j][i] = sw;
            /*temp[i][j]=e_vecs[j][i];*/
            }

      tqli(n,e_vals,e_vecs,fv1,matz,toler);
      /*tqli(n,e_vals,temp,fv1,matz,toler);*/

      for (i=0; i < n; i++)
         for (j=0; j < i; j++){
            sw = e_vecs[i][j];
            e_vecs[i][j] = e_vecs[j][i];
            e_vecs[j][i] = sw;
            /*e_vecs[i][j]=temp[j][i];*/
            }

      if (ascend_order)
        eigsort(e_vals,e_vecs,n);
      else
        eigsort(e_vals,e_vecs,-n);

      free(fv1);
      /*free_matrix(temp,n);*/
}

}
