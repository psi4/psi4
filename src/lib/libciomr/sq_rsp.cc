/*!
** \file sq_rsp.cc
** \ingroup (CIOMR)
*/

#include <psifiles.h>
#include "includes.h"
#include <libciomr/libciomr.h>

extern "C" {

extern void tred2(int n,double** a,double* d,double* e,int matz);
extern void tqli(int n, double *d, double **z, double *e, int matz, double toler);
  
/* translation into c of a translation into FORTRAN77 of the EISPACK */
/* matrix diagonalization routines */

/*!
** sq_rsp: diagomalize a symmetric square matrix ('array').
**
**   \param nm     = rows of matrix
**   \param n      = columns of matrix
**   \param nv     = number of elements in lower triangle (n*(n+1)/2)
**   \param array  = matrix to diagonalize
**   \param e_vals = array to hold eigenvalues
**   \param matz   = 0 (no eigenvectors, eigenvals in ascending order)
**                 = 1 (eigenvectors and eigenvalues in ascending order)
**                 = 2 (no eigenvectors, eigenvalues in descending order)
**                 = 3 (eigenvectors and eigenvalues in descending order)
**   \param e_vecs = matrix of eigenvectors (one column for each eigvector)
**   \param toler  = tolerance for eigenvalues?  Often 1.0E-14.
**
** \ingroup (CIOMR)
*/
void sq_rsp(int nm, int n, double **array, double *e_vals, int matz,
            double **e_vecs, double toler)
   {
      int i, j, ii, ij, ierr;
      int ascend_order;
      double *fv1, **temp;
      double zero = 0.0;
      double one = 1.0;

/* Modified by Ed - matz can have the values 0 through 3 */
      
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
      temp = (double **) init_matrix(n,n);

      if (n > nm) {
         ierr = 10*n;
         fprintf(stderr,"n = %d is greater than nm = %d in rsp\n",n,nm);
         exit(PSI_RETURN_FAILURE);
         }

      for (i=0; i < n; i++) {
         for (j=0; j < n; j++) {
            e_vecs[i][j] = array[i][j];
            }
          }

      tred2(n,e_vecs,e_vals,fv1,matz);

      for (i=0; i < n; i++)
         for (j=0; j < n; j++)
            temp[i][j]=e_vecs[j][i];

      tqli(n,e_vals,temp,fv1,matz,toler);

      for (i=0; i < n; i++)
         for (j=0; j < n; j++)
            e_vecs[i][j]=temp[j][i];

      if (ascend_order)
        eigsort(e_vals,e_vecs,n);
      else
        eigsort(e_vals,e_vecs,(-1)*n);

      free(fv1);
      free_matrix(temp,n);

      }

} /* extern "C" */
