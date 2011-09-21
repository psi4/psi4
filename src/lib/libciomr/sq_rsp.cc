/*!
** \file
** \brief Diagnoalize a symmetrix square matrix
** \ingroup CIOMR
*/

#include <psifiles.h>
#include "libciomr.h"
#include <cstdlib>
#include "libqt/qt.h"

namespace psi {

extern void tred2(int n,double** a,double* d,double* e,int matz);
extern void tqli(int n, double *d, double **z, double *e, int matz,
  double toler);

/* translation into c of a translation into FORTRAN77 of the EISPACK */
/* matrix diagonalization routines */

  /**
  *  WARNING: Psi 3 Fortran routine sq_rsp  deprecated
  *  by Robert Parrish, robparrish@gmail.com
  *
  *  sq_rsp now calls the LAPACK method DSYEV for
  *  numerical stability, speed, and threading
  *
  *  the signature of this method remains the same
  *
  *  June 22, 2010
  **/

/*!
** sq_rsp(): diagomalize a symmetric square matrix ('array').
**
** \param nm     = rows of matrix
** \param n      = columns of matrix
** \param nv     = number of elements in lower triangle (n*(n+1)/2)
** \param array  = matrix to diagonalize
** \param e_vals = array to hold eigenvalues
** \param matz   = 0 (no eigenvectors, eigenvals in ascending order)
**               = 1 (eigenvectors and eigenvalues in ascending order)
**               = 2 (no eigenvectors, eigenvalues in descending order)
**               = 3 (eigenvectors and eigenvalues in descending order)
** \param e_vecs = matrix of eigenvectors (one column for each eigvector)
** \param toler  = tolerance for eigenvalues?  Often 1.0E-14.
**
** \ingroup CIOMR
*/
void sq_rsp(int /*nm*/, int n, double **array, double *e_vals, int matz,
  double **e_vecs, double /*toler*/)
{

  if ((matz > 3) || (matz < 0)) matz = 0;
  //Do you want eigenvectors?
  bool eigenvectors = (matz == 1 || matz == 3)?true:false;
  //Ascending or Descending?
  bool ascending = (matz == 0 || matz == 1)?true:false;

  //Get Eigenvectors (Use 'V')
  if (eigenvectors) {

    //First Temp array (the lda isn't so fly in libciomr)
    double** Temp_sqrsp = block_matrix(n,n);

    //Copy array to Temp_sqrsp (for loops required)
    for (int r = 0; r<n; r++)
        for (int c = 0; c<n; c++)
            Temp_sqrsp[r][c] = array[r][c];

    //fprintf("  Initial Matrix:\n");
    //print_mat(e_vecs,n,n,outfile);
    //printf("%d",n);

    //Get scratch vector and call DSYEV
    //The eigenvectors are placed in e_vecs in ascending order
    int lwork_sqrsp = 3*n;
    double* work_sqrsp = init_array(lwork_sqrsp);
    C_DSYEV('V','U',n,&Temp_sqrsp[0][0],n,&e_vals[0],&work_sqrsp[0],lwork_sqrsp);
    free(work_sqrsp);

    //fprintf(outfile,"  Eigenvectors:\n");
    //print_mat(e_vecs,n,n,outfile);


    //printf("\n  Eigenvalues:\n");
    //for (int r = 0; r<n; r++)
    //    printf("  r = %d, %14.10f\n",r+1,e_vals[r]);

    //LAPACK stores eigenvectors in rows, we need them in columns
    //This overhead is why you should always call DSYEV!
    double **T_sqrsp_2 = block_matrix(n,n);
    C_DCOPY(n*n,&Temp_sqrsp[0][0],1,&T_sqrsp_2[0][0],1);
    for (int r = 0; r<n; r++)
        C_DCOPY(n,&T_sqrsp_2[r][0],1,&Temp_sqrsp[0][r],n);
    free_block(T_sqrsp_2);

    //If descending is required, the canonical order must be reversed
    //Sort is stable
    if (!ascending) {
        double* Temp_sqrsp_col = init_array(n);
        double w_Temp_sqrsp;

        for (int c = 0; c<n/2; c++) {

            //Swap eigenvectors
            C_DCOPY(n,&Temp_sqrsp[0][c],n,&Temp_sqrsp_col[0],1);
            C_DCOPY(n,&Temp_sqrsp[0][n-c-1],n,&Temp_sqrsp[0][c],n);
            C_DCOPY(n,&Temp_sqrsp_col[0],1,&Temp_sqrsp[0][n-c-1],n);

            //Swap eigenvalues
            w_Temp_sqrsp = e_vals[c];
            e_vals[c] = e_vals[n-c-1];
            e_vals[n-c-1] = w_Temp_sqrsp;

        }

        free(Temp_sqrsp_col);
    }
    //Copy from Temp_sqrsp to e_vecs (for loops required)
    for (int r = 0; r<n; r++)
        for (int c = 0; c<n; c++)
            e_vecs[r][c] = Temp_sqrsp[r][c];

    free_block(Temp_sqrsp);
    //fprintf(outfile,"  Eigenvectors (After sort):\n");
    //print_mat(e_vecs,n,n,outfile);

    //printf("\n  Eigenvalues (After sort):\n");
    //for (int r = 0; r<n; r++)
    //    printf("  r = %d, %14.10f\n",r+1,e_vals[r]);

    //No Eigenvectors (Use 'N')
  } else {
    //e_vecs might not be initialized
    //Use a Temp array
    double** Temp_sqrsp = block_matrix(n,n);
    //Copy array to Temp
    //must use for loops because of bloody irrep blocking in CSCF
    for (int r = 0; r<n; r++)
        for (int c = 0; c<n; c++)
            Temp_sqrsp[r][c] = array[r][c];

    //Form scratch array and call DSYEV
    //'N' in parameter 1 to DSYEV terminates
    //the algorithm after eigenvectors are found
    //Canonical order is ascending in LAPACK
    int lwork_sqrsp = 3*n;
    double* work_sqrsp = init_array(lwork_sqrsp);
    C_DSYEV('N','U',n,&Temp_sqrsp[0][0],n,&e_vals[0],&work_sqrsp[0],lwork_sqrsp);
    free(work_sqrsp);
    free_block(Temp_sqrsp);

    //If descending is required, the canonical order must be reversed
    //Sort is stable
    if (!ascending) {
        double w_Temp_sqrsp;

        //Eigenvalues only
        for (int c = 0; c<n/2; c++) {

            w_Temp_sqrsp = e_vals[c];
            e_vals[c] = e_vals[n-c-1];
            e_vals[n-c-1] = w_Temp_sqrsp;

        }

    }
  }
  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  //
  //  DEPRECATED METHOD AND ASSOCIATED CALLS
  //
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
/**
  int i, j, ierr;
  int ascend_order;
  double *fv1, **temp;

  // Modified by Ed - matz can have the values 0 through 3

  if ((matz > 3) || (matz < 0)) {
    matz = 0;
    ascend_order = 1;
  }
  else
    if (matz < 2)
      ascend_order = 1;	// Eigenvalues in ascending order
  else {
    matz -= 2;
    ascend_order = 0;	// Eigenvalues in descending order
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
**/
}

}

