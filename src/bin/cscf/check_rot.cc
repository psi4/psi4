/*! \file 
    \ingroup (CSCF)
    \brief Enter brief description of file here 
*/

/*! \defgroup CSCF Add a description of the group CSCF */

/*
** CHECK_ROT.C
** 
** Check the MO rotation performed in rotate_vector() to make sure that
** columns of the C matrix haven't swapped.  If so, swap them (and their
** eigenvalues) back to the correct ordering.
**
** David Sherrill and Daniel Crawford, July 1995
**
*/

#include <stdio.h>
#include <libciomr/libciomr.h>
#define EXTERN
#include "includes.h"
#include "common.h"

namespace psi { namespace cscf {

void swap_vectors(double **c, int nn, int j, int i);

void check_rot(int nn, int num_mo, double **cold, double **cnew, double *smat_pac, 
      double *fock_evals, int irrep)
{
  int i,j,jmaxcoeff,swapped;
  double maxcoeff=0.0, tval;
  double **smat;
  double **tmp1, **tmp2;
  static int printed=0;
  int count=0;

  smat = init_matrix(nn,nn);
  tmp1 = init_matrix(nn,nn);
  tmp2 = init_matrix(nn,nn);

  tri_to_sq(smat_pac,smat,nn);

  do {
      swapped = 0;
/*      mxmb(cnew,nn,1,smat,1,nn,tmp1,1,nn,nn,nn,nn);
      mxmb(tmp1,1,nn,cold,1,nn,tmp2,1,nn,nn,nn,nn);*/
      mmult(cnew,1,smat,0,tmp1,0,num_mo,nn,nn,0);
      mmult(tmp1,0,cold,0,tmp2,0,num_mo,nn,num_mo,0);

/*       fprintf(outfile, "C_New Matrix:\n");
         print_mat(cnew,nn,nn,outfile);
         fprintf(outfile, "C_New * S * C_old:\n");
         print_mat(tmp2,nn,nn,outfile); */
         

      for(i=0; i < num_mo; i++) {
        maxcoeff = 0.0;
        for(j=0; j < num_mo; j++) {
          if(fabs(tmp2[j][i]) > maxcoeff) { 
            maxcoeff = fabs(tmp2[j][i]);
            jmaxcoeff = j;
          }
        }
        if (maxcoeff < 0.75) {
           fprintf(outfile, "\n   Warning!  Diagonality check C'SC gives ");
           fprintf(outfile, "a maximum element of\n   %lf for (%d,%d)\n",
                   maxcoeff, jmaxcoeff+1, i+1);
           fprintf(outfile, "   Won't perform MO swapping for this column\n");
           }
        else if (jmaxcoeff != i) {
            swap_vectors(cnew,nn,jmaxcoeff,i);
            tval = fock_evals[jmaxcoeff];
            fock_evals[jmaxcoeff] = fock_evals[i];
            fock_evals[i] = tval;
            swapped = 1;  count++;
            if (!printed) {
               fprintf(outfile, 
                  "\n   Warning!  MO rotation swapped columns of C matrix\n");
               printed = 1;
               }
            fprintf(outfile, 
               "   Swapping back columns %d and %d for irrep %d\n",
               jmaxcoeff+1,i+1,irrep+1);
            break;
         }
      }
  } while(swapped && count < 50);

  if (count == 50) {
     fprintf(outfile, "(check_rot): This is bad.  Tried 50 swaps for ");
     fprintf(outfile, "irrep %d!\n", irrep+1);
     fprintf(outfile, "   You may want to set check_rot = false\n");
     }

  free_matrix(smat,nn);
  free_matrix(tmp1,nn);
  free_matrix(tmp2,nn);
  
  /*  
  fprintf(outfile, "C_New Matrix(after phase change):\n");
  print_mat(cnew,nn,nn,outfile); 
  */
}


void swap_vectors(double **c, int nn, int j, int i)
{
  int k;
  double tmp;

  for(k=0; k < nn; k++)  {
      tmp = c[k][j];
      c[k][j] = c[k][i];
      c[k][i] = tmp;
    }
}

}} // namespace psi::cscf
