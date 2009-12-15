/*! \file
    \ingroup DETCAS
    \brief Enter brief description of file here 
*/
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
#include <libqt/qt.h>
#include <cmath>
#include "globaldefs.h"
#include "globals.h"

namespace psi { namespace detcas {

#define MO_HESS_MIN 1.0E-2


/*
** calc_orb_step()
**
** This function calculates the step in theta space for the orbitals
** given the orbital gradient and an approximate orbital Hessian
**
** C. David Sherrill
** April 1998
*/
void calc_orb_step(int npairs, double *grad, double *hess_diag, double *theta)
{

  int pair;
  double numer, denom;

  for (pair=0; pair<npairs; pair++) {
    numer = grad[pair];
    denom = hess_diag[pair];
    if (denom < 0.0) {
      fprintf(outfile, "Warning: MO Hessian denominator negative\n");
      denom = -denom;
    }
    if (denom < MO_HESS_MIN) {
      fprintf(outfile, "Warning: MO Hessian denominator too small\n");
      denom = MO_HESS_MIN;
    } 
    theta[pair] =  - numer / denom;
  }

}


/*
** calc_orb_step_full()
**
** This function calculates the step in theta space for the orbitals
** given the orbital gradient and a square matrix orbital Hessian
**
** C. David Sherrill
** September 2003
*/
void calc_orb_step_full(int npairs, double *grad, double **hess, double *theta)
{
  double **hess_inv;
  double **hess_copy; /* for testing! */
  int i,j;
  double tval;
  int solved;
  double *BVector;
  int *pivots;
  double hess_det = 1.0;
  int *indx;
  double biggest_step;

  hess_copy = block_matrix(npairs, npairs);
  indx = init_int_array(npairs);
 
  for (i=0; i<npairs; i++) {
    for (j=0; j<npairs; j++) {
      hess_copy[i][j] = hess[i][j];
    }
  }

  ludcmp(hess_copy,npairs,indx,&hess_det);
  for (j=0;j<npairs;j++){
    hess_det *= hess_copy[j][j];
  }
  fprintf(outfile,"The determinant of the hessian is %8.3E\n",hess_det);
  fflush(outfile);

  /* 
     if the orbital Hessian is not positive definite, we may have some
     trouble converging the orbitals.  Guarantee it's positive definite
     by levelshifting
  */
  if (Params.level_shift) {
    while (hess_det < Params.determ_min) {
      fprintf(outfile,"Level shifting the hessian by %8.3E\n",Params.shift);
      for (i=0;i<npairs;i++) {
        hess[i][i] += Params.shift;
      }
      for (i=0;i<npairs;i++) {
        for (j=0;j<npairs;j++) {
          hess_copy[i][j]=hess[i][j];
        }
      }
      ludcmp(hess_copy,npairs,indx,&hess_det);
      for (j=0;j<npairs;j++){
        hess_det *= hess_copy[j][j];
      }
      fprintf(outfile,"The determinant of the hessian is %8.3E\n",hess_det);
    }
    fprintf(outfile,"Determinant of the hessian is greater than %8.3E\n",
      Params.determ_min);
  }


  /* let's re-copy hess into hess_copy because we ludcmp'd hess_copy */
  for (i=0;i<npairs;i++) {
    for (j=0;j<npairs;j++) {
      hess_copy[i][j]=hess[i][j];
    }
  }
 
  if (!Params.invert_hessian) { /* solve H delta = - g */
    fprintf(outfile,"Solving system of linear equations for orbital step...");
    BVector = init_array(npairs);
    pivots = init_int_array((npairs * (npairs - 1))/2);
    for(i=0;i<npairs;i++){
      BVector[i] = -grad[i];
      theta[i] = 0.0;
    }
    solved = C_DGESV(npairs,1,&(hess_copy[0][0]),npairs,pivots,
      BVector,npairs);
    if (solved == 0) {
      fprintf(outfile,"equations solved!\n");
      for(i=0;i<npairs;i++) {
        theta[i] = BVector[i];
      }
    }
    else {
      fprintf(outfile,"FAILED TO SOLVE FOR THETA VALUES\n");
      fprintf(outfile,"DGESV returned error %5d \n",solved);
      exit(PSI_RETURN_FAILURE);
    }
    free(BVector);
    free(pivots);
  } /* end solution of linear equations H delta = -g */

  else { /* direct inversion of orbital Hessian */
    fprintf(outfile,"Attempting to directly invert the Hessian matrix\n");
    hess_inv = block_matrix(npairs,npairs);

    /* note: this will destroy hessian matrix; don't use it again later! */
    invert_matrix(hess_copy,hess_inv,npairs,outfile);

    /* debug check */
    mmult(hess_inv,0,hess,0,hess_copy,0,npairs,npairs,npairs,0);
    fprintf(outfile, "Hessian * Hessian inverse = \n");
    print_mat(hess_copy,npairs,npairs,outfile); 
    fprintf(outfile, "\n");
  
    /* step = - B^{-1} * g */
    zero_arr(theta,npairs);
    /* the line below seems to have trouble unless I take out the -1
       which should be there, and even then it's not really working */
    /*
    C_DGEMV('n',npairs,npairs,-1.0,hess_inv[0],npairs,grad,1,0.0,theta,1);
    */

    for (i=0; i<npairs; i++) {
      tval = 0.0;
      for (j=0; j<npairs; j++) {
        tval += hess_inv[i][j] * grad[j];
      }
      theta[i] = - tval;
    }
    free_block(hess_inv);
  } /* end direct inversion of Hessian */

  /* make sure the step is not too big */
  biggest_step = 0.0;
  for (i=0; i<npairs; i++) {
    tval = theta[i];
    if (fabs(tval) > biggest_step) biggest_step = fabs(tval);
  }
  fprintf(outfile,"\nLargest step in theta space is %12.6lf \n", biggest_step);
  if (biggest_step > Params.step_max) {
    fprintf(outfile, "Scaling the step\n");
    for (i=0;i<npairs;i++) {
      theta[i] = theta[i] * Params.step_max / biggest_step;
    }
  }

  free_block(hess_copy);
  free(indx);
}


/*
** calc_orb_step_bfgs()
**
** This function calculates the step in theta space for the orbitals
** given the orbital gradient and a square matrix orbital Hessian INVERSE.
** With the inverse already available, this is very straightforward.
**
** C. David Sherrill
** March 2004
*/
void calc_orb_step_bfgs(int npairs, double *grad, double **hess, double *theta)
{

  int i, j;
  double tval, biggest_step;

  for (i=0; i<npairs; i++) {
    tval = 0.0;
    for (j=0; j<npairs; j++) {
      tval += hess[i][j] * grad[j];
    }
    theta[i] = - tval;
  }

  /* make sure the step is not too big */
  biggest_step = 0.0;
  for (i=0; i<npairs; i++) {
    tval = theta[i];
    if (fabs(tval) > biggest_step) biggest_step = fabs(tval);
  }
  fprintf(outfile,"\nLargest step in theta space is %12.6lf \n", biggest_step);
  if (biggest_step > Params.step_max) {
    fprintf(outfile, "Largest allowed step %12.6lf --- scaling the step\n",
      Params.step_max);
    for (i=0;i<npairs;i++) {
      theta[i] = theta[i] * Params.step_max / biggest_step;
    }
  }

}


/*
** print_step
**
** This function prints out the information for a given orbital iteration
*/
void print_step(int npairs, int steptype)
{
  FILE *sumfile;
  char sumfile_name[] = "file14.dat";
  int i, entries, iter, *nind;
  double *rmsgrad, *scaled_rmsgrad, *energies, energy;
  char **comments;

  /* open ascii file, get number of entries already in it */

  ffile_noexit(&sumfile,sumfile_name,2);
  if (sumfile == NULL) { /* the file doesn't exist yet */
    entries = 0;
    if (Params.print_lvl)
      fprintf(outfile, "\nPreparing new file %s\n", sumfile_name);
  }
  else {
    if (fscanf(sumfile, "%d", &entries) != 1) {
      fprintf(outfile,"(print_step): Trouble reading num entries in file %s\n",
        sumfile_name);
      fclose(sumfile);
      return;
    }
  }

  rmsgrad = init_array(entries+1);
  scaled_rmsgrad = init_array(entries+1);
  energies= init_array(entries+1);
  nind = init_int_array(entries+1);
  comments = (char **) malloc ((entries+1) * sizeof (char *));
  for (i=0; i<entries+1; i++) {
    comments[i] = (char *) malloc (MAX_COMMENT * sizeof(char));
  }

  for (i=0; i<entries; i++) {
    fscanf(sumfile, "%d %d %lf %lf %lf %s", &iter, &(nind[i]), 
           &(scaled_rmsgrad[i]), &(rmsgrad[i]), &(energies[i]), comments[i]);
  }

  chkpt_init(PSIO_OPEN_OLD);
  if (chkpt_exist("State averaged energy")) {
    energy = chkpt_rd_e_labeled("State averged energy");
  }
  else
    energy = chkpt_rd_etot();
  chkpt_close();

  scaled_rmsgrad[entries] = CalcInfo.scaled_mo_grad_rms;
  rmsgrad[entries] = CalcInfo.mo_grad_rms;
  energies[entries] = energy;
  nind[entries] = npairs;

  if (steptype == 0) 
    strcpy(comments[entries], "CONV");
  else if (steptype == 1)
    strcpy(comments[entries], "NR");
  else if (steptype == 2)
    strcpy(comments[entries], "DIIS"); 
  else {
    fprintf(outfile, "(print_step): Unrecognized steptype %d\n", steptype);
    strcpy(comments[entries], "?");
  }

  if (entries) fclose(sumfile);

  /* now open file for writing, write out old info plus new */
  ffile_noexit(&sumfile,"file14.dat",0);
  if (sumfile == NULL) {
    fprintf(outfile, "(print_step): Unable to open file %s\n", sumfile_name);
  }
  else {
    entries++;
    fprintf(sumfile, "%5d\n", entries);
    for (i=0; i<entries; i++) {
      fprintf(sumfile, "%5d %5d %14.9lf %14.9lf %20.12lf %9s\n", i+1, nind[i], 
              scaled_rmsgrad[i], rmsgrad[i], energies[i], comments[i]);
    }
    fclose(sumfile);
  }

  free(scaled_rmsgrad);
  free(rmsgrad);
  free(energies);
  free(nind);
  for (i=0; i<entries; i++)
    free(comments[i]);
  free(comments);

}

}} // end namespace psi::detcas

