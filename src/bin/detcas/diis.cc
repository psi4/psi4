/*! \file
    \ingroup DETCAS
    \brief Enter brief description of file here 
*/
/*
** DIIS.C
** 
** Routines for Direct Inversion of the Iterative Subspace Interpolation
** of orbital rotation angles, P. Pulay, Chem. Phys. Lett. 73, 393 (1980).
**
** C. David Sherrill
** University of California, Berkeley
** May 1998
*/

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <libciomr/libciomr.h>
#include "globaldefs.h"
#include "globals.h"

namespace psi { namespace detcas {

#define DIIS_MIN_DET 1.0E-16

/*
** diis()
**
** This top-level routine manages all the DIIS stuff.  
**
** Parameters:
**   veclen = length of the vectors to extrapolate
**   vec    = new vector to add
**   errvec = new error vector to add
**
** Returns:
**   1 if DIIS step taken, otherwise 0
*/
int diis(int veclen, double *vec, double *errvec)
{
  int i, j, k;
  int num_vecs, new_num_vecs, offset, diis_iter, do_diis;
  double **vecs, **errvecs, **bmat, *bvec, tval, det, scale_factor;
  FILE *fp;

  
  /* add the vector and error vector to subspace */
  ffileb_noexit(&fp,"diis.dat",2);
  if (fp != NULL) {

    if (fread(&num_vecs, sizeof(int), 1, fp) != 1) {
      fprintf(outfile, "(diis): Error reading number of diis vectors.\n");
      return(0);
    }

    if (fread(&diis_iter, sizeof(int), 1, fp) != 1) {
      fprintf(outfile, "(diis): Error reading diis iteration number.\n");
      return(0);
    }

    vecs = block_matrix(num_vecs+1, veclen);
    errvecs = block_matrix(num_vecs+1, veclen);

    for (i=0; i<num_vecs; i++) {

      if (fread(vecs[i], sizeof(double), veclen, fp) != veclen) { 
        fprintf(outfile, "(diis): Error reading diis vector %d\n", i);
        return(0);
      }

      if (fread(errvecs[i], sizeof(double), veclen, fp) != veclen) { 
        fprintf(outfile, "(diis): Error reading diis error vector %d\n", i);
        return(0);
      }

    }

    fclose(fp);
  }

  else {
    num_vecs = 0;
    diis_iter = 0;
    vecs = block_matrix(num_vecs+1, veclen);
    errvecs = block_matrix(num_vecs+1, veclen);
  } 

  /* will we do a diis this time? */
  diis_iter++;
  num_vecs++;
  if ((diis_iter % Params.diis_freq) || (num_vecs < Params.diis_min_vecs)) 
    do_diis = 0; 
  else
    do_diis = 1;

  for (i=0; i<veclen; i++) vecs[num_vecs-1][i] = vec[i];
  for (i=0; i<veclen; i++) errvecs[num_vecs-1][i] = errvec[i];

  offset = 0;
  if (num_vecs > Params.diis_max_vecs) 
    offset = num_vecs - Params.diis_max_vecs;

  if (Params.print_lvl > 2) 
    fprintf(outfile, "Diis: iter %2d, vecs %d, do_diis %d, offset %d\n", 
            diis_iter, num_vecs, do_diis, offset);

  new_num_vecs = num_vecs - offset;

  /* write out the diis info */
  ffileb_noexit(&fp,"diis.dat",0);
  if (fp == NULL) {
    fprintf(outfile, "(diis): Error opening diis.dat\n");
    return(0);
  } 

  if (fwrite(&new_num_vecs, sizeof(int), 1, fp) != 1) {
    fprintf(outfile, "(diis): Error writing number of diis vectors.\n");
    return(0);
  }

  if (fwrite(&diis_iter, sizeof(int), 1, fp) != 1) {
    fprintf(outfile, "(diis): Error writing diis iteration number.\n");
    return(0);
  }

  for (i=offset; i<num_vecs; i++) {
    if (fwrite(vecs[i], sizeof(double), veclen, fp) != veclen) {
      fprintf(outfile, "(diis): Error writing vector %d.\n", i);
    }
    if (fwrite(errvecs[i], sizeof(double), veclen, fp) != veclen) {
      fprintf(outfile, "(diis): Error writing error vector %d.\n", i);
    }
  }

  fclose(fp);


  /* don't take a diis step if it's not time */
  if (!do_diis) {
    free_block(vecs);
    free_block(errvecs);
    return(0); 
  }

  /* form diis matrix, solve equations */
  if (Params.print_lvl) 
    fprintf(outfile, "Attempting a DIIS step\n");

  bmat = block_matrix(new_num_vecs+1, new_num_vecs+1);
  bvec = init_array(new_num_vecs+1);

  for (i=0; i<new_num_vecs; i++) {
    for (j=0; j<=i; j++) {
      tval = 0.0;
      for (k=0; k<veclen; k++) {
        tval += errvecs[i+offset][k] * errvecs[j+offset][k];
      }
      bmat[i][j] = tval;
      bmat[j][i] = tval;
    }
  }

  for (i=0; i<new_num_vecs; i++) {
    bmat[i][new_num_vecs] = -1.0;
    bmat[new_num_vecs][i] = -1.0;
  }

  bmat[new_num_vecs][new_num_vecs] = 0.0;
  bvec[new_num_vecs] = -1.0;

  if (Params.print_lvl > 2) {
    fprintf(outfile, "DIIS B Matrix:\n");
    print_mat(bmat, new_num_vecs+1, new_num_vecs+1, outfile);
  }

  /* scale the B matrix */
  scale_factor = 1.0 / bmat[0][0];
  for (i=0; i<new_num_vecs; i++) {
    for (j=0; j<new_num_vecs; j++) {
      bmat[i][j] = bmat[i][j] * scale_factor;
    }
  }

  if (Params.print_lvl > 2) {
    fprintf(outfile, "DIIS B Matrix:\n");
    print_mat(bmat, new_num_vecs+1, new_num_vecs+1, outfile);
  }

  /* now solve the linear equations */ 
  flin(bmat, bvec, new_num_vecs+1, 1, &det); 

  if (fabs(det) < DIIS_MIN_DET) {
    fprintf(outfile, "Warning: diis matrix near-singular\n");
    fprintf(outfile, "Determinant is %6.3E\n", det);
  }

  if (Params.print_lvl > 3) {
    fprintf(outfile, "\nCoefficients of DIIS extrapolant:\n");
    for (i=0; i<new_num_vecs; i++) {
      fprintf(outfile, "%12.6lf\n", bvec[i]);
    }
    fprintf(outfile, "Lambda:\n");
    fprintf(outfile, "%12.6lf\n", bvec[i]);
  }

  /* get extrapolated vector */
  zero_arr(CalcInfo.theta_cur, veclen);
  for (i=0; i<new_num_vecs; i++) {
    for (j=0; j<veclen; j++) {
      CalcInfo.theta_cur[j] += bvec[i] * vecs[i+offset][j];
    }
  }

  free_block(vecs);
  free_block(errvecs);
  free_block(bmat); 
  free(bvec);

  return(1);
}

}} // end namespace psi::detcas

