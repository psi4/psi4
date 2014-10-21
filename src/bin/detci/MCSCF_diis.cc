/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

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
#include <psi4-dec.h>
#include <libpsio/psio.hpp>
#include <libpsio/psio.h>
#include <psifiles.h>
#include "MCSCF_globaldefs.h"
#include "MCSCF_globals.h"

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
  char diis_char[80];
  
  /* add the vector and error vector to subspace */
  if (psio_tocentry_exists(PSIF_DETCAS, "Num vectors")){ 
    psio_open(PSIF_DETCAS, PSIO_OPEN_OLD);

    psio_read_entry(PSIF_DETCAS, "Num vectors", (char *) &(num_vecs),
                    sizeof(int));
    psio_read_entry(PSIF_DETCAS, "Iteration number", (char *) &(diis_iter),
                    sizeof(int));

    vecs = block_matrix(num_vecs+1, veclen);
    errvecs = block_matrix(num_vecs+1, veclen);

    for (i=0; i<num_vecs; i++) {
      sprintf(diis_char, "DIIS vector %3d", i);
      psio_read_entry(PSIF_DETCAS, diis_char, (char *) vecs[i],
                      veclen*sizeof(double));

      sprintf(diis_char, "DIIS error vector %3d", i);
      psio_read_entry(PSIF_DETCAS, diis_char, (char *) errvecs[i],
                      veclen*sizeof(double));

    }

    psio_close(PSIF_DETCAS, 1);
  }
  else {
    num_vecs = 0;
    diis_iter = 0;
    vecs = block_matrix(num_vecs+1, veclen);
    errvecs = block_matrix(num_vecs+1, veclen);
  } 
  /* end diis read */

  /* will we do a diis this time? */
  diis_iter++;
  num_vecs++;
  if ((diis_iter % Parameters.diis_freq) || (num_vecs < Parameters.diis_min_vecs)) 
    do_diis = 0; 
  else
    do_diis = 1;

  for (i=0; i<veclen; i++) vecs[num_vecs-1][i] = vec[i];
  for (i=0; i<veclen; i++) errvecs[num_vecs-1][i] = errvec[i];

  offset = 0;
  if (num_vecs > Parameters.diis_max_vecs) 
    offset = num_vecs - Parameters.diis_max_vecs;

  if (Parameters.print_lvl > 2) 
    outfile->Printf("Diis: iter %2d, vecs %d, do_diis %d, offset %d\n", 
            diis_iter, num_vecs, do_diis, offset);

  new_num_vecs = num_vecs - offset;

  /* write out the diis info */
  psio_open(PSIF_DETCAS, PSIO_OPEN_OLD);

  psio_write_entry(PSIF_DETCAS, "Num vectors", (char *) &(new_num_vecs),
                  sizeof(int));
  psio_write_entry(PSIF_DETCAS, "Iteration number", (char *) &(diis_iter),
                  sizeof(int));

  for (i=offset; i<num_vecs; i++) {
    sprintf(diis_char, "DIIS vector %3d", i);
    psio_write_entry(PSIF_DETCAS, diis_char, (char *) vecs[i],
                    veclen*sizeof(double));

    sprintf(diis_char, "DIIS error vector %3d", i);
    psio_write_entry(PSIF_DETCAS, diis_char, (char *) errvecs[i],
                    veclen*sizeof(double));

  }

  psio_close(PSIF_DETCAS, 1);
  /* end write out*/  

  /* don't take a diis step if it's not time */
  if (!do_diis) {
    free_block(vecs);
    free_block(errvecs);
    return(0); 
  }

  /* form diis matrix, solve equations */
  if (Parameters.print_lvl) 
    outfile->Printf("Attempting a DIIS step\n");

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

  if (Parameters.print_lvl > 2) {
    outfile->Printf("DIIS B Matrix:\n");
    print_mat(bmat, new_num_vecs+1, new_num_vecs+1, "outfile");
  }

  /* scale the B matrix */
  scale_factor = 1.0 / bmat[0][0];
  for (i=0; i<new_num_vecs; i++) {
    for (j=0; j<new_num_vecs; j++) {
      bmat[i][j] = bmat[i][j] * scale_factor;
    }
  }

  if (Parameters.print_lvl > 2) {
    outfile->Printf("DIIS B Matrix:\n");
    print_mat(bmat, new_num_vecs+1, new_num_vecs+1, "outfile");
  }

  /* now solve the linear equations */ 
  flin(bmat, bvec, new_num_vecs+1, 1, &det); 

  if (fabs(det) < DIIS_MIN_DET) {
    outfile->Printf("Warning: diis matrix near-singular\n");
    outfile->Printf("Determinant is %6.3E\n", det);
  }

  if (Parameters.print_lvl > 3) {
    outfile->Printf("\nCoefficients of DIIS extrapolant:\n");
    for (i=0; i<new_num_vecs; i++) {
      outfile->Printf("%12.6lf\n", bvec[i]);
    }
    outfile->Printf("Lambda:\n");
    outfile->Printf("%12.6lf\n", bvec[i]);
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

