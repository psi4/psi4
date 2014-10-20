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
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
#include <libqt/qt.h>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>
#include "CAS_globaldefs.h"
#include "CAS_globals.h"
#include "psi4-dec.h"

namespace psi { namespace detcas {

void rotate_test(int dim, int npairs, int *p_arr, int *q_arr, 
                 double *theta_arr);


void rotate_orbs_irrep(int irrep, int dim, double **mo_coeffs,
                      int npairs, int *p_arr, int *q_arr, double *theta_arr)
{

  int p, q, i, j, pair;
  double **tmpmat, theta, sintheta, costheta;

  if (Params.print_lvl > 3) {
    outfile->Printf("Thetas for irrep %d\n", irrep);
    for (pair=0; pair<npairs; pair++) {
      outfile->Printf("Pair (%2d,%2d) = %12.6lf\n",
              p_arr[pair], q_arr[pair], theta_arr[pair]);
    }
    outfile->Printf("\n");
    //fflush(outfile);
  }

  /* apply the transformation C = C^0 U, using the fact that */
  /* C is written in terms of SOs, so it's in square blocks  */

  /* first copy the MO coefficients in to a block matrix     */
  tmpmat = block_matrix(dim,dim);
  for (i=0; i<dim; i++) {
    for (j=0; j<dim; j++) {
      tmpmat[i][j] = mo_coeffs[i][j];
    }
  }
 
  /* print new coefficients */
  if (Params.print_mos) {
    outfile->Printf("\n\tOld molecular orbitals for irrep %s\n", 
      CalcInfo.labels[irrep]);
    print_mat(tmpmat, dim, dim, "outfile");
  }

  /* now apply U, as a series of Givens rotation matrices */
  for (pair=0; pair<npairs; pair++) {
    p = *p_arr++; 
    q = *q_arr++;
    theta = *theta_arr++;
    theta = -theta; /* the DROT routine rotates around the other way    */
                    /* compared to our def of G(theta) in the VBD paper */
    costheta = cos(theta);
    sintheta = sin(theta);
    C_DROT(dim,&(tmpmat[0][q]),dim,&(tmpmat[0][p]),dim,costheta,sintheta);
  }

  /* print new coefficients */
  if (Params.print_mos) {
    outfile->Printf("\n\tNew molecular orbitals for irrep %s\n", 
      CalcInfo.labels[irrep]);
    print_mat(tmpmat, dim, dim, "outfile");
  }


  /* write the new block of MO coefficients to file30 */
  chkpt_init(PSIO_OPEN_OLD);
  chkpt_wt_scf_irrep(tmpmat, irrep);
  chkpt_close();

  free_block(tmpmat);
}


/*
** rotate_test
**
** This does about the same thing as rotate_orbs_irrep except that
** it does it on a unit matrix so the results can be checked easily
*/
void rotate_test(int dim, int npairs, int *p_arr, int *q_arr, 
                 double *theta_arr)
{

  int p, q, i, j, pair;
  double **tmpmat, theta, sintheta, costheta;

  /* set up a unit matrix */
  tmpmat = block_matrix(dim,dim);
  for (i=0; i<dim; i++) {
    tmpmat[i][i] = 1.0;
  }
 
  /* print new coefficients */
  outfile->Printf("\n\tOld molecular orbitals\n");
  print_mat(tmpmat, dim, dim, "outfile");


  /* now apply U, as a series of Givens rotation matrices */
  for (pair=0; pair<npairs; pair++) {
    p = *p_arr++; 
    q = *q_arr++;
    theta = *theta_arr++;
    theta = -theta;  /* DROT rotates around the other way */
    costheta = cos(theta);
    sintheta = sin(theta);
    outfile->Printf("\nApplying rotation (%2d,%2d) = %12.6lf\n", p, q, theta);
    outfile->Printf("Cos(theta)=%12.6lf, Sin(theta)=%12.6lf\n", 
            costheta, sintheta);
    C_DROT(dim,&(tmpmat[0][q]),dim,&(tmpmat[0][p]),dim,costheta,sintheta);
    outfile->Printf("\n\tMatrix after transformation:\n");
    print_mat(tmpmat, dim, dim, "outfile");
  }

  /* print new coefficients */
  outfile->Printf("\n\tNew molecular orbitals\n");
  print_mat(tmpmat, dim, dim, "outfile");

  free_block(tmpmat);

}

}} // end namespace psi::detcas

