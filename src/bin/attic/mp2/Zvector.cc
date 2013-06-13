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
    \ingroup MP2
    \brief Enter brief description of file here 
*/
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include <libciomr/libciomr.h>
#include <cmath>
#define EXTERN
#include "globals.h"

namespace psi{ namespace mp2{

void rhf_Zvector(void);
void uhf_Zvector(void);
void rhf_sf_Zvector(void);
void uhf_sf_Zvector(void);

void Zvector(void)
{
  if(params.gradient) {
    if(params.ref == 0) rhf_sf_Zvector();
    else if(params.ref == 2) uhf_sf_Zvector();
  }
  else {
    if(params.ref == 0) rhf_Zvector();
    else if(params.ref == 2) uhf_Zvector();
  }
}

void rhf_Zvector(void)
{
  dpdfile2 L;
  dpdfile2 D;
  dpdbuf4 A;
  double **Z;
  int h, nirreps;
  int a, i, num_ai, count;
  int I, B;

  nirreps = mo.nirreps;

  dpd_->file2_init(&L, PSIF_CC_OEI, 0, 1, 0, "LAI");
  dpd_->file2_mat_init(&L);
  dpd_->file2_mat_rd(&L);
  num_ai = 0;
  for(h=0; h < nirreps; h++)
    num_ai += L.params->rowtot[h]*L.params->coltot[h];

  Z = block_matrix(1,num_ai);
  for(h=0,count=0; h < nirreps; h++)
    for(a=0; a < L.params->rowtot[h]; a++)
      for(i=0; i < L.params->coltot[h]; i++) 
	Z[0][count++] = -L.matrix[h][a][i];

  dpd_->file2_mat_close(&L);
  dpd_->file2_close(&L);

  dpd_->buf4_init(&A, PSIF_CC_MISC, 0, 11, 11, 11, 11, 0, "A(AI,BJ)");
  dpd_->buf4_mat_irrep_init(&A, 0);
  dpd_->buf4_mat_irrep_rd(&A, 0);

  pople(A.matrix[0], Z[0], num_ai, 1, 1e-12, outfile, 0);

  dpd_->buf4_mat_irrep_close(&A, 0);
  dpd_->buf4_close(&A);

  dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 0, "DAI");
  dpd_->file2_mat_init(&D);
  for(h=0,count=0; h < nirreps; h++)
    for(a=0; a < D.params->rowtot[h]; a++)
      for(i=0; i < D.params->coltot[h]; i++) 
	D.matrix[h][a][i] = Z[0][count++];

  dpd_->file2_mat_wrt(&D);
  dpd_->file2_mat_close(&D);
  dpd_->file2_close(&D);

  free_block(Z);
}

void uhf_Zvector(void)
{

}

void rhf_sf_Zvector(void)
{
  dpdbuf4 A;
  dpdfile2 X1, D;
  double **Z;
  int num_ai, h, nirreps, a, i, count;

  nirreps = mo.nirreps;

  /* Place all the elements of the orbital rotation gradient, X into a
     linear array, Z */
  dpd_->file2_init(&X1, PSIF_CC_MISC, 0, 1, 0, "X(A,I)");
  dpd_->file2_mat_init(&X1);
  dpd_->file2_mat_rd(&X1);
  num_ai = 0;
  for(h=0; h < nirreps; h++)
    num_ai += X1.params->rowtot[h]*X1.params->coltot[h];

  Z = block_matrix(1,num_ai);
  for(h=0,count=0; h < nirreps; h++)
    for(a=0; a < X1.params->rowtot[h]; a++)
      for(i=0; i < X1.params->coltot[h]; i++) 
	Z[0][count++] = -X1.matrix[h][a][i];

  dpd_->file2_mat_close(&X1);
  dpd_->file2_close(&X1);

  /* Now, grab only irrep 0 of the orbital Hessian */
  dpd_->buf4_init(&A, PSIF_CC_MISC, 0, 11, 11, 11, 11, 0, "A(EM,AI)");
  dpd_->buf4_mat_irrep_init(&A, 0);
  dpd_->buf4_mat_irrep_rd(&A, 0);

  /* Trying out Matt's Pople code --- way to go, Matt! */
  pople(A.matrix[0], Z[0], num_ai, 1, 1e-12, outfile, 0);

  dpd_->buf4_mat_irrep_close(&A, 0);
  dpd_->buf4_close(&A);

  /* Build the orbital component of Dai --- we'll build these as separate
     spin cases for future simplicity (e.g., UHF-based codes)*/

  dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 0, "D(orb)(A,I)");
  dpd_->file2_mat_init(&D);
  for(h=0,count=0; h < nirreps; h++)
    for(a=0; a < D.params->rowtot[h]; a++)
      for(i=0; i < D.params->coltot[h]; i++) 
	D.matrix[h][a][i] = Z[0][count++];
  dpd_->file2_mat_wrt(&D);
  dpd_->file2_mat_close(&D);
  dpd_->file2_close(&D);

  dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 0, "D(orb)(a,i)");
  dpd_->file2_mat_init(&D);
  for(h=0,count=0; h < nirreps; h++)
    for(a=0; a < D.params->rowtot[h]; a++) 
      for(i=0; i < D.params->coltot[h]; i++) 
	D.matrix[h][a][i] = Z[0][count++];
  dpd_->file2_mat_wrt(&D);
  dpd_->file2_mat_close(&D);
  dpd_->file2_close(&D);

  free_block(Z);
}

void uhf_sf_Zvector(void)
{

}

}} /* End namespaces */
