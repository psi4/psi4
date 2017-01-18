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

/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here
*/
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

/* relax_I_ROHF(): Add the ROHF orbital-response contributions from
** the one-electron density matrix to the I(I,J) and I(I,A) blocks of
** the Lagrangian.  These terms arise from the first-order CPHF
** equations.  I *think* the following code is general enough to deal
** with both RHF and ROHF cases. */

void relax_I_ROHF(void)
{
  dpdfile2 I, D, f;
  dpdbuf4 E;
  int h, nirreps, i, j, e, *occpi, *virtpi, *openpi;

  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi;
  virtpi = moinfo.virtpi;
  openpi = moinfo.openpi;

  /*** occupied-virtual relaxation terms */

  /* I(I,A) = I'(I,A) - sum_M f(I,M) D(orb)(A,M) */
  global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 1, "I'IA");
  global_dpd_->file2_copy(&I, PSIF_CC_OEI, "I(I,A)");
  global_dpd_->file2_close(&I);
  global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 1, "I(I,A)");
  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 0, "D(orb)(A,I)");
  global_dpd_->file2_init(&f, PSIF_CC_OEI, 0, 0, 0, "fIJ");
  global_dpd_->contract222(&f, &D, &I, 0, 0, -1.0, 1.0);
  global_dpd_->file2_close(&f);
  global_dpd_->file2_close(&D);
  global_dpd_->file2_close(&I);

  /* I(i,a) = I'(i,a) - sum_m f(i,m) D(orb)(a,m) */
  global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 1, "I'ia");
  global_dpd_->file2_copy(&I, PSIF_CC_OEI, "I(i,a)");
  global_dpd_->file2_close(&I);
  global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 1, "I(i,a)");
  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 0, "D(orb)(a,i)");
  global_dpd_->file2_init(&f, PSIF_CC_OEI, 0, 0, 0, "fij");
  global_dpd_->contract222(&f, &D, &I, 0, 0, -1.0, 1.0);
  global_dpd_->file2_close(&f);
  global_dpd_->file2_close(&D);
  global_dpd_->file2_close(&I);

  /*** occupied-occupied relaxtion terms */

  /* I(I,J) <-- I'(I,J) - sum_E,M D(orb)(E,M) [<EI||MJ> + <EJ||MI>]
                      - 2 sum_e,m D(orb)(e,m) <eI|mJ> */
  global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'IJ");
  global_dpd_->file2_copy(&I, PSIF_CC_OEI, "I(I,J)");
  global_dpd_->file2_close(&I);
  global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I(I,J)");
  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 0, "D(orb)(A,I)");
  global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 1, "E <ai|jk>");
  global_dpd_->dot13(&D, &E, &I, 0, 0, -1.0, 1.0);
  global_dpd_->dot13(&D, &E, &I, 0, 1, -1.0, 1.0);
  global_dpd_->buf4_close(&E);
  global_dpd_->file2_close(&D);
  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 0, "D(orb)(a,i)");
  global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
  global_dpd_->dot13(&D, &E, &I, 0, 0, -2.0, 1.0);
  global_dpd_->buf4_close(&E);
  global_dpd_->file2_close(&D);

  /* I(I,J) <-- - 2 sum_E  f(E,I) D(orb)(E,J) (J,socc)

   Note that this same term is not needed in the I(i,j) block since J
   is required to be a singly occupied orbital */
  global_dpd_->file2_mat_init(&I);
  global_dpd_->file2_mat_rd(&I);

  global_dpd_->file2_init(&f, PSIF_CC_OEI, 0, 0, 1, "fIA");
  global_dpd_->file2_mat_init(&f);
  global_dpd_->file2_mat_rd(&f);
  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 0, "D(orb)(A,I)");
  global_dpd_->file2_mat_init(&D);
  global_dpd_->file2_mat_rd(&D);

  for(h=0; h < nirreps; h++) {
      for(i=0; i < occpi[h]; i++)
	  for(j=(occpi[h] - openpi[h]); j < occpi[h]; j++)
	      for(e=0; e < virtpi[h]; e++)
		  I.matrix[h][i][j] -= 2 * f.matrix[h][i][e] * D.matrix[h][e][j];
    }

  global_dpd_->file2_mat_close(&D);
  global_dpd_->file2_close(&D);
  global_dpd_->file2_mat_close(&f);
  global_dpd_->file2_close(&f);

  global_dpd_->file2_mat_wrt(&I);
  global_dpd_->file2_mat_close(&I);
  global_dpd_->file2_close(&I);

  /* I(i,j) <-- I'(i,j) - sum_e,m D(orb)(e,m) [<ei||mj> + <ej||mi>]
                      - 2 sum_E,M D(orb)(E,M) <Ei|Mj> */
  global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'ij");
  global_dpd_->file2_copy(&I, PSIF_CC_OEI, "I(i,j)");
  global_dpd_->file2_close(&I);
  global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I(i,j)");
  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 0, "D(orb)(a,i)");
  global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 1, "E <ai|jk>");
  global_dpd_->dot13(&D, &E, &I, 0, 0, -1.0, 1.0);
  global_dpd_->dot13(&D, &E, &I, 0, 1, -1.0, 1.0);
  global_dpd_->buf4_close(&E);
  global_dpd_->file2_close(&D);
  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 0, "D(orb)(A,I)");
  global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
  global_dpd_->dot13(&D, &E, &I, 0, 0, -2.0, 1.0);
  global_dpd_->buf4_close(&E);
  global_dpd_->file2_close(&D);
  global_dpd_->file2_close(&I);

  /* Clean the I(i,j) block yet again */
  global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I(i,j)");
  global_dpd_->file2_mat_init(&I);
  global_dpd_->file2_mat_rd(&I);

  for(h=0; h < nirreps; h++)
      for(i=0; i < occpi[h]; i++)
	  for(j=0; j < occpi[h]; j++)
	      if((i >= (occpi[h] - openpi[h])) ||
		 (j >= (occpi[h] - openpi[h])) )
		  I.matrix[h][i][j] = 0.0;

  global_dpd_->file2_mat_wrt(&I);
  global_dpd_->file2_mat_close(&I);
  global_dpd_->file2_close(&I);
}

}} // namespace psi::ccdensity
