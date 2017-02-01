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
    \ingroup CCENERGY
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "MOInfo.h"
#include "Params.h"
#include "ccwave.h"

namespace psi { namespace ccenergy {

/* Computes a modified D1 diagnostic developed by T.J. Lee, but not yet
 * published.
 * */

double CCEnergyWavefunction::new_d1diag_t1_rohf(void)
{
  int h, nirreps, i, j;
  int nclsd, nuocc, nopen;
  double **T1_hp, **T1_hx, **T1_xp, **T1_sq;
  double *E, **C;
  double max_hp=0.0, max_xp=0.0, max_hx=0.0, max;
  dpdfile2 T1_a, T1_b;

  nirreps = moinfo_.nirreps;

  global_dpd_->file2_init(&T1_a, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->file2_mat_init(&T1_a);
  global_dpd_->file2_mat_rd(&T1_a);

  global_dpd_->file2_init(&T1_b, PSIF_CC_OEI, 0, 0, 1, "tia");
  global_dpd_->file2_mat_init(&T1_b);
  global_dpd_->file2_mat_rd(&T1_b);

  for(h=0; h < nirreps; h++) {
    nclsd = moinfo_.clsdpi[h];
    nuocc = moinfo_.uoccpi[h];
    nopen = moinfo_.openpi[h];

    if(nclsd && nuocc) {
      T1_hp = block_matrix(nclsd, nuocc);
      for (i=0; i < nclsd; i++)
	for (j=0; j < nuocc; j++)
	  T1_hp[i][j] = (T1_a.matrix[h][i][j] + T1_b.matrix[h][i][j])/2.;

      T1_sq = block_matrix(nclsd, nclsd);
      C_DGEMM('n','t',nclsd,nclsd,nuocc,1.0,&(T1_hp[0][0]),nuocc,
	     &(T1_hp[0][0]),nuocc,0.0,&(T1_sq[0][0]),nclsd);

      E = init_array(nclsd);
      C = block_matrix(nclsd, nclsd);
      sq_rsp(nclsd, nclsd, T1_sq, E, 0, C, 1e-12);
      for(i=0; i < nclsd; i++) if(E[i] > max_hp) max_hp = E[i];
      free(E);
      free_block(C);
      free_block(T1_sq);
      free_block(T1_hp);
    }

    if(nclsd && nopen) {
      T1_hx = block_matrix(nclsd, nopen);
      for (i=0; i < nclsd; i++)
	for (j=0; j < nopen; j++)
	  T1_hx[i][j] = T1_b.matrix[h][i][nuocc+j]/sqrt(2.);

      T1_sq = block_matrix(nclsd, nclsd);
      C_DGEMM('n','t',nclsd,nclsd,nopen,1.0,&(T1_hx[0][0]),nopen,
	     &(T1_hx[0][0]),nopen,0.0,&(T1_sq[0][0]),nclsd);

      E = init_array(nclsd);
      C = block_matrix(nclsd, nclsd);
      sq_rsp(nclsd, nclsd, T1_sq, E, 0, C, 1e-12);
      for(i=0; i < nclsd; i++) if(E[i] > max_hx) max_hx = E[i];
      free(E);
      free_block(C);
      free_block(T1_sq);
      free_block(T1_hx);
    }

    if(nopen && nuocc) {
      T1_xp = block_matrix(nopen, nuocc);
      for (i=0; i < nopen; i++)
	for (j=0; j < nuocc; j++)
	  T1_xp[i][j] = T1_a.matrix[h][nclsd+i][j]/sqrt(2.);

      T1_sq = block_matrix(nopen, nopen);
      C_DGEMM('n','t',nopen,nopen,nuocc,1.0,&(T1_xp[0][0]),nuocc,
	     &(T1_xp[0][0]),nuocc,0.0,&(T1_sq[0][0]),nopen);

      E = init_array(nopen);
      C = block_matrix(nopen, nopen);
      sq_rsp(nopen, nopen, T1_sq, E, 0, C, 1e-12);
      for(i=0; i < nopen; i++) if(E[i] > max_xp) max_xp = E[i];
      free(E);
      free_block(C);
      free_block(T1_sq);
      free_block(T1_xp);
    }
  }

  global_dpd_->file2_mat_close(&T1_a);
  global_dpd_->file2_close(&T1_a);

  global_dpd_->file2_mat_close(&T1_b);
  global_dpd_->file2_close(&T1_b);

  max_hp = sqrt(max_hp);
  max_hx = sqrt(max_hx);
  max_xp = sqrt(max_xp);

  /*
  outfile->Printf( "ND1: hp=%8.6f hx=%8.6f xp=%8.6f\n", max_hp, max_hx, max_xp);
  */

  max = max_hp;
  if (max_hx > max) max = max_hx;
  if (max_xp > max) max = max_xp;

  return max;
}

double CCEnergyWavefunction::new_d1diag(void)
{
  double norm = 0.0;

  if(params_.ref == 0) { /** RHF **/
    norm = d1diag_t1_rhf();
  }
  else if (params_.ref == 1) { /** ROHF **/
    norm = new_d1diag_t1_rohf();
  }
  return norm;
}
}} // namespace psi::ccenergy
