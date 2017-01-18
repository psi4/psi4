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

/* lag(): Build the orbital Lagrangian, I'pq, defined in spin-orbitals
** as:
**
** I'pq = sum_r fpr (Dqr + Drq) + sum_rs <pr||qs> (Drs + Dsr) (q,occ)
**        + sum_rst <pr||st> Gqrst + 2 fpq (q,occ)
**
** The orbital-response component of the gradient (for non-correlated
** orbitals) is defined as
**
** dE/dx <--- sum_pq I'pq U(x)pq
**
** where U(x)pq is the usual CPHF coefficient.  Note, however, that
** the final expression we want involves not CPHF coefficients, but
** overlap derivatives.  For example, in the occupied-occupied and
** virtual-virtual blocks, a choice of non-canonical perturbed
** orbitals allows the assignments
**
** U(x)ij = -1/2 S(x)ij      and      U(x)ab = -1/2 S(x)ab
**
** to be made.  We also choose to incorporate the -1/2 prefactor
** into the Largrangian itself so the final orbital response
** expression will appear as
**
** dE/dx <--- sum_pq Ipq S(x)pq
**
** where Ipq is the "relaxed" Lagrangian (see relax_I.c).
**
** The final set of loops force the appropriate open-shell terms to
** zero for ROHF refs. (See the description of the treatment of
** open-shells in the ROHF-CCSD code as discussed in CCSORT for an
** explanation of why this is necessary.) */

void Iij(struct RHO_Params rho_params);
void Iab(struct RHO_Params rho_params);
void Iai(struct RHO_Params rho_params);
void Iia(struct RHO_Params rho_params);

void lag(struct RHO_Params rho_params)
{
  int h, nirreps, i, j, a, b;
  int *occpi, *virtpi, *openpi;
  dpdfile2 I;

  Iij(rho_params);
  Iab(rho_params);
  Iai(rho_params);
  Iia(rho_params);

  /* Multiply all I'pq components by -1/2 for compatibility with the
     final gradient expression */

  if(params.ref == 0) { /** RHF **/
    global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'IJ");
    global_dpd_->file2_scm(&I, -0.5);
    global_dpd_->file2_close(&I);
    global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 1, "I'AB");
    global_dpd_->file2_scm(&I, -0.5);
    global_dpd_->file2_close(&I);
    global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 1, "I'IA");
    global_dpd_->file2_scm(&I, -0.5);
    global_dpd_->file2_close(&I);
    global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 0, "I'AI");
    global_dpd_->file2_scm(&I, -0.5);
    global_dpd_->file2_close(&I);
  }
  else if(params.ref == 1) { /** ROHF **/

    global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'IJ");
    global_dpd_->file2_scm(&I, -0.5);
    global_dpd_->file2_close(&I);
    global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'ij");
    global_dpd_->file2_scm(&I, -0.5);
    global_dpd_->file2_close(&I);
    global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 1, "I'AB");
    global_dpd_->file2_scm(&I, -0.5);
    global_dpd_->file2_close(&I);
    global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 1, "I'ab");
    global_dpd_->file2_scm(&I, -0.5);
    global_dpd_->file2_close(&I);
    global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 1, "I'IA");
    global_dpd_->file2_scm(&I, -0.5);
    global_dpd_->file2_close(&I);
    global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 1, "I'ia");
    global_dpd_->file2_scm(&I, -0.5);
    global_dpd_->file2_close(&I);
    global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 0, "I'AI");
    global_dpd_->file2_scm(&I, -0.5);
    global_dpd_->file2_close(&I);
    global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 0, "I'ai");
    global_dpd_->file2_scm(&I, -0.5);
    global_dpd_->file2_close(&I);
  }
  else if(params.ref == 2) { /** UHF **/

    global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'IJ");
    global_dpd_->file2_scm(&I, -0.5);
    global_dpd_->file2_close(&I);
    global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 2, 2, "I'ij");
    global_dpd_->file2_scm(&I, -0.5);
    global_dpd_->file2_close(&I);
    global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 1, "I'AB");
    global_dpd_->file2_scm(&I, -0.5);
    global_dpd_->file2_close(&I);
    global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 3, 3, "I'ab");
    global_dpd_->file2_scm(&I, -0.5);
    global_dpd_->file2_close(&I);
    global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 1, "I'IA");
    global_dpd_->file2_scm(&I, -0.5);
    global_dpd_->file2_close(&I);
    global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 2, 3, "I'ia");
    global_dpd_->file2_scm(&I, -0.5);
    global_dpd_->file2_close(&I);
    global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 0, "I'AI");
    global_dpd_->file2_scm(&I, -0.5);
    global_dpd_->file2_close(&I);
    global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 3, 2, "I'ai");
    global_dpd_->file2_scm(&I, -0.5);
    global_dpd_->file2_close(&I);

  }

  /* Now go through all terms involving open-shell orbitals and force
     the appropriate spin cases to zero. */

  if(params.ref == 1) { /** ROHF **/
    nirreps = moinfo.nirreps;
    occpi = moinfo.occpi; virtpi = moinfo.virtpi; openpi = moinfo.openpi;


    global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'ij");
    global_dpd_->file2_mat_init(&I);
    global_dpd_->file2_mat_rd(&I);
    for(h=0; h < nirreps; h++) {
      for(i=(occpi[h]-openpi[h]); i < occpi[h]; i++) {
	for(j=(occpi[h]-openpi[h]); j < occpi[h]; j++) {
	  I.matrix[h][i][j] = 0.0;
	}
      }
      for(i=(occpi[h]-openpi[h]); i < occpi[h]; i++) {
	for(j=0; j < occpi[h]; j++) {
	  I.matrix[h][i][j] = 0.0;
	}
      }
      for(i=0; i < occpi[h]; i++) {
	for(j=(occpi[h]-openpi[h]); j < occpi[h]; j++) {
	  I.matrix[h][i][j] = 0.0;
	}
      }
    }
    global_dpd_->file2_mat_wrt(&I);
    global_dpd_->file2_mat_close(&I);
    global_dpd_->file2_close(&I);

    global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 1, "I'AB");
    global_dpd_->file2_mat_init(&I);
    global_dpd_->file2_mat_rd(&I);
    for(h=0; h < nirreps; h++) {
      for(a=(virtpi[h]-openpi[h]); a < virtpi[h]; a++) {
	for(b=(virtpi[h]-openpi[h]); b < virtpi[h]; b++) {
	  I.matrix[h][a][b] = 0.0;
	}
      }
      for(a=(virtpi[h]-openpi[h]); a < virtpi[h]; a++) {
	for(b=0; b < virtpi[h]; b++) {
	  I.matrix[h][a][b] = 0.0;
	}
      }
      for(a=0; a < virtpi[h]; a++) {
	for(b=(virtpi[h]-openpi[h]); b < virtpi[h]; b++) {
	  I.matrix[h][a][b] = 0.0;
	}
      }
    }
    global_dpd_->file2_mat_wrt(&I);
    global_dpd_->file2_mat_close(&I);
    global_dpd_->file2_close(&I);

    global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 1, "I'ab");
    global_dpd_->file2_close(&I);

    global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 0, "I'AI");
    global_dpd_->file2_mat_init(&I);
    global_dpd_->file2_mat_rd(&I);
    for(h=0; h < nirreps; h++) {
      for(a=(virtpi[h]-openpi[h]); a < virtpi[h]; a++) {
	for(i=0; i < occpi[h]; i++) {
	  I.matrix[h][a][i] = 0.0;
	}
      }
    }
    global_dpd_->file2_mat_wrt(&I);
    global_dpd_->file2_mat_close(&I);
    global_dpd_->file2_close(&I);

    global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 0, "I'ai");
    global_dpd_->file2_mat_init(&I);
    global_dpd_->file2_mat_rd(&I);
    for(h=0; h < nirreps; h++) {
      for(a=0; a < virtpi[h]; a++) {
	for(i=(occpi[h] - openpi[h]); i < occpi[h]; i++) {
	  I.matrix[h][a][i] = 0.0;
	}
      }
    }
    global_dpd_->file2_mat_wrt(&I);
    global_dpd_->file2_mat_close(&I);
    global_dpd_->file2_close(&I);

    global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 1, "I'IA");
    global_dpd_->file2_mat_init(&I);
    global_dpd_->file2_mat_rd(&I);
    for(h=0; h < nirreps; h++) {
      for(i=0; i < occpi[h]; i++) {
	for(a=(virtpi[h] - openpi[h]); a < virtpi[h]; a++) {
	  I.matrix[h][i][a] = 0.0;
	}
      }
    }
    global_dpd_->file2_mat_wrt(&I);
    global_dpd_->file2_mat_close(&I);
    global_dpd_->file2_close(&I);

    global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 1, "I'ia");
    global_dpd_->file2_mat_init(&I);
    global_dpd_->file2_mat_rd(&I);
    for(h=0; h < nirreps; h++) {
      for(i=(occpi[h] - openpi[h]); i < occpi[h]; i++) {
	for(a=0; a < virtpi[h]; a++) {
	  I.matrix[h][i][a] = 0.0;
	}
      }
    }
    global_dpd_->file2_mat_wrt(&I);
    global_dpd_->file2_mat_close(&I);
    global_dpd_->file2_close(&I);

  }
}

}} // namespace psi::ccdensity
