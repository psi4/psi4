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
#include <cstdio>
#include "psi4/libdpd/dpd.h"
#include <cmath>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libiwl/iwl.h"
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

/* SORTI_ROHF(): Place all the components of the ROHF Lagrangian into
** a large matrix, I (moinfo.I), which we also symmetrize by computing
** Ipq = 1/2 (Ipq + Iqp).  This matrix is later written to disk in
** dump() for subsequent backtransformation.  Note that some of the
** components of the Lagrangian computed into the IIJ, Iij, IIA, and
** Iia matrices remain non-symmetric (e.g., IIJ neq IJI).  I re-used
** my sortone.c code here, so don't let some of the variable names
** confuse you. */

void sortI_ROHF(void)
{
  int h, nirreps, nmo, nfzv, nfzc, nclsd, nopen;
  int row, col, i, j, I, J, a, b, A, B, p, q;
  int *occpi, *virtpi, *occ_off, *vir_off;
  int *occ_sym, *vir_sym, *openpi;
  int *qt_occ, *qt_vir;
  double **O, chksum, value;
  dpdfile2 D;

  nmo = moinfo.nmo;
  nfzc = moinfo.nfzc;
  nfzv = moinfo.nfzv;
  nclsd = moinfo.nclsd;
  nopen = moinfo.nopen;
  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi; virtpi = moinfo.virtpi;
  occ_off = moinfo.occ_off; vir_off = moinfo.vir_off;
  occ_sym = moinfo.occ_sym; vir_sym = moinfo.vir_sym;
  openpi = moinfo.openpi;
  qt_occ = moinfo.qt_occ; qt_vir = moinfo.qt_vir;

  O = block_matrix(nmo,nmo);

  /* Sort alpha components first */
  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, "I(I,J)");
  global_dpd_->file2_mat_init(&D);
  global_dpd_->file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
      for(i=0; i < occpi[h]; i++) {
          I = qt_occ[occ_off[h] + i];
          for(j=0; j < occpi[h]; j++) {
              J = qt_occ[occ_off[h] + j];
              O[I][J] += D.matrix[h][i][j];
            }
        }
    }
  global_dpd_->file2_mat_close(&D);
  global_dpd_->file2_close(&D);

  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, "I'AB");
  global_dpd_->file2_mat_init(&D);
  global_dpd_->file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
      for(a=0; a < (virtpi[h] - openpi[h]); a++) {
          A = qt_vir[vir_off[h] + a];
          for(b=0; b < (virtpi[h] - openpi[h]); b++) {
              B = qt_vir[vir_off[h] + b];

              O[A][B] += D.matrix[h][a][b];
            }
        }
    }
  global_dpd_->file2_mat_close(&D);
  global_dpd_->file2_close(&D);

  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, "I(I,A)");
  global_dpd_->file2_mat_init(&D);
  global_dpd_->file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
      for(i=0; i < occpi[h]; i++) {
          I = qt_occ[occ_off[h] + i];
          for(a=0; a < (virtpi[h] - openpi[h]); a++) {
              A = qt_vir[vir_off[h] + a];

              O[A][I] += D.matrix[h][i][a];
	      O[I][A] += D.matrix[h][i][a];
            }
        }
    }
  global_dpd_->file2_mat_close(&D);
  global_dpd_->file2_close(&D);

  /* Sort beta components */
  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, "I(i,j)");
  global_dpd_->file2_mat_init(&D);
  global_dpd_->file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
      for(i=0; i < (occpi[h] - openpi[h]); i++) {
          I = qt_occ[occ_off[h] + i];
          for(j=0; j < (occpi[h] - openpi[h]); j++) {
              J = qt_occ[occ_off[h] + j];
              O[I][J] += D.matrix[h][i][j];
            }
        }
    }
  global_dpd_->file2_mat_close(&D);
  global_dpd_->file2_close(&D);

  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, "I'ab");
  global_dpd_->file2_mat_init(&D);
  global_dpd_->file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
      for(a=0; a < virtpi[h]; a++) {
          A = qt_vir[vir_off[h] + a];
          for(b=0; b < virtpi[h]; b++) {
              B = qt_vir[vir_off[h] + b];

              O[A][B] += D.matrix[h][a][b];
            }
        }
    }
  global_dpd_->file2_mat_close(&D);
  global_dpd_->file2_close(&D);

  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, "I(i,a)");
  global_dpd_->file2_mat_init(&D);
  global_dpd_->file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
      for(i=0; i < (occpi[h] - openpi[h]); i++) {
          I = qt_occ[occ_off[h] + i];
          for(a=0; a < virtpi[h]; a++) {
              A = qt_vir[vir_off[h] + a];

              O[A][I] += D.matrix[h][i][a];
	      O[I][A] += D.matrix[h][i][a];
            }
        }
    }
  global_dpd_->file2_mat_close(&D);
  global_dpd_->file2_close(&D);

  /* Symmetrize the Lagrangian */
  for(p=0; p < (nmo-nfzv); p++) {
      for(q=0; q < p; q++) {
          value = 0.5*(O[p][q] + O[q][p]);
          O[p][q] = O[q][p] = value;
        }
    }

  /* Multiply the Lagrangian by -2.0 for the final energy derivative
     expression */
  for(p=0; p < (nmo-nfzv); p++) {
      for(q=0; q < (nmo-nfzv); q++) {
	  O[p][q] *= -2.0;
	}
    }

  moinfo.I = O;
}

}} // namespace psi::ccdensity
