/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2018 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
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
#include "psi4/libpsio/psio.h"
#include "psi4/libiwl/iwl.h"
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {
#include "psi4/psifiles.h"

/*
** sortone_rohf(): Place all the components of the 1pdm into a large
** matrix, O (moinfo.opdm), which we also symmetrize by computing Opq
** = 1/2 (Opq + Oqp).  This matrix is later written to disk in dump()
** for subsequent backtransformation.  Note that the components of the
** 1pdm computed into the DIJ, Dij, DAB, Dab, DAI, Dai, DIA, and Dia
** matrices remain non-symmetric (e.g., DIJ neq DJI).
**
** This version doesn't work with frozen orbitals yet.
**
** TDC, 1/03
*/

void sortone_ROHF(struct RHO_Params rho_params)
{
  int h, nirreps, nmo, nfzv, nfzc, nclsd, nopen;
  int row, col, i, j, I, J, a, b, A, B, p, q;
  int *occpi, *virtpi, *occ_off, *vir_off;
  int *occ_sym, *vir_sym, *openpi;
  int *qt_occ, *qt_vir;
  double chksum, value;
  dpdfile2 D;
  psio_address next;

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

  /* O = block_matrix(nmo-nfzc,nmo-nfzc); */
  double **O_a = block_matrix(nmo, nmo);
  double **O_b = block_matrix(nmo, nmo);

  /* Sort A components first */
  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, rho_params.DIJ_lbl);
  global_dpd_->file2_mat_init(&D);
  global_dpd_->file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
      for(i=0; i < occpi[h]; i++) {
          I = qt_occ[occ_off[h] + i];
          for(j=0; j < occpi[h]; j++) {
              J = qt_occ[occ_off[h] + j];
              O_a[I][J] += D.matrix[h][i][j];
            }
        }
    }
  global_dpd_->file2_mat_close(&D);
  global_dpd_->file2_close(&D);

  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, rho_params.DAB_lbl);
  global_dpd_->file2_mat_init(&D);
  global_dpd_->file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
      for(a=0; a < (virtpi[h] - openpi[h]); a++) {
          A = qt_vir[vir_off[h] + a];
          for(b=0; b < (virtpi[h] - openpi[h]); b++) {
              B = qt_vir[vir_off[h] + b];

              O_a[A][B] += D.matrix[h][a][b];
            }
        }
    }
  global_dpd_->file2_mat_close(&D);
  global_dpd_->file2_close(&D);

  /* Note that this component of the density is stored occ-vir */
  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
  global_dpd_->file2_mat_init(&D);
  global_dpd_->file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
      for(i=0; i < occpi[h]; i++) {
          I = qt_occ[occ_off[h] + i];
          for(a=0; a < (virtpi[h] - openpi[h]); a++) {
              A = qt_vir[vir_off[h] + a];

              O_a[A][I] += D.matrix[h][i][a];
            }
        }
    }
  global_dpd_->file2_mat_close(&D);
  global_dpd_->file2_close(&D);

  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
  global_dpd_->file2_mat_init(&D);
  global_dpd_->file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
      for(i=0; i < occpi[h]; i++) {
          I = qt_occ[occ_off[h] + i];
          for(a=0; a < (virtpi[h] - openpi[h]); a++) {
              A = qt_vir[vir_off[h] + a];

              O_a[I][A] += D.matrix[h][i][a];
            }
        }
    }
  global_dpd_->file2_mat_close(&D);
  global_dpd_->file2_close(&D);

  /* Sort B components */
  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, rho_params.Dij_lbl);
  global_dpd_->file2_mat_init(&D);
  global_dpd_->file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
      for(i=0; i < (occpi[h] - openpi[h]); i++) {
          I = qt_occ[occ_off[h] + i];
          for(j=0; j < (occpi[h] - openpi[h]); j++) {
              J = qt_occ[occ_off[h] + j];
              O_b[I][J] += D.matrix[h][i][j];
            }
        }
    }
  global_dpd_->file2_mat_close(&D);
  global_dpd_->file2_close(&D);

  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, rho_params.Dab_lbl);
  global_dpd_->file2_mat_init(&D);
  global_dpd_->file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
      for(a=0; a < virtpi[h]; a++) {
          A = qt_vir[vir_off[h] + a];
          for(b=0; b < virtpi[h]; b++) {
              B = qt_vir[vir_off[h] + b];

              O_b[A][B] += D.matrix[h][a][b];
            }
        }
    }
  global_dpd_->file2_mat_close(&D);
  global_dpd_->file2_close(&D);

  /* Note that this component of the density is stored occ-vir */
  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.Dai_lbl);
  global_dpd_->file2_mat_init(&D);
  global_dpd_->file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
      for(i=0; i < (occpi[h] - openpi[h]); i++) {
          I = qt_occ[occ_off[h] + i];
          for(a=0; a < virtpi[h]; a++) {
              A = qt_vir[vir_off[h] + a];

              O_b[A][I] += D.matrix[h][i][a];
            }
        }
    }
  global_dpd_->file2_mat_close(&D);
  global_dpd_->file2_close(&D);

  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.Dia_lbl);
  global_dpd_->file2_mat_init(&D);
  global_dpd_->file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
      for(i=0; i < (occpi[h] - openpi[h]); i++) {
          I = qt_occ[occ_off[h] + i];
          for(a=0; a < virtpi[h]; a++) {
              A = qt_vir[vir_off[h] + a];

              O_b[I][A] += D.matrix[h][i][a];
            }
        }
    }
  global_dpd_->file2_mat_close(&D);
  global_dpd_->file2_close(&D);

  /* Symmetrize the onepdm */

  for(p=0; p < nmo; p++) {
      for(q=0; q < p; q++) {
          value = 0.5 * (O_a[p][q] + O_a[q][p]);
          O_a[p][q] = O_a[q][p] = value;

          value = 0.5 * (O_b[p][q] + O_b[q][p]);
          O_b[p][q] = O_b[q][p] = value;
        }
    }
  moinfo.opdm_a = O_a;
  moinfo.opdm_b = O_b;

  /*
  for (i=0;i<(nmo-nfzv);++i) {
    for (j=0;j<(nmo-nfzv); ++j) {
      outfile->Printf("%15.10lf\n", O[i][j]);
    }
    outfile->Printf("\n");
  }
  outfile->Printf("End of sortone\n");
	*/
  /*
  for (i=0;i<2;++i) {
    for (j=0;j<5;++j)
      outfile->Printf("%15.10lf",O[i][j]);
    outfile->Printf("\n");
  }
  print_mat(O,nmo-nfzv,nmo-nfzv,outfile);
  */

  double **O = block_matrix(nmo, nmo);
  for(p=0; p < nmo; p++)
    for(q=0; q < nmo; q++)
      O[p][q] = O_a[p][q] + O_b[p][q];
  moinfo.opdm = O;
}

}} // namespace psi::ccdensity
