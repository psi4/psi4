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
    \ingroup CCEOM
    \brief Enter brief description of file here
*/
/*
 root_following: returns index of vector in EOM_Cxxx with maximum overlap with
 vector "CCSD Cxxx" in CC3_MISC
 */

#include <cstdio>
#include <cmath>
#include <cstring>
#include "psi4/libciomr/libciomr.h"
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cceom {

int overlap(int C_irr, int current) {
  dpdfile2 R1, R1A, R1B;
  dpdfile2 R1_old, R1A_old, R1B_old;
  dpdbuf4 R2, R2AA, R2BB, R2AB;
  dpdbuf4 R2_old, R2AA_old, R2BB_old, R2AB_old;
  char lbl[32];

  outfile->Printf( "Overlap of EOM state %d with saved wfns:\n");

  // Current wfn
  if(params.eom_ref == 0) {
    sprintf(lbl, "%s %d %d", "RIA", C_irr, current);
    global_dpd_->file2_init(&R1, PSIF_CC_RAMPS, C_irr, 0, 1, lbl);
    sprintf(lbl, "%s %d %d", "RIjAb", C_irr, current);
    global_dpd_->buf4_init(&R2, PSIF_CC_RAMPS, C_irr, 0, 5, 0, 5, 0, lbl);
  }
  else if(params.eom_ref == 1) {
    sprintf(lbl, "%s %d %d", "RIA", C_irr, current);
    global_dpd_->file2_init(&R1A, PSIF_CC_RAMPS, C_irr, 0, 1, lbl);
    sprintf(lbl, "%s %d %d", "Ria", C_irr, current);
    global_dpd_->file2_init(&R1B, PSIF_CC_RAMPS, C_irr, 0, 1, lbl);

    sprintf(lbl, "%s %d %d", "RIJAB", C_irr, current);
    global_dpd_->buf4_init(&R2AA, PSIF_CC_RAMPS, C_irr, 2, 7, 2, 7, 0, lbl);
    sprintf(lbl, "%s %d %d", "Rijab", C_irr, current);
    global_dpd_->buf4_init(&R2BB, PSIF_CC_RAMPS, C_irr, 2, 7, 2, 7, 0, lbl);
    sprintf(lbl, "%s %d %d", "RIjAb", C_irr, current);
    global_dpd_->buf4_init(&R2AB, PSIF_CC_RAMPS, C_irr, 0, 5, 0, 5, 0, lbl);
  }
  else if(params.eom_ref == 2) {
    sprintf(lbl, "%s %d %d", "RIA", C_irr, current);
    global_dpd_->file2_init(&R1A, PSIF_CC_RAMPS, C_irr, 0, 1, lbl);
    sprintf(lbl, "%s %d %d", "Ria", C_irr, current);
    global_dpd_->file2_init(&R1B, PSIF_CC_RAMPS, C_irr, 2, 3, lbl);

    sprintf(lbl, "%s %d %d", "RIJAB", C_irr, current);
    global_dpd_->buf4_init(&R2AA, PSIF_CC_RAMPS, C_irr, 2, 7, 2, 7, 0, lbl);
    sprintf(lbl, "%s %d %d", "Rijab", C_irr, current);
    global_dpd_->buf4_init(&R2BB, PSIF_CC_RAMPS, C_irr, 12, 17, 12, 17, 0, lbl);
    sprintf(lbl, "%s %d %d", "RIjAb", C_irr, current);
    global_dpd_->buf4_init(&R2AB, PSIF_CC_RAMPS, C_irr, 22, 28, 22, 28, 0, lbl);
  }

  // Stored wfns
  for(int i=0; i < eom_params.cs_per_irrep[C_irr]; i++) {

    if(params.eom_ref == 0) {
      sprintf(lbl, "%s %d %d", "RIA_old", C_irr, i);
      global_dpd_->file2_init(&R1_old, PSIF_CC_RAMPS, C_irr, 0, 1, lbl);
      sprintf(lbl, "%s %d %d", "RIjAb_old", C_irr, i);
      global_dpd_->buf4_init(&R2_old, PSIF_CC_RAMPS, C_irr, 0, 5, 0, 5, 0, lbl);
    }
    else if(params.eom_ref == 1) {
      sprintf(lbl, "%s %d %d", "RIA_old", C_irr, i);
      global_dpd_->file2_init(&R1A_old, PSIF_CC_RAMPS, C_irr, 0, 1, lbl);
      sprintf(lbl, "%s %d %d", "Ria_old", C_irr, i);
      global_dpd_->file2_init(&R1B_old, PSIF_CC_RAMPS, C_irr, 0, 1, lbl);

      sprintf(lbl, "%s %d %d", "RIJAB_old", C_irr, i);
      global_dpd_->buf4_init(&R2AA_old, PSIF_CC_RAMPS, C_irr, 2, 7, 2, 7, 0, lbl);
      sprintf(lbl, "%s %d %d", "Rijab_old", C_irr, i);
      global_dpd_->buf4_init(&R2BB_old, PSIF_CC_RAMPS, C_irr, 2, 7, 2, 7, 0, lbl);
      sprintf(lbl, "%s %d %d", "RIjAb_old", C_irr, i);
      global_dpd_->buf4_init(&R2AB_old, PSIF_CC_RAMPS, C_irr, 0, 5, 0, 5, 0, lbl);
    }
    else if(params.eom_ref == 2) {
      sprintf(lbl, "%s %d %d", "RIA_old", C_irr, i);
      global_dpd_->file2_init(&R1A_old, PSIF_CC_RAMPS, C_irr, 0, 1, lbl);
      sprintf(lbl, "%s %d %d", "Ria_old", C_irr, i);
      global_dpd_->file2_init(&R1B_old, PSIF_CC_RAMPS, C_irr, 2, 3, lbl);

      sprintf(lbl, "%s %d %d", "RIJAB_old", C_irr, i);
      global_dpd_->buf4_init(&R2AA_old, PSIF_CC_RAMPS, C_irr, 2, 7, 2, 7, 0, lbl);
      sprintf(lbl, "%s %d %d", "Rijab_old", C_irr, i);
      global_dpd_->buf4_init(&R2BB_old, PSIF_CC_RAMPS, C_irr, 12, 17, 12, 17, 0, lbl);
      sprintf(lbl, "%s %d %d", "RIjAb_old", C_irr, i);
      global_dpd_->buf4_init(&R2AB_old, PSIF_CC_RAMPS, C_irr, 22, 28, 22, 28, 0, lbl);
    }

    double overlap;
    if (params.eom_ref == 0) {
      overlap = 2.0 * global_dpd_->file2_dot(&R1, &R1_old);
      overlap += global_dpd_->buf4_dot(&R2, &R2_old);
    }
    else if (params.eom_ref == 1 || params.eom_ref == 2) {
      overlap = global_dpd_->file2_dot(&R1A, &R1A_old);
      overlap += global_dpd_->file2_dot(&R1B, &R1B_old);
      overlap += global_dpd_->buf4_dot(&R2AA, &R2AA_old);
      overlap += global_dpd_->buf4_dot(&R2BB, &R2BB_old);
      overlap += global_dpd_->buf4_dot(&R2AB, &R2AB_old);
    }

    outfile->Printf( "State %d --> %5.3f\n", i, fabs(overlap));

    if(params.eom_ref == 0) {
      global_dpd_->file2_close(&R1_old);
      global_dpd_->buf4_close(&R2_old);
    }
    else if(params.eom_ref == 1 || params.eom_ref == 2) {
      global_dpd_->file2_close(&R1A_old);
      global_dpd_->file2_close(&R1B_old);
      global_dpd_->buf4_close(&R2AA_old);
      global_dpd_->buf4_close(&R2BB_old);
      global_dpd_->buf4_close(&R2AB_old);
    }
  }

  if(params.eom_ref == 0) {
    global_dpd_->file2_close(&R1);
    global_dpd_->buf4_close(&R2);
  }
  else if(params.eom_ref == 1 || params.eom_ref == 2) {
    global_dpd_->file2_close(&R1A);
    global_dpd_->file2_close(&R1B);
    global_dpd_->buf4_close(&R2AA);
    global_dpd_->buf4_close(&R2BB);
    global_dpd_->buf4_close(&R2AB);
  }

}

void overlap_stash(int C_irr)
{
  dpdfile2 R1, R1A, R1B;
  dpdbuf4 R2, R2AA, R2BB, R2AB;
  char lbl[32];

  for(int i=0; i < eom_params.cs_per_irrep[C_irr]; i++) {

    if(params.eom_ref == 0) {
      sprintf(lbl, "%s %d %d", "RIA", C_irr, i);
      global_dpd_->file2_init(&R1, PSIF_CC_RAMPS, C_irr, 0, 1, lbl);
      sprintf(lbl, "%s %d %d", "RIA_old", C_irr, i);
      global_dpd_->file2_copy(&R1, PSIF_CC_RAMPS, lbl);

      sprintf(lbl, "%s %d %d", "RIjAb", C_irr, i);
      global_dpd_->buf4_init(&R2, PSIF_CC_RAMPS, C_irr, 0, 5, 0, 5, 0, lbl);
      sprintf(lbl, "%s %d %d", "RIjAb_old", C_irr, i);
      global_dpd_->buf4_copy(&R2, PSIF_CC_RAMPS, lbl);
    }
    else if(params.eom_ref == 1) {
      sprintf(lbl, "%s %d %d", "RIA", C_irr, i);
      global_dpd_->file2_init(&R1A, PSIF_CC_RAMPS, C_irr, 0, 1, lbl);
      sprintf(lbl, "%s %d %d", "RIA_old", C_irr, i);
      global_dpd_->file2_copy(&R1A, PSIF_CC_RAMPS, lbl);

      sprintf(lbl, "%s %d %d", "Ria", C_irr, i);
      global_dpd_->file2_init(&R1B, PSIF_CC_RAMPS, C_irr, 0, 1, lbl);
      sprintf(lbl, "%s %d %d", "Ria_old", C_irr, i);
      global_dpd_->file2_copy(&R1B, PSIF_CC_RAMPS, lbl);

      sprintf(lbl, "%s %d %d", "RIJAB", C_irr, i);
      global_dpd_->buf4_init(&R2AA, PSIF_CC_RAMPS, C_irr, 2, 7, 2, 7, 0, lbl);
      sprintf(lbl, "%s %d %d", "RIJAB_old", C_irr, i);
      global_dpd_->buf4_copy(&R2AA, PSIF_CC_RAMPS, lbl);

      sprintf(lbl, "%s %d %d", "Rijab", C_irr, i);
      global_dpd_->buf4_init(&R2BB, PSIF_CC_RAMPS, C_irr, 2, 7, 2, 7, 0, lbl);
      sprintf(lbl, "%s %d %d", "Rijab_old", C_irr, i);
      global_dpd_->buf4_copy(&R2BB, PSIF_CC_RAMPS, lbl);

      sprintf(lbl, "%s %d %d", "RIjAb", C_irr, i);
      global_dpd_->buf4_init(&R2AB, PSIF_CC_RAMPS, C_irr, 0, 5, 0, 5, 0, lbl);
      sprintf(lbl, "%s %d %d", "RIjAb_old", C_irr, i);
      global_dpd_->buf4_copy(&R2AB, PSIF_CC_RAMPS, lbl);
    }
    else if(params.eom_ref == 2) {
      sprintf(lbl, "%s %d %d", "RIA", C_irr, i);
      global_dpd_->file2_init(&R1A, PSIF_CC_RAMPS, C_irr, 0, 1, lbl);
      sprintf(lbl, "%s %d %d", "RIA_old", C_irr, i);
      global_dpd_->file2_copy(&R1A, PSIF_CC_RAMPS, lbl);

      sprintf(lbl, "%s %d %d", "Ria", C_irr, i);
      global_dpd_->file2_init(&R1B, PSIF_CC_RAMPS, C_irr, 2, 3, lbl);
      sprintf(lbl, "%s %d %d", "Ria_old", C_irr, i);
      global_dpd_->file2_copy(&R1B, PSIF_CC_RAMPS, lbl);

      sprintf(lbl, "%s %d %d", "RIJAB", C_irr, i);
      global_dpd_->buf4_init(&R2AA, PSIF_CC_RAMPS, C_irr, 2, 7, 2, 7, 0, lbl);
      sprintf(lbl, "%s %d %d", "RIJAB_old", C_irr, i);
      global_dpd_->buf4_copy(&R2AA, PSIF_CC_RAMPS, lbl);

      sprintf(lbl, "%s %d %d", "Rijab", C_irr, i);
      global_dpd_->buf4_init(&R2BB, PSIF_CC_RAMPS, C_irr, 12, 17, 12, 17, 0, lbl);
      sprintf(lbl, "%s %d %d", "Rijab_old", C_irr, i);
      global_dpd_->buf4_copy(&R2BB, PSIF_CC_RAMPS, lbl);

      sprintf(lbl, "%s %d %d", "RIjAb", C_irr, i);
      global_dpd_->buf4_init(&R2AB, PSIF_CC_RAMPS, C_irr, 22, 28, 22, 28, 0, lbl);
      sprintf(lbl, "%s %d %d", "RIjAb_old", C_irr, i);
      global_dpd_->buf4_copy(&R2AB, PSIF_CC_RAMPS, lbl);
    }
  }
}

}} // namespace psi::cceom
