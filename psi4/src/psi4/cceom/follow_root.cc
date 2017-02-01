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
#include "psi4/libciomr/libciomr.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cceom {

int follow_root(int L, double **alpha, int C_irr) {
  dpdfile2 CME, Cme, RME, Rme;
  dpdbuf4 CMNEF, Cmnef, CMnEf, RMNEF, Rmnef, RMnEf;
  char lbl[32];
  double *CR_overlap, tval;
  int i,j;

  CR_overlap = init_array(L);

  /* open CCSD vector "R" from CC3_MISC */
  if (params.eom_ref == 0) {
    global_dpd_->file2_init(&RME, PSIF_CC3_MISC, C_irr, 0, 1, "CCSD CME");
    global_dpd_->buf4_init(&RMnEf, PSIF_CC3_MISC, C_irr, 0, 5, 0, 5, 0, "CCSD CMnEf");
  }
  else if (params.eom_ref == 1) {
    global_dpd_->file2_init(&RME, PSIF_CC3_MISC, C_irr, 0, 1, "CCSD CME");
    global_dpd_->file2_init(&Rme, PSIF_CC3_MISC, C_irr, 0, 1, "CCSD Cme");
    global_dpd_->buf4_init(&RMNEF, PSIF_CC3_MISC, C_irr, 2, 7, 2, 7, 0, "CCSD CMNEF");
    global_dpd_->buf4_init(&Rmnef, PSIF_CC3_MISC, C_irr, 2, 7, 2, 7, 0, "CCSD Cmnef");
    global_dpd_->buf4_init(&RMnEf, PSIF_CC3_MISC, C_irr, 0, 5, 0, 5, 0, "CCSD CMnEf");
  }
  else if (params.eom_ref == 2) {
    global_dpd_->file2_init(&RME, PSIF_CC3_MISC, C_irr, 0, 1, "CCSD CME");
    global_dpd_->file2_init(&Rme, PSIF_CC3_MISC, C_irr, 2, 3, "CCSD Cme");
    global_dpd_->buf4_init(&RMNEF, PSIF_CC3_MISC, C_irr, 2, 7, 2, 7, 0, "CCSD CMNEF");
    global_dpd_->buf4_init(&Rmnef, PSIF_CC3_MISC, C_irr, 12, 17, 12, 17, 0, "CCSD Cmnef");
    global_dpd_->buf4_init(&RMnEf, PSIF_CC3_MISC, C_irr, 22, 28, 22, 28, 0, "CCSD CMnEf");
  }

  /* loop over trial C vectors */
  for (i=0; i<L; ++i) {

    /* read C vector from EOM_Cxxx */
    if (params.eom_ref == 0) {
      sprintf(lbl, "%s %d", "CME", i);
      global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, lbl);
      sprintf(lbl, "%s %d", "CMnEf", i);
      global_dpd_->buf4_init(&CMnEf, PSIF_EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, lbl);
    }
    else if (params.eom_ref == 1) {
      sprintf(lbl, "%s %d", "CME", i);
      global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, lbl);
      sprintf(lbl, "%s %d", "Cme", i);
      global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 0, 1, lbl);
      sprintf(lbl, "%s %d", "CMNEF", i);
      global_dpd_->buf4_init(&CMNEF, PSIF_EOM_CMNEF, C_irr, 2, 7, 2, 7, 0, lbl);
      sprintf(lbl, "%s %d", "Cmnef", i);
      global_dpd_->buf4_init(&Cmnef, PSIF_EOM_Cmnef, C_irr, 2, 7, 2, 7, 0, lbl);
      sprintf(lbl, "%s %d", "CMnEf", i);
      global_dpd_->buf4_init(&CMnEf, PSIF_EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, lbl);
    }
    else if (params.eom_ref == 2) {
      sprintf(lbl, "%s %d", "CME", i);
      global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, lbl);
      sprintf(lbl, "%s %d", "Cme", i);
      global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 2, 3, lbl);
      sprintf(lbl, "%s %d", "CMNEF", i);
      global_dpd_->buf4_init(&CMNEF, PSIF_EOM_CMNEF, C_irr, 2, 7, 2, 7, 0, lbl);
      sprintf(lbl, "%s %d", "Cmnef", i);
      global_dpd_->buf4_init(&Cmnef, PSIF_EOM_Cmnef, C_irr, 12, 17, 12, 17, 0, lbl);
      sprintf(lbl, "%s %d", "CMnEf", i);
      global_dpd_->buf4_init(&CMnEf, PSIF_EOM_CMnEf, C_irr, 22, 28, 22, 28, 0, lbl);
    }

    /* dot C vector with R vector */
    tval = 0.0;
    if (params.eom_ref == 0) {
      tval = 2.0 * global_dpd_->file2_dot(&CME, &RME);
      tval += global_dpd_->buf4_dot(&CMnEf, &RMnEf);
    }
    else if (params.eom_ref == 1) {
      tval = global_dpd_->file2_dot(&CME, &RME);
      tval += global_dpd_->file2_dot(&Cme, &Rme);
      tval += global_dpd_->buf4_dot(&CMNEF, &RMNEF);
      tval += global_dpd_->buf4_dot(&Cmnef, &Rmnef);
      tval += global_dpd_->buf4_dot(&CMnEf, &RMnEf);
    }
    else if (params.eom_ref == 2) {
      tval = global_dpd_->file2_dot(&CME, &RME);
      tval += global_dpd_->file2_dot(&Cme, &Rme);
      tval += global_dpd_->buf4_dot(&CMNEF, &RMNEF);
      tval += global_dpd_->buf4_dot(&Cmnef, &Rmnef);
      tval += global_dpd_->buf4_dot(&CMnEf, &RMnEf);
    }

    /* loop over roots and add in overlap */
    for (j=0; j<L; ++j)
      CR_overlap[j] += alpha[i][j] * tval;

    if (params.eom_ref == 0) {
      global_dpd_->file2_close(&CME);
      global_dpd_->buf4_close(&CMnEf);
    }
    else {
      global_dpd_->file2_close(&CME);
      global_dpd_->file2_close(&Cme);
      global_dpd_->buf4_close(&CMNEF);
      global_dpd_->buf4_close(&Cmnef);
      global_dpd_->buf4_close(&CMnEf);
    }
  }

  if (params.eom_ref == 0) {
    global_dpd_->file2_close(&RME);
    global_dpd_->buf4_close(&RMnEf);
  }
  else {
    global_dpd_->file2_close(&RME);
    global_dpd_->file2_close(&Rme);
    global_dpd_->buf4_close(&RMNEF);
    global_dpd_->buf4_close(&Rmnef);
    global_dpd_->buf4_close(&RMnEf);
  }

  outfile->Printf("Overlaps of Rs with EOM CCSD eigenvector:\n");
  for(i=0;i<L;++i) {
    outfile->Printf("\t %d  %12.6lf\n", i, CR_overlap[i]);
  }

  /* return index with greatest overlap */
  tval = -1.0;

  for(i=0;i<L;++i) {
    if ( fabs(CR_overlap[i]) > tval) {
      tval = fabs(CR_overlap[i]);
      j = i;
    }
  }

  outfile->Printf("follow_root returning: %d\n", j);
  return j;
}

}} // namespace psi::cceom
