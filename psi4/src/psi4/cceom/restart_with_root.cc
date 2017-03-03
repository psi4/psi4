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
 restart_with_root: copies C's from position prop_root to position 0 
 in EOM_Cxxx files 

 also put copy in CC3_MISC for root_following
*/

#include <cstdio>
#include <cmath>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cceom {

void restart_with_root(int prop_root, int C_irr) {
  dpdfile2 CME, Cme;
  dpdbuf4 CMNEF, Cmnef, CMnEf;
  char lbl[32];

  outfile->Printf("Copying root %d to start of EOM_Cxxx files.\n",prop_root+1);

  if (params.eom_ref == 0) {
    sprintf(lbl, "CME %d", prop_root);
    global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, lbl);
    global_dpd_->file2_copy(&CME, PSIF_EOM_CME, "CME 0");
    global_dpd_->file2_close(&CME);

    sprintf(lbl, "CMnEf %d", prop_root);
    global_dpd_->buf4_init(&CMnEf, PSIF_EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_copy(&CMnEf, PSIF_EOM_CMnEf, "CMnEf 0");
    global_dpd_->buf4_close(&CMnEf);
  }
  else if (params.eom_ref == 1) {
    sprintf(lbl, "CME %d", prop_root);
    global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, lbl);
    global_dpd_->file2_copy(&CME, PSIF_EOM_CME, "CME 0");
    global_dpd_->file2_close(&CME);

    sprintf(lbl, "Cme %d", prop_root);
    global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 0, 1, lbl);
    global_dpd_->file2_copy(&Cme, PSIF_EOM_Cme, "Cme 0");
    global_dpd_->file2_close(&Cme);

    sprintf(lbl, "CMNEF %d", prop_root);
    global_dpd_->buf4_init(&CMNEF, PSIF_EOM_CMNEF, C_irr, 2, 7, 2, 7, 0, lbl);
    global_dpd_->buf4_copy(&CMNEF, PSIF_EOM_CMNEF, "CMNEF 0");
    global_dpd_->buf4_close(&CMNEF);

    sprintf(lbl, "Cmnef %d", prop_root);
    global_dpd_->buf4_init(&Cmnef, PSIF_EOM_Cmnef, C_irr, 2, 7, 2, 7, 0, lbl);
    global_dpd_->buf4_copy(&Cmnef, PSIF_EOM_Cmnef, "Cmnef 0");
    global_dpd_->buf4_close(&Cmnef);

    sprintf(lbl, "CMnEf %d", prop_root);
    global_dpd_->buf4_init(&CMnEf, PSIF_EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_copy(&CMnEf, PSIF_EOM_CMnEf, "CMnEf 0");
    global_dpd_->buf4_close(&CMnEf);
  }
  else {
    sprintf(lbl, "CME %d", prop_root);
    global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, lbl);
    global_dpd_->file2_copy(&CME, PSIF_EOM_CME, "CME 0");
    global_dpd_->file2_close(&CME);

    sprintf(lbl, "Cme %d", prop_root);
    global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 2, 3, lbl);
    global_dpd_->file2_copy(&Cme, PSIF_EOM_Cme, "Cme 0");
    global_dpd_->file2_close(&Cme);

    sprintf(lbl, "CMNEF %d", prop_root);
    global_dpd_->buf4_init(&CMNEF, PSIF_EOM_CMNEF, C_irr, 2, 7, 2, 7, 0, lbl);
    global_dpd_->buf4_copy(&CMNEF, PSIF_EOM_CMNEF, "CMNEF 0");
    global_dpd_->buf4_close(&CMNEF);

    sprintf(lbl, "Cmnef %d", prop_root);
    global_dpd_->buf4_init(&Cmnef, PSIF_EOM_Cmnef, C_irr, 12, 17, 12, 17, 0, lbl);
    global_dpd_->buf4_copy(&Cmnef, PSIF_EOM_Cmnef, "Cmnef 0");
    global_dpd_->buf4_close(&Cmnef);

    sprintf(lbl, "CMnEf %d", prop_root);
    global_dpd_->buf4_init(&CMnEf, PSIF_EOM_CMnEf, C_irr, 22, 28, 22, 28, 0, lbl);
    global_dpd_->buf4_copy(&CMnEf, PSIF_EOM_CMnEf, "CMnEf 0");
    global_dpd_->buf4_close(&CMnEf);
  }
  return;
}

/*
 save_C_ccsd: copies C's from position prop_root to CC3_MISC file
*/

void save_C_ccsd(int prop_root, int C_irr) {
  dpdfile2 CME, Cme;
  dpdbuf4 CMNEF, Cmnef, CMnEf;
  char lbl[32];

  outfile->Printf("Copying root %d to CC3_MISC file.\n",prop_root+1);

  if (params.eom_ref == 0) {
    sprintf(lbl, "CME %d", prop_root);
    global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, lbl);
    global_dpd_->file2_copy(&CME, PSIF_CC3_MISC, "CCSD CME");
    global_dpd_->file2_close(&CME);

    sprintf(lbl, "CMnEf %d", prop_root);
    global_dpd_->buf4_init(&CMnEf, PSIF_EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_copy(&CMnEf, PSIF_CC3_MISC, "CCSD CMnEf");
    global_dpd_->buf4_close(&CMnEf);
  }
  else if (params.eom_ref == 1) {
    sprintf(lbl, "CME %d", prop_root);
    global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, lbl);
    global_dpd_->file2_copy(&CME, PSIF_CC3_MISC, "CCSD CME");
    global_dpd_->file2_close(&CME);

    sprintf(lbl, "Cme %d", prop_root);
    global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 0, 1, lbl);
    global_dpd_->file2_copy(&Cme, PSIF_CC3_MISC, "CCSD Cme");
    global_dpd_->file2_close(&Cme);

    sprintf(lbl, "CMNEF %d", prop_root);
    global_dpd_->buf4_init(&CMNEF, PSIF_EOM_CMNEF, C_irr, 2, 7, 2, 7, 0, lbl);
    global_dpd_->buf4_copy(&CMNEF, PSIF_CC3_MISC, "CCSD CMNEF");
    global_dpd_->buf4_close(&CMNEF);

    sprintf(lbl, "Cmnef %d", prop_root);
    global_dpd_->buf4_init(&Cmnef, PSIF_EOM_Cmnef, C_irr, 2, 7, 2, 7, 0, lbl);
    global_dpd_->buf4_copy(&Cmnef, PSIF_CC3_MISC, "CCSD Cmnef");
    global_dpd_->buf4_close(&Cmnef);

    sprintf(lbl, "CMnEf %d", prop_root);
    global_dpd_->buf4_init(&CMnEf, PSIF_EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_copy(&CMnEf, PSIF_CC3_MISC, "CCSD CMnEf");
    global_dpd_->buf4_close(&CMnEf);
  }
  else {
    sprintf(lbl, "CME %d", prop_root);
    global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, lbl);
    global_dpd_->file2_copy(&CME, PSIF_CC3_MISC, "CCSD CME");
    global_dpd_->file2_close(&CME);

    sprintf(lbl, "Cme %d", prop_root);
    global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 2, 3, lbl);
    global_dpd_->file2_copy(&Cme, PSIF_CC3_MISC, "CCSD Cme");
    global_dpd_->file2_close(&Cme);

    sprintf(lbl, "CMNEF %d", prop_root);
    global_dpd_->buf4_init(&CMNEF, PSIF_EOM_CMNEF, C_irr, 2, 7, 2, 7, 0, lbl);
    global_dpd_->buf4_copy(&CMNEF, PSIF_CC3_MISC, "CCSD CMNEF");
    global_dpd_->buf4_close(&CMNEF);

    sprintf(lbl, "Cmnef %d", prop_root);
    global_dpd_->buf4_init(&Cmnef, PSIF_EOM_Cmnef, C_irr, 12, 17, 12, 17, 0, lbl);
    global_dpd_->buf4_copy(&Cmnef, PSIF_CC3_MISC, "CCSD Cmnef");
    global_dpd_->buf4_close(&Cmnef);

    sprintf(lbl, "CMnEf %d", prop_root);
    global_dpd_->buf4_init(&CMnEf, PSIF_EOM_CMnEf, C_irr, 22, 28, 22, 28, 0, lbl);
    global_dpd_->buf4_copy(&CMnEf, PSIF_CC3_MISC, "CCSD CMnEf");
    global_dpd_->buf4_close(&CMnEf);
  }
  return;
}

}} // namespace psi::cceom
