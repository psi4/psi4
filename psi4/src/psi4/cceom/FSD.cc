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
#include <cstdio>
#include <cmath>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cceom {

/* This function computes the H-bar singles-doubles block contribution
   to a Sigma vector stored at Sigma plus 'i' */

void FSD(int i, int C_irr) {
  dpdfile2 Fme, FME;
  dpdfile2 SIA, Sia;
  dpdbuf4 CMNEF, Cmnef, CMnEf, CmNeF;
  char lbl[32];

  if (params.eom_ref == 0) {  /* RHF */
    sprintf(lbl, "%s %d", "SIA", i);
    global_dpd_->file2_init(&SIA, PSIF_EOM_SIA, C_irr, 0, 1, lbl);
    global_dpd_->file2_init(&FME, PSIF_CC_OEI, H_IRR, 0, 1, "FME");
    global_dpd_->buf4_init(&CMnEf, PSIF_EOM_TMP, C_irr, 0, 5, 0, 5, 0, "2CMnEf - CMnfE");
    global_dpd_->dot24(&FME,&CMnEf,&SIA, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&CMnEf);
    global_dpd_->file2_close(&FME);
    global_dpd_->file2_close(&SIA);
  }

  else if (params.eom_ref == 1) { /* ROHF */
    sprintf(lbl, "%s %d", "SIA", i);
    global_dpd_->file2_init(&SIA, PSIF_EOM_SIA, C_irr, 0, 1, lbl);
    sprintf(lbl, "%s %d", "Sia", i);
    global_dpd_->file2_init(&Sia, PSIF_EOM_Sia, C_irr, 0, 1, lbl);

    /* SIA += FME*CIMAE + Fme*CImAe */
    global_dpd_->file2_init(&FME, PSIF_CC_OEI, H_IRR, 0, 1, "FME");
    sprintf(lbl, "%s %d", "CMNEF", i);
    global_dpd_->buf4_init(&CMNEF, PSIF_EOM_CMNEF, C_irr, 0, 5, 2, 7, 0, lbl);
    global_dpd_->dot24(&FME,&CMNEF,&SIA, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&CMNEF);
    global_dpd_->file2_close(&FME);

    global_dpd_->file2_init(&Fme, PSIF_CC_OEI, H_IRR, 0, 1, "Fme");
    sprintf(lbl, "%s %d", "CMnEf", i);
    global_dpd_->buf4_init(&CMnEf, PSIF_EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, lbl);
    global_dpd_->dot24(&Fme,&CMnEf,&SIA, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&CMnEf);
    global_dpd_->file2_close(&Fme);

    /* Sia += Fme*Cimae + FME*CiMaE */
    global_dpd_->file2_init(&Fme, PSIF_CC_OEI, H_IRR, 0, 1, "Fme");
    sprintf(lbl, "%s %d", "Cmnef", i);
    global_dpd_->buf4_init(&Cmnef, PSIF_EOM_Cmnef, C_irr, 0, 5, 2, 7, 0, lbl);
    global_dpd_->dot24(&Fme,&Cmnef,&Sia, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&Cmnef);
    global_dpd_->file2_close(&Fme);

    global_dpd_->file2_init(&FME, PSIF_CC_OEI, H_IRR, 0, 1, "FME");
    global_dpd_->buf4_init(&CmNeF, PSIF_EOM_TMP, C_irr, 0, 5, 0, 5, 0, "CmNeF");
    global_dpd_->dot24(&FME,&CmNeF,&Sia, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&CmNeF);
    global_dpd_->file2_close(&FME);

    global_dpd_->file2_close(&Sia);
    global_dpd_->file2_close(&SIA);
  }

  else { /* UHF */

    sprintf(lbl, "%s %d", "SIA", i);
    global_dpd_->file2_init(&SIA, PSIF_EOM_SIA, C_irr, 0, 1, lbl);
    sprintf(lbl, "%s %d", "Sia", i);
    global_dpd_->file2_init(&Sia, PSIF_EOM_Sia, C_irr, 2, 3, lbl);

    /* SIA += FME*CIMAE + Fme*CImAe */
    global_dpd_->file2_init(&FME, PSIF_CC_OEI, H_IRR, 0, 1, "FME");
    sprintf(lbl, "%s %d", "CMNEF", i);
    global_dpd_->buf4_init(&CMNEF, PSIF_EOM_CMNEF, C_irr, 0, 5, 2, 7, 0, lbl);
    global_dpd_->dot24(&FME,&CMNEF,&SIA, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&CMNEF);
    global_dpd_->file2_close(&FME);

    global_dpd_->file2_init(&Fme, PSIF_CC_OEI, H_IRR, 2, 3, "Fme");
    sprintf(lbl, "%s %d", "CMnEf", i);
    global_dpd_->buf4_init(&CMnEf, PSIF_EOM_CMnEf, C_irr, 22, 28, 22, 28, 0, lbl);
    global_dpd_->dot24(&Fme,&CMnEf,&SIA, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&CMnEf);
    global_dpd_->file2_close(&Fme);

    /* Sia += Fme*Cimae + FME*CiMaE */
    global_dpd_->file2_init(&Fme, PSIF_CC_OEI, H_IRR, 2, 3, "Fme");
    sprintf(lbl, "%s %d", "Cmnef", i);
    global_dpd_->buf4_init(&Cmnef, PSIF_EOM_Cmnef, C_irr, 10, 15, 12, 17, 0, lbl);
    global_dpd_->dot24(&Fme,&Cmnef,&Sia, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&Cmnef);
    global_dpd_->file2_close(&Fme);

    global_dpd_->file2_init(&FME, PSIF_CC_OEI, H_IRR, 0, 1, "FME");
    global_dpd_->buf4_init(&CmNeF, PSIF_EOM_TMP, C_irr, 23, 29, 23, 29, 0, "CmNeF");
    global_dpd_->dot24(&FME,&CmNeF,&Sia, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&CmNeF);
    global_dpd_->file2_close(&FME);

    global_dpd_->file2_close(&Sia);
    global_dpd_->file2_close(&SIA);
  }

#ifdef EOM_DEBUG
  check_sum("FSD    ",i,C_irr);
#endif

  return;
}

}} // namespace psi::cceom
