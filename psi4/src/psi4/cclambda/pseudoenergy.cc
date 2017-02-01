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
    \ingroup CCLAMBDA
    \brief Enter brief description of file here
*/
#include <cstdio>
#include "psi4/libdpd/dpd.h"
#include <cmath>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cclambda {

double pseudoenergy(struct L_Params L_params)
{
  double LIJAB_energy, Lijab_energy, LIjAb_energy;
  double LIA_energy=0.0, Lia_energy=0.0, tval;
  dpdbuf4 LIJAB, Lijab, LIjAb, D;
  dpdfile2 Lia, LIA, Fme, FME;
  int L_irr;
  L_irr = L_params.irrep;

  if ( L_params.ground || ((L_params.irrep ==0)&&(fabs(L_params.R0)>1e-10)) ) {
    if(params.ref == 0) { /** RHF **/

      Lia_energy = 0.0;
      global_dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");
      global_dpd_->file2_init(&LIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "LIA");
      LIA_energy = global_dpd_->file2_dot(&FME,&LIA);
      global_dpd_->file2_close(&LIA);
      global_dpd_->file2_close(&FME);

      LIJAB_energy = 0.0;
      Lijab_energy = 0.0;
      global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
      global_dpd_->buf4_init(&LIjAb, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
      LIjAb_energy = global_dpd_->buf4_dot(&D, &LIjAb);
      global_dpd_->buf4_close(&LIjAb);
      global_dpd_->buf4_close(&D);
    }
    else if(params.ref == 1) { /** ROHF **/
      global_dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 0, 1, "Fme");
      global_dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");
      global_dpd_->file2_init(&Lia, PSIF_CC_LAMBDA, L_irr, 0, 1, "Lia");
      global_dpd_->file2_init(&LIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "LIA");

      LIA_energy = global_dpd_->file2_dot(&FME,&LIA);
      Lia_energy = global_dpd_->file2_dot(&Fme,&Lia);

      global_dpd_->file2_close(&Lia);
      global_dpd_->file2_close(&LIA);
      global_dpd_->file2_close(&Fme);
      global_dpd_->file2_close(&FME);

      global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
      global_dpd_->buf4_init(&LIJAB, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
      LIJAB_energy = global_dpd_->buf4_dot(&D, &LIJAB);
      global_dpd_->buf4_close(&LIJAB);
      global_dpd_->buf4_init(&Lijab, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "Lijab");
      Lijab_energy = global_dpd_->buf4_dot(&D, &Lijab);
      global_dpd_->buf4_close(&Lijab);
      global_dpd_->buf4_close(&D);

      global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
      global_dpd_->buf4_init(&LIjAb, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
      LIjAb_energy = global_dpd_->buf4_dot(&D, &LIjAb);
      global_dpd_->buf4_close(&LIjAb);
      global_dpd_->buf4_close(&D);
    }
    else if(params.ref == 2) { /** UHF **/

      global_dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 2, 3, "Fme");
      global_dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");
      global_dpd_->file2_init(&Lia, PSIF_CC_LAMBDA, L_irr, 2, 3, "Lia");
      global_dpd_->file2_init(&LIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "LIA");

      LIA_energy = global_dpd_->file2_dot(&FME,&LIA);
      Lia_energy = global_dpd_->file2_dot(&Fme,&Lia);

      global_dpd_->file2_close(&Lia);
      global_dpd_->file2_close(&LIA);
      global_dpd_->file2_close(&Fme);
      global_dpd_->file2_close(&FME);

      global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 2, 7, 2, 7, 0, "D <IJ||AB> (I>J,A>B)");
      global_dpd_->buf4_init(&LIJAB, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
      LIJAB_energy = global_dpd_->buf4_dot(&D, &LIJAB);
      global_dpd_->buf4_close(&LIJAB);
      global_dpd_->buf4_close(&D);

      global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 12, 17, 12, 17, 0, "D <ij||ab> (i>j,a>b)");
      global_dpd_->buf4_init(&Lijab, PSIF_CC_LAMBDA, L_irr, 12, 17, 12, 17, 0, "Lijab");
      Lijab_energy = global_dpd_->buf4_dot(&D, &Lijab);
      global_dpd_->buf4_close(&Lijab);
      global_dpd_->buf4_close(&D);

      global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
      global_dpd_->buf4_init(&LIjAb, PSIF_CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "LIjAb");
      LIjAb_energy = global_dpd_->buf4_dot(&D, &LIjAb);
      global_dpd_->buf4_close(&LIjAb);
      global_dpd_->buf4_close(&D);
    }
    /*
      outfile->Printf( "One A Energy = %20.14f\n", LIA_energy);
      outfile->Printf( "One B Energy = %20.14f\n", Lia_energy);
      outfile->Printf( "Two AA Energy = %20.14f\n", LIJAB_energy);
      outfile->Printf( "Two BB Energy = %20.14f\n", Lijab_energy);
      outfile->Printf( "Two AB Energy = %20.14f\n", LIjAb_energy);
    */
    return (LIJAB_energy + Lijab_energy + LIjAb_energy);
  }
  else { /* since pseudoenergy is 0 lets compute norm instead */
    if (params.ref <= 1) { /* RHF or ROHF */
      global_dpd_->file2_init(&Lia, PSIF_CC_LAMBDA, L_irr, 0, 1, "Lia");
      global_dpd_->file2_init(&LIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "LIA");
      LIA_energy = global_dpd_->file2_dot_self(&LIA);
      Lia_energy = global_dpd_->file2_dot_self(&Lia);
      global_dpd_->file2_close(&Lia);
      global_dpd_->file2_close(&LIA);
      global_dpd_->buf4_init(&LIJAB, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
      LIJAB_energy = global_dpd_->buf4_dot_self(&LIJAB);
      global_dpd_->buf4_close(&LIJAB);
      global_dpd_->buf4_init(&Lijab, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "Lijab");
      Lijab_energy = global_dpd_->buf4_dot_self(&Lijab);
      global_dpd_->buf4_close(&Lijab);
      global_dpd_->buf4_init(&LIjAb, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
      LIjAb_energy = global_dpd_->buf4_dot_self(&LIjAb);
      global_dpd_->buf4_close(&LIjAb);
      tval = LIA_energy + Lia_energy + LIJAB_energy + Lijab_energy + LIjAb_energy;
      tval = sqrt(tval);
      return tval;
    }
    else if (params.ref == 2) { /* UHF */
      global_dpd_->file2_init(&LIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "LIA");
      global_dpd_->file2_init(&Lia, PSIF_CC_LAMBDA, L_irr, 2, 3, "Lia");
      LIA_energy = global_dpd_->file2_dot_self(&LIA);
      Lia_energy = global_dpd_->file2_dot_self(&Lia);
      global_dpd_->file2_close(&Lia);
      global_dpd_->file2_close(&LIA);
      global_dpd_->buf4_init(&LIJAB, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
      LIJAB_energy = global_dpd_->buf4_dot_self(&LIJAB);
      global_dpd_->buf4_close(&LIJAB);
      global_dpd_->buf4_init(&Lijab, PSIF_CC_LAMBDA, L_irr, 12, 17, 12, 17, 0, "Lijab");
      Lijab_energy = global_dpd_->buf4_dot_self(&Lijab);
      global_dpd_->buf4_close(&Lijab);
      global_dpd_->buf4_init(&LIjAb, PSIF_CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "LIjAb");
      LIjAb_energy = global_dpd_->buf4_dot_self(&LIjAb);
      global_dpd_->buf4_close(&LIjAb);
      tval = LIA_energy + Lia_energy + LIJAB_energy + Lijab_energy + LIjAb_energy;
      tval = sqrt(tval);
      return tval;
    }
  }

  return 0.0;
}

}} // namespace psi::cclambda
