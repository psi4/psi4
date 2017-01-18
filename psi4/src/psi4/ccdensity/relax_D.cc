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

/* relax_D(): Add the orbital-response contributions to the Dia block
** of the one-electron density matrix:
**
** D(A,I) = D(amp)(A,I) + D(orb)(A,I)
**
** D(I,A) = D(amp)(I,A) + D(orb)(A,I)
**
** */

void relax_D(struct RHO_Params rho_params)
{
  dpdfile2 D1, D2, I1, I2;

  if(params.ref == 0) {
    global_dpd_->file2_init(&D1, PSIF_CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
    global_dpd_->file2_init(&D2, PSIF_CC_OEI, 0, 1, 0, "D(orb)(A,I)");
    global_dpd_->file2_axpy(&D2, &D1, 1.0, 1);
    global_dpd_->file2_close(&D2);
    global_dpd_->file2_close(&D1);

    global_dpd_->file2_init(&D1, PSIF_CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
    global_dpd_->file2_init(&D2, PSIF_CC_OEI, 0, 1, 0, "D(orb)(A,I)");
    global_dpd_->file2_axpy(&D2, &D1, 1.0, 1);
    global_dpd_->file2_close(&D2);
    global_dpd_->file2_close(&D1);
    /* Add the contributions of dependent pairs (i,j) and (a,b) to the density due 
     * to the use of canonical perturbed orbitals:  
     * Dij += (I'ij - I'ji)/(fii - fjj) :  Dab += (I'ab - I'ba)/(faa - fbb)  */
    if (params.wfn == "CCSD_T" && params.dertype ==1){
    global_dpd_->file2_init(&D1, PSIF_CC_OEI, 0, 0, 0, rho_params.DIJ_lbl);
    global_dpd_->file2_init(&I1, PSIF_CC_TMP, 0, 0, 0, "delta_I/delta_f_IJ");
    global_dpd_->file2_init(&D2, PSIF_CC_OEI, 0, 1, 1, rho_params.DAB_lbl);
    global_dpd_->file2_init(&I2, PSIF_CC_TMP, 0, 1, 1, "delta_I/delta_f_AB");
    global_dpd_->file2_axpy(&I1, &D1, 1.0, 0);
    global_dpd_->file2_axpy(&I2, &D2, 1.0, 0);
    global_dpd_->file2_close(&D1);
    global_dpd_->file2_close(&D2);
    global_dpd_->file2_close(&I1);
    global_dpd_->file2_close(&I2);
   }
 }
  else if(params.ref == 1) {
    global_dpd_->file2_init(&D1, PSIF_CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
    global_dpd_->file2_init(&D2, PSIF_CC_OEI, 0, 1, 0, "D(orb)(A,I)");
    global_dpd_->file2_axpy(&D2, &D1, 1.0, 1);
    global_dpd_->file2_close(&D2);
    global_dpd_->file2_close(&D1);

    global_dpd_->file2_init(&D1, PSIF_CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
    global_dpd_->file2_init(&D2, PSIF_CC_OEI, 0, 1, 0, "D(orb)(A,I)");
    global_dpd_->file2_axpy(&D2, &D1, 1.0, 1);
    global_dpd_->file2_close(&D2);
    global_dpd_->file2_close(&D1);

    global_dpd_->file2_init(&D1, PSIF_CC_OEI, 0, 0, 1, rho_params.Dai_lbl);
    global_dpd_->file2_init(&D2, PSIF_CC_OEI, 0, 1, 0, "D(orb)(a,i)");
    global_dpd_->file2_axpy(&D2, &D1, 1.0, 1);
    global_dpd_->file2_close(&D2);
    global_dpd_->file2_close(&D1);

    global_dpd_->file2_init(&D1, PSIF_CC_OEI, 0, 0, 1, rho_params.Dia_lbl);
    global_dpd_->file2_init(&D2, PSIF_CC_OEI, 0, 1, 0, "D(orb)(a,i)");
    global_dpd_->file2_axpy(&D2, &D1, 1.0, 1);
    global_dpd_->file2_close(&D2);
    global_dpd_->file2_close(&D1);
  }
  else if(params.ref == 2) {

    global_dpd_->file2_init(&D1, PSIF_CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
    global_dpd_->file2_init(&D2, PSIF_CC_OEI, 0, 1, 0, "D(orb)(A,I)");
    global_dpd_->file2_axpy(&D2, &D1, 1.0, 1);
    global_dpd_->file2_close(&D2);
    global_dpd_->file2_close(&D1);

    global_dpd_->file2_init(&D1, PSIF_CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
    global_dpd_->file2_init(&D2, PSIF_CC_OEI, 0, 1, 0, "D(orb)(A,I)");
    global_dpd_->file2_axpy(&D2, &D1, 1.0, 1);
    global_dpd_->file2_close(&D2);
    global_dpd_->file2_close(&D1);

    global_dpd_->file2_init(&D1, PSIF_CC_OEI, 0, 2, 3, rho_params.Dai_lbl);
    global_dpd_->file2_init(&D2, PSIF_CC_OEI, 0, 3, 2, "D(orb)(a,i)");
    global_dpd_->file2_axpy(&D2, &D1, 1.0, 1);
    global_dpd_->file2_close(&D2);
    global_dpd_->file2_close(&D1);

    global_dpd_->file2_init(&D1, PSIF_CC_OEI, 0, 2, 3, rho_params.Dia_lbl);
    global_dpd_->file2_init(&D2, PSIF_CC_OEI, 0, 3, 2, "D(orb)(a,i)");
    global_dpd_->file2_axpy(&D2, &D1, 1.0, 1);
    global_dpd_->file2_close(&D2);
    global_dpd_->file2_close(&D1);


  }
}

}} // namespace psi::ccdensity
