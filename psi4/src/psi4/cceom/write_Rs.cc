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
 *   write_Rs writes out all of the converged R's to RAMPS for the current irrep
 */

#include <cstdio>
#include <string>
#include <cmath>
#include "psi4/libpsio/psio.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cceom {

void write_Rs(int C_irr, double *evals, int *converged) {
  int i;
  dpdfile2 CME, Cme;
  dpdbuf4 CMNEF, Cmnef, CMnEf;
  int R_index = -1;
  char C_lbl[32], R_lbl[32], E_lbl[32];
  double etot, expectation_val, C0;

  for (i=0; i<eom_params.cs_per_irrep[C_irr]; ++i) {
    if (!converged[i]) continue; /* this root did not converge */
    ++R_index;

    if (C_irr == eom_params.prop_sym) {
      if (i == eom_params.prop_root) {
        if (!params.full_matrix) {
          etot = evals[eom_params.prop_root]+moinfo.ecc+moinfo.eref;
        }
        else {
          etot = evals[eom_params.prop_root]+moinfo.eref;
        }
        psio_write_entry(PSIF_CC_INFO, "Total energy", (char *) &etot, sizeof(double));
        outfile->Printf("Energy written to CC_INFO:Etot %15.10lf\n", etot);
        psio_write_entry(PSIF_CC_INFO, "States per irrep", (char *) eom_params.states_per_irrep, moinfo.nirreps * sizeof(int));
        outfile->Printf("States per irrep written to CC_INFO.\n");
      }
    }
    /* cclambda expects excitation energies */
    if (!params.full_matrix) {
      etot = evals[i];
    }
    else {
      psio_read_entry(PSIF_CC_HBAR, "Reference expectation value",
  		      (char *) &(expectation_val), sizeof(double));
      etot = evals[i] - expectation_val;
    }

    if(params.wfn == "EOM_CC2") {
      sprintf(E_lbl, "EOM CC2 Energy for root %d %d", C_irr, R_index);
      psio_write_entry(PSIF_CC_INFO, E_lbl, (char *) &etot, sizeof(double));
    }
    else if(params.wfn == "EOM_CCSD") {
      sprintf(E_lbl, "EOM CCSD Energy for root %d %d", C_irr, R_index);
      psio_write_entry(PSIF_CC_INFO, E_lbl, (char *) &etot, sizeof(double));
    }
    else if(params.wfn == "EOM_CC3") {
      sprintf(E_lbl, "EOM CC3 Energy for root %d %d", C_irr, R_index);
      psio_write_entry(PSIF_CC_INFO, E_lbl, (char *) &etot, sizeof(double));
    }

    sprintf(C_lbl, "CME %d", i);
    sprintf(R_lbl, "RIA %d %d", C_irr, R_index);

    global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, C_lbl);
    global_dpd_->file2_copy(&CME, PSIF_CC_RAMPS, R_lbl);
    global_dpd_->file2_close(&CME);

    if (params.full_matrix) {
      sprintf(C_lbl, "C0 %d", i);
      psio_read_entry(PSIF_EOM_CME, C_lbl, (char *) &C0, sizeof(double));
      sprintf(R_lbl, "R0 %d %d", C_irr, R_index);
      psio_write_entry(PSIF_CC_RAMPS, R_lbl, (char *) &C0, sizeof(double));
    }

    sprintf(C_lbl, "CMnEf %d", i);
    sprintf(R_lbl, "RIjAb %d %d", C_irr, R_index);
    if (params.eom_ref <= 1)
      global_dpd_->buf4_init(&CMnEf, PSIF_EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, C_lbl);
    else if (params.eom_ref == 2)
      global_dpd_->buf4_init(&CMnEf, PSIF_EOM_CMnEf, C_irr, 22, 28, 22, 28, 0, C_lbl);
    global_dpd_->buf4_copy(&CMnEf, PSIF_CC_RAMPS, R_lbl);
    global_dpd_->buf4_close(&CMnEf);

    if(params.eom_ref > 0) {
      sprintf(C_lbl, "Cme %d", i);
      sprintf(R_lbl, "Ria %d %d", C_irr, R_index);
      if (params.eom_ref == 1)
         global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 0, 1, C_lbl);
      else if (params.eom_ref == 2)
         global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 2, 3, C_lbl);
      global_dpd_->file2_copy(&Cme, PSIF_CC_RAMPS, R_lbl);
      global_dpd_->file2_close(&Cme);

      sprintf(C_lbl, "CMNEF %d", i);
      sprintf(R_lbl, "RIJAB %d %d", C_irr, R_index);

      global_dpd_->buf4_init(&CMNEF, PSIF_EOM_CMNEF, C_irr, 2, 7, 2, 7, 0, C_lbl);
      global_dpd_->buf4_copy(&CMNEF, PSIF_CC_RAMPS, R_lbl);
      global_dpd_->buf4_close(&CMNEF);

      sprintf(C_lbl, "Cmnef %d", i);
      sprintf(R_lbl, "Rijab %d %d", C_irr, R_index);

      if (params.eom_ref == 1)
        global_dpd_->buf4_init(&Cmnef, PSIF_EOM_Cmnef, C_irr, 2, 7, 2, 7, 0, C_lbl);
      else if (params.eom_ref ==2)
        global_dpd_->buf4_init(&Cmnef, PSIF_EOM_Cmnef, C_irr, 12, 17, 12, 17, 0, C_lbl);
        global_dpd_->buf4_copy(&Cmnef, PSIF_CC_RAMPS, R_lbl);
        global_dpd_->buf4_close(&Cmnef);
    }
  }
}

}} // namespace psi::cceom
