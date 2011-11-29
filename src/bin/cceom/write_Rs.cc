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
#include <libchkpt/chkpt.h>
#include <libpsio/psio.h>
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
        chkpt_init(PSIO_OPEN_OLD);
        if (!params.full_matrix) {
			    etot = evals[eom_params.prop_root]+moinfo.ecc+moinfo.eref;
		    } else {
			    etot = evals[eom_params.prop_root]+moinfo.eref;
			  }
        chkpt_wt_etot(etot);
        fprintf(outfile,"Energy written to chkpt:Etot %15.10lf\n", etot);
        chkpt_wt_statespi(eom_params.states_per_irrep);
        fprintf(outfile,"States per irrep written to chkpt.\n");
        chkpt_close();
			}
    }
    /* cclambda expects excitation energies */
    if (!params.full_matrix) {
      etot = evals[i];
    } 
    else {
      psio_read_entry(CC_HBAR, "Reference expectation value",
  		      (char *) &(expectation_val), sizeof(double));
      etot = evals[i] - expectation_val;
    }

    if(params.wfn == "EOM_CC2") {
      sprintf(E_lbl, "EOM CC2 Energy for root %d %d", C_irr, R_index);
      psio_write_entry(CC_INFO, E_lbl, (char *) &etot, sizeof(double));
    }
    else if(params.wfn == "EOM_CCSD") {
      sprintf(E_lbl, "EOM CCSD Energy for root %d %d", C_irr, R_index);
      psio_write_entry(CC_INFO, E_lbl, (char *) &etot, sizeof(double));
    }
    else if(params.wfn == "EOM_CC3") {
      sprintf(E_lbl, "EOM CC3 Energy for root %d %d", C_irr, R_index);
      psio_write_entry(CC_INFO, E_lbl, (char *) &etot, sizeof(double));
    }

    sprintf(C_lbl, "CME %d", i);
    sprintf(R_lbl, "RIA %d %d", C_irr, R_index);

    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, C_lbl);
    dpd_file2_copy(&CME, CC_RAMPS, R_lbl);
    dpd_file2_close(&CME);

    if (params.full_matrix) {
      sprintf(C_lbl, "C0 %d", i);
      psio_read_entry(EOM_CME, C_lbl, (char *) &C0, sizeof(double));
      sprintf(R_lbl, "R0 %d %d", C_irr, R_index);
      psio_write_entry(CC_RAMPS, R_lbl, (char *) &C0, sizeof(double));
    }

    sprintf(C_lbl, "CMnEf %d", i);
    sprintf(R_lbl, "RIjAb %d %d", C_irr, R_index);
    if (params.eom_ref <= 1)
      dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, C_lbl);
    else if (params.eom_ref == 2)
      dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 22, 28, 22, 28, 0, C_lbl);
    dpd_buf4_copy(&CMnEf, CC_RAMPS, R_lbl);
    dpd_buf4_close(&CMnEf);
          
    if(params.eom_ref > 0) {
      sprintf(C_lbl, "Cme %d", i);
      sprintf(R_lbl, "Ria %d %d", C_irr, R_index);
      if (params.eom_ref == 1)
         dpd_file2_init(&Cme, EOM_Cme, C_irr, 0, 1, C_lbl);
      else if (params.eom_ref == 2)
         dpd_file2_init(&Cme, EOM_Cme, C_irr, 2, 3, C_lbl);
      dpd_file2_copy(&Cme, CC_RAMPS, R_lbl);
      dpd_file2_close(&Cme);

      sprintf(C_lbl, "CMNEF %d", i);
      sprintf(R_lbl, "RIJAB %d %d", C_irr, R_index);

      dpd_buf4_init(&CMNEF, EOM_CMNEF, C_irr, 2, 7, 2, 7, 0, C_lbl);
      dpd_buf4_copy(&CMNEF, CC_RAMPS, R_lbl);
      dpd_buf4_close(&CMNEF);

      sprintf(C_lbl, "Cmnef %d", i);
      sprintf(R_lbl, "Rijab %d %d", C_irr, R_index);

      if (params.eom_ref == 1)
        dpd_buf4_init(&Cmnef, EOM_Cmnef, C_irr, 2, 7, 2, 7, 0, C_lbl);
      else if (params.eom_ref ==2)
        dpd_buf4_init(&Cmnef, EOM_Cmnef, C_irr, 12, 17, 12, 17, 0, C_lbl);
        dpd_buf4_copy(&Cmnef, CC_RAMPS, R_lbl);
        dpd_buf4_close(&Cmnef);
    }
  }
}

}} // namespace psi::cceom
