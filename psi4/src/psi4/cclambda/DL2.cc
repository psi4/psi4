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
#include <string.h>
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cclambda {

    void DL2(struct L_Params L_params)
    {
      dpdbuf4 D, Dold, X2;
      int L_irr;
      L_irr = L_params.irrep;

      if (L_params.ground) {
	/* RHS = <ij||ab> */
	if(params.ref == 0) { /** RHF **/
	  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
	  global_dpd_->buf4_copy(&D, PSIF_CC_LAMBDA, "New LIjAb");
	  global_dpd_->buf4_close(&D);
         // Add T3 contribution to CCSD(T) lambda equations
          if(params.wfn == "CCSD_T") {
            global_dpd_->buf4_init(&D, PSIF_CC_MISC, 0, 0, 5, 0, 5, 0, "SIjAb(T)");
            global_dpd_->buf4_init(&X2, PSIF_CC_LAMBDA, 0, 0, 5, 0, 5, 0, "New LIjAb");
            global_dpd_->buf4_axpy(&D, &X2, 1.0);
            global_dpd_->buf4_close(&X2);
            global_dpd_->buf4_close(&D);
           }
	}
	else if(params.ref == 1) { /** ROHF **/
	  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
	  global_dpd_->buf4_copy(&D, PSIF_CC_LAMBDA, "New LIJAB");
	  global_dpd_->buf4_copy(&D, PSIF_CC_LAMBDA, "New Lijab");
	  global_dpd_->buf4_close(&D);

	  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
	  global_dpd_->buf4_copy(&D, PSIF_CC_LAMBDA, "New LIjAb");
	  global_dpd_->buf4_close(&D);
	}
	else if(params.ref == 2) { /** UHF **/
	  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 2, 7, 2, 7, 0, "D <IJ||AB> (I>J,A>B)");
	  global_dpd_->buf4_copy(&D, PSIF_CC_LAMBDA, "New LIJAB");
	  global_dpd_->buf4_close(&D);

	  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 12, 17, 12, 17, 0, "D <ij||ab> (i>j,a>b)");
	  global_dpd_->buf4_copy(&D, PSIF_CC_LAMBDA, "New Lijab");
	  global_dpd_->buf4_close(&D);

	  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
	  global_dpd_->buf4_copy(&D, PSIF_CC_LAMBDA, "New LIjAb");
	  global_dpd_->buf4_close(&D);

	  /* If CCSD(T) gradient, add T3 contributions */
	  if(params.wfn == "CCSD_T") {
	    global_dpd_->buf4_init(&D, PSIF_CC_MISC, 0, 2, 7, 2, 7, 0, "SIJAB");
	    global_dpd_->buf4_init(&X2, PSIF_CC_LAMBDA, 0, 2, 7, 2, 7, 0, "New LIJAB");
	    global_dpd_->buf4_axpy(&D, &X2, 1);
	    global_dpd_->buf4_close(&X2);
	    global_dpd_->buf4_close(&D);

	    global_dpd_->buf4_init(&D, PSIF_CC_MISC, 0, 12, 17, 12, 17, 0, "Sijab");
	    global_dpd_->buf4_init(&X2, PSIF_CC_LAMBDA, 0, 12, 17, 12, 17, 0, "New Lijab");
	    global_dpd_->buf4_axpy(&D, &X2, 1);
	    global_dpd_->buf4_close(&X2);
	    global_dpd_->buf4_close(&D);

	    global_dpd_->buf4_init(&D, PSIF_CC_MISC, 0, 22, 28, 22, 28, 0, "SIjAb");
	    global_dpd_->buf4_init(&X2, PSIF_CC_LAMBDA, 0, 22, 28, 22, 28, 0, "New LIjAb");
	    global_dpd_->buf4_axpy(&D, &X2, 1);
	    global_dpd_->buf4_close(&X2);
	    global_dpd_->buf4_close(&D);
	  }
	}
      }
      /* excited state - no inhomogeneous term, first term is E*L */
      else if (!params.zeta) {
	if (params.ref == 0) { /* RHF */
	  global_dpd_->buf4_init(&D, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
	  global_dpd_->buf4_init(&Dold, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
	  global_dpd_->buf4_axpy(&Dold, &D, -1.0 * L_params.cceom_energy);
	  global_dpd_->buf4_close(&Dold);
	  global_dpd_->buf4_close(&D);
	}
	else if (params.ref == 1 ) { /* ROHF */
	  global_dpd_->buf4_init(&D, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB");
	  global_dpd_->buf4_init(&Dold, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
	  global_dpd_->buf4_axpy(&Dold, &D, -1.0 * L_params.cceom_energy);
	  global_dpd_->buf4_close(&Dold);
	  global_dpd_->buf4_close(&D);
	  global_dpd_->buf4_init(&D, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New Lijab");
	  global_dpd_->buf4_init(&Dold, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "Lijab");
	  global_dpd_->buf4_axpy(&Dold, &D, -1.0 * L_params.cceom_energy);
	  global_dpd_->buf4_close(&Dold);
	  global_dpd_->buf4_close(&D);
	  global_dpd_->buf4_init(&D, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
	  global_dpd_->buf4_init(&Dold, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
	  global_dpd_->buf4_axpy(&Dold, &D, -1.0 * L_params.cceom_energy);
	  global_dpd_->buf4_close(&Dold);
	  global_dpd_->buf4_close(&D);
	}
	else { /** UHF **/
	  /* do nothing - TDC did not change to increments for the UHF case */
	}
      }
      /* solving zeta equations, homogeneous term is Xi, zero out files */
      else {
	if (params.ref == 0) { /* RHF */
	  global_dpd_->buf4_init(&X2, PSIF_EOM_XI, L_irr, 0, 5, 0, 5, 0, "XIjAb");
	  global_dpd_->buf4_copy(&X2, PSIF_CC_LAMBDA, "New LIjAb");
	  global_dpd_->buf4_close(&X2);
	}
	else if (params.ref == 1 ) { /* ROHF */
	  global_dpd_->buf4_init(&X2, PSIF_EOM_XI, L_irr, 2, 7, 2, 7, 0, "XIJAB");
	  global_dpd_->buf4_copy(&X2, PSIF_CC_LAMBDA, "New LIJAB");
	  global_dpd_->buf4_close(&X2);
	  global_dpd_->buf4_init(&X2, PSIF_EOM_XI, L_irr, 2, 7, 2, 7, 0, "Xijab");
	  global_dpd_->buf4_copy(&X2, PSIF_CC_LAMBDA, "New Lijab");
	  global_dpd_->buf4_close(&X2);
	  global_dpd_->buf4_init(&X2, PSIF_EOM_XI, L_irr, 0, 5, 0, 5, 0, "XIjAb");
	  global_dpd_->buf4_copy(&X2, PSIF_CC_LAMBDA, "New LIjAb");
	  global_dpd_->buf4_close(&X2);
	}
	else { /** UHF **/
	  global_dpd_->buf4_init(&X2, PSIF_EOM_XI, L_irr, 2, 7, 2, 7, 0, "XIJAB");
	  global_dpd_->buf4_copy(&X2, PSIF_CC_LAMBDA, "New LIJAB");
	  global_dpd_->buf4_close(&X2);
	  global_dpd_->buf4_init(&X2, PSIF_EOM_XI, L_irr, 12, 17, 12, 17, 0, "Xijab");
	  global_dpd_->buf4_copy(&X2, PSIF_CC_LAMBDA, "New Lijab");
	  global_dpd_->buf4_close(&X2);
	  global_dpd_->buf4_init(&X2, PSIF_EOM_XI, L_irr, 22, 28, 22, 28, 0, "XIjAb");
	  global_dpd_->buf4_copy(&X2, PSIF_CC_LAMBDA, "New LIjAb");
	  global_dpd_->buf4_close(&X2);
	}
      }
    }

  }} // namespace psi::cclambda
