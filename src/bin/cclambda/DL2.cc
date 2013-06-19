/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

/*! \file
    \ingroup CCLAMBDA
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <string.h>
#include <libdpd/dpd.h>
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
	  dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
	  dpd_->buf4_copy(&D, PSIF_CC_LAMBDA, "New LIjAb");
	  dpd_->buf4_close(&D);
	}
	else if(params.ref == 1) { /** ROHF **/
	  dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
	  dpd_->buf4_copy(&D, PSIF_CC_LAMBDA, "New LIJAB");
	  dpd_->buf4_copy(&D, PSIF_CC_LAMBDA, "New Lijab");
	  dpd_->buf4_close(&D);
  
	  dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
	  dpd_->buf4_copy(&D, PSIF_CC_LAMBDA, "New LIjAb");
	  dpd_->buf4_close(&D);
	}
	else if(params.ref == 2) { /** UHF **/
	  dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 2, 7, 2, 7, 0, "D <IJ||AB> (I>J,A>B)");
	  dpd_->buf4_copy(&D, PSIF_CC_LAMBDA, "New LIJAB");
	  dpd_->buf4_close(&D);
  
	  dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 12, 17, 12, 17, 0, "D <ij||ab> (i>j,a>b)");
	  dpd_->buf4_copy(&D, PSIF_CC_LAMBDA, "New Lijab");
	  dpd_->buf4_close(&D);
  
	  dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
	  dpd_->buf4_copy(&D, PSIF_CC_LAMBDA, "New LIjAb");
	  dpd_->buf4_close(&D);

	  /* If CCSD(T) gradient, add T3 contributions */
	  if(params.wfn == "CCSD_T") {
	    dpd_->buf4_init(&D, PSIF_CC_MISC, 0, 2, 7, 2, 7, 0, "SIJAB");
	    dpd_->buf4_init(&X2, PSIF_CC_LAMBDA, 0, 2, 7, 2, 7, 0, "New LIJAB");
	    dpd_->buf4_axpy(&D, &X2, 1);
	    dpd_->buf4_close(&X2);
	    dpd_->buf4_close(&D);

	    dpd_->buf4_init(&D, PSIF_CC_MISC, 0, 12, 17, 12, 17, 0, "Sijab");
	    dpd_->buf4_init(&X2, PSIF_CC_LAMBDA, 0, 12, 17, 12, 17, 0, "New Lijab");
	    dpd_->buf4_axpy(&D, &X2, 1);
	    dpd_->buf4_close(&X2);
	    dpd_->buf4_close(&D);
	   
	    dpd_->buf4_init(&D, PSIF_CC_MISC, 0, 22, 28, 22, 28, 0, "SIjAb");
	    dpd_->buf4_init(&X2, PSIF_CC_LAMBDA, 0, 22, 28, 22, 28, 0, "New LIjAb");
	    dpd_->buf4_axpy(&D, &X2, 1);
	    dpd_->buf4_close(&X2);
	    dpd_->buf4_close(&D);
	  }
	}
      }
      /* excited state - no inhomogeneous term, first term is E*L */
      else if (!params.zeta) {
	if (params.ref == 0) { /* RHF */
	  dpd_->buf4_init(&D, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
	  dpd_->buf4_init(&Dold, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
	  dpd_->buf4_axpy(&Dold, &D, -1.0 * L_params.cceom_energy);
	  dpd_->buf4_close(&Dold);
	  dpd_->buf4_close(&D);
	}
	else if (params.ref == 1 ) { /* ROHF */
	  dpd_->buf4_init(&D, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB");
	  dpd_->buf4_init(&Dold, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
	  dpd_->buf4_axpy(&Dold, &D, -1.0 * L_params.cceom_energy);
	  dpd_->buf4_close(&Dold);
	  dpd_->buf4_close(&D);
	  dpd_->buf4_init(&D, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New Lijab");
	  dpd_->buf4_init(&Dold, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "Lijab");
	  dpd_->buf4_axpy(&Dold, &D, -1.0 * L_params.cceom_energy);
	  dpd_->buf4_close(&Dold);
	  dpd_->buf4_close(&D);
	  dpd_->buf4_init(&D, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
	  dpd_->buf4_init(&Dold, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
	  dpd_->buf4_axpy(&Dold, &D, -1.0 * L_params.cceom_energy);
	  dpd_->buf4_close(&Dold);
	  dpd_->buf4_close(&D);
	}
	else { /** UHF **/
	  /* do nothing - TDC did not change to increments for the UHF case */
	}
      }
      /* solving zeta equations, homogeneous term is Xi, zero out files */
      else {
	if (params.ref == 0) { /* RHF */
	  dpd_->buf4_init(&X2, PSIF_EOM_XI, L_irr, 0, 5, 0, 5, 0, "XIjAb");
	  dpd_->buf4_copy(&X2, PSIF_CC_LAMBDA, "New LIjAb");
	  dpd_->buf4_close(&X2);
	}
	else if (params.ref == 1 ) { /* ROHF */
	  dpd_->buf4_init(&X2, PSIF_EOM_XI, L_irr, 2, 7, 2, 7, 0, "XIJAB");
	  dpd_->buf4_copy(&X2, PSIF_CC_LAMBDA, "New LIJAB");
	  dpd_->buf4_close(&X2);
	  dpd_->buf4_init(&X2, PSIF_EOM_XI, L_irr, 2, 7, 2, 7, 0, "Xijab");
	  dpd_->buf4_copy(&X2, PSIF_CC_LAMBDA, "New Lijab");
	  dpd_->buf4_close(&X2);
	  dpd_->buf4_init(&X2, PSIF_EOM_XI, L_irr, 0, 5, 0, 5, 0, "XIjAb");
	  dpd_->buf4_copy(&X2, PSIF_CC_LAMBDA, "New LIjAb");
	  dpd_->buf4_close(&X2);
	}
	else { /** UHF **/
	  dpd_->buf4_init(&X2, PSIF_EOM_XI, L_irr, 2, 7, 2, 7, 0, "XIJAB");
	  dpd_->buf4_copy(&X2, PSIF_CC_LAMBDA, "New LIJAB");
	  dpd_->buf4_close(&X2);
	  dpd_->buf4_init(&X2, PSIF_EOM_XI, L_irr, 12, 17, 12, 17, 0, "Xijab");
	  dpd_->buf4_copy(&X2, PSIF_CC_LAMBDA, "New Lijab");
	  dpd_->buf4_close(&X2);
	  dpd_->buf4_init(&X2, PSIF_EOM_XI, L_irr, 22, 28, 22, 28, 0, "XIjAb");
	  dpd_->buf4_copy(&X2, PSIF_CC_LAMBDA, "New LIjAb");
	  dpd_->buf4_close(&X2);
	}
      }
    }

  }} // namespace psi::cclambda
