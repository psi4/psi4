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
	  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
	  dpd_buf4_copy(&D, CC_LAMBDA, "New LIjAb");
	  dpd_buf4_close(&D);
	}
	else if(params.ref == 1) { /** ROHF **/
	  dpd_buf4_init(&D, CC_DINTS, 0, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
	  dpd_buf4_copy(&D, CC_LAMBDA, "New LIJAB");
	  dpd_buf4_copy(&D, CC_LAMBDA, "New Lijab");
	  dpd_buf4_close(&D);
  
	  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
	  dpd_buf4_copy(&D, CC_LAMBDA, "New LIjAb");
	  dpd_buf4_close(&D);
	}
	else if(params.ref == 2) { /** UHF **/
	  dpd_buf4_init(&D, CC_DINTS, 0, 2, 7, 2, 7, 0, "D <IJ||AB> (I>J,A>B)");
	  dpd_buf4_copy(&D, CC_LAMBDA, "New LIJAB");
	  dpd_buf4_close(&D);
  
	  dpd_buf4_init(&D, CC_DINTS, 0, 12, 17, 12, 17, 0, "D <ij||ab> (i>j,a>b)");
	  dpd_buf4_copy(&D, CC_LAMBDA, "New Lijab");
	  dpd_buf4_close(&D);
  
	  dpd_buf4_init(&D, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
	  dpd_buf4_copy(&D, CC_LAMBDA, "New LIjAb");
	  dpd_buf4_close(&D);

	  /* If CCSD(T) gradient, add T3 contributions */
	  if(!strcmp(params.wfn,"CCSD_T")) {
	    dpd_buf4_init(&D, CC_MISC, 0, 2, 7, 2, 7, 0, "SIJAB");
	    dpd_buf4_init(&X2, CC_LAMBDA, 0, 2, 7, 2, 7, 0, "New LIJAB");
	    dpd_buf4_axpy(&D, &X2, 1);
	    dpd_buf4_close(&X2);
	    dpd_buf4_close(&D);

	    dpd_buf4_init(&D, CC_MISC, 0, 12, 17, 12, 17, 0, "Sijab");
	    dpd_buf4_init(&X2, CC_LAMBDA, 0, 12, 17, 12, 17, 0, "New Lijab");
	    dpd_buf4_axpy(&D, &X2, 1);
	    dpd_buf4_close(&X2);
	    dpd_buf4_close(&D);
	   
	    dpd_buf4_init(&D, CC_MISC, 0, 22, 28, 22, 28, 0, "SIjAb");
	    dpd_buf4_init(&X2, CC_LAMBDA, 0, 22, 28, 22, 28, 0, "New LIjAb");
	    dpd_buf4_axpy(&D, &X2, 1);
	    dpd_buf4_close(&X2);
	    dpd_buf4_close(&D);
	  }
	}
      }
      /* excited state - no inhomogeneous term, first term is E*L */
      else if (!params.zeta) {
	if (params.ref == 0) { /* RHF */
	  dpd_buf4_init(&D, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
	  dpd_buf4_init(&Dold, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
	  dpd_buf4_axpy(&Dold, &D, -1.0 * L_params.cceom_energy);
	  dpd_buf4_close(&Dold);
	  dpd_buf4_close(&D);
	}
	else if (params.ref == 1 ) { /* ROHF */
	  dpd_buf4_init(&D, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB");
	  dpd_buf4_init(&Dold, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
	  dpd_buf4_axpy(&Dold, &D, -1.0 * L_params.cceom_energy);
	  dpd_buf4_close(&Dold);
	  dpd_buf4_close(&D);
	  dpd_buf4_init(&D, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New Lijab");
	  dpd_buf4_init(&Dold, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "Lijab");
	  dpd_buf4_axpy(&Dold, &D, -1.0 * L_params.cceom_energy);
	  dpd_buf4_close(&Dold);
	  dpd_buf4_close(&D);
	  dpd_buf4_init(&D, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
	  dpd_buf4_init(&Dold, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
	  dpd_buf4_axpy(&Dold, &D, -1.0 * L_params.cceom_energy);
	  dpd_buf4_close(&Dold);
	  dpd_buf4_close(&D);
	}
	else { /** UHF **/
	  /* do nothing - TDC did not change to increments for the UHF case */
	}
      }
      /* solving zeta equations, homogeneous term is Xi, zero out files */
      else {
	if (params.ref == 0) { /* RHF */
	  dpd_buf4_init(&X2, EOM_XI, L_irr, 0, 5, 0, 5, 0, "XIjAb");
	  dpd_buf4_copy(&X2, CC_LAMBDA, "New LIjAb");
	  dpd_buf4_close(&X2);
	}
	else if (params.ref == 1 ) { /* ROHF */
	  dpd_buf4_init(&X2, EOM_XI, L_irr, 2, 7, 2, 7, 0, "XIJAB");
	  dpd_buf4_copy(&X2, CC_LAMBDA, "New LIJAB");
	  dpd_buf4_close(&X2);
	  dpd_buf4_init(&X2, EOM_XI, L_irr, 2, 7, 2, 7, 0, "Xijab");
	  dpd_buf4_copy(&X2, CC_LAMBDA, "New Lijab");
	  dpd_buf4_close(&X2);
	  dpd_buf4_init(&X2, EOM_XI, L_irr, 0, 5, 0, 5, 0, "XIjAb");
	  dpd_buf4_copy(&X2, CC_LAMBDA, "New LIjAb");
	  dpd_buf4_close(&X2);
	}
	else { /** UHF **/
	  dpd_buf4_init(&X2, EOM_XI, L_irr, 2, 7, 2, 7, 0, "XIJAB");
	  dpd_buf4_copy(&X2, CC_LAMBDA, "New LIJAB");
	  dpd_buf4_close(&X2);
	  dpd_buf4_init(&X2, EOM_XI, L_irr, 12, 17, 12, 17, 0, "Xijab");
	  dpd_buf4_copy(&X2, CC_LAMBDA, "New Lijab");
	  dpd_buf4_close(&X2);
	  dpd_buf4_init(&X2, EOM_XI, L_irr, 22, 28, 22, 28, 0, "XIjAb");
	  dpd_buf4_copy(&X2, CC_LAMBDA, "New LIjAb");
	  dpd_buf4_close(&X2);
	}
      }
    }

  }} // namespace psi::cclambda
