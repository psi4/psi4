/*! \file
    \ingroup CCENERGY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccenergy {

void Z_build(void)
{
  dpdbuf4 ZIJMA, Zijma, ZIjMa, ZIjmA, ZIjAm, ZMaIj, ZmAIj, Z;
  dpdbuf4 tauIJAB, tauijab, tauIjAb, tauIjbA, F_anti, F, tau;
  int Gmb, Gij, mb, nrows, ncols;

#ifdef TIME_CCENERGY
  timer_on("Z");
#endif

  if(params.ref == 0) { /** RHF **/
    /* ZMbIj = <Mb|Ef> * tau(Ij,Ef) */
    /* OOC code added 3/23/05  -TDC */
    dpd_buf4_init(&Z, CC_MISC, 0, 10, 0, 10, 0, 0, "ZMbIj");
    dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    dpd_buf4_init(&tau, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
    dpd_contract444(&F, &tau, &Z, 0, 0, 1, 0);
/*     for(Gmb=0; Gmb < moinfo.nirreps; Gmb++) { */
/*       Gij = Gmb;  /\* tau is totally symmetric *\/ */
/*       dpd_buf4_mat_irrep_init(&tau, Gij); */
/*       dpd_buf4_mat_irrep_rd(&tau, Gij); */

/*       dpd_buf4_mat_irrep_init(&Z, Gmb); */

/*       dpd_buf4_mat_irrep_row_init(&F, Gmb); */

/*       for(mb=0; mb < F.params->rowtot[Gmb]; mb++) { */
/* 	dpd_buf4_mat_irrep_row_rd(&F, Gmb, mb); */

/* 	nrows = tau.params->rowtot[Gij]; */
/* 	ncols = tau.params->coltot[Gij]; */
/* 	if(nrows && ncols) */
/* 	  C_DGEMV('n',nrows,ncols,1.0,tau.matrix[Gij][0],ncols,F.matrix[Gmb][0],1, */
/* 		  0.0,Z.matrix[Gmb][mb],1); */
/*       } */

/*       dpd_buf4_mat_irrep_row_close(&F, Gmb); */
/*       dpd_buf4_mat_irrep_wrt(&Z, Gmb); */
/*       dpd_buf4_mat_irrep_close(&Z, Gmb); */
/*       dpd_buf4_mat_irrep_close(&tau, Gij); */
/*     } */
    dpd_buf4_close(&tau);
    dpd_buf4_close(&F);
    dpd_buf4_close(&Z);  
  }
  else if(params.ref == 1) { /** ROHF **/
    dpd_buf4_init(&ZIJMA, CC_MISC, 0, 2, 10, 2, 10, 0, "ZIJMA");
    dpd_buf4_init(&Zijma, CC_MISC, 0, 2, 10, 2, 10, 0, "Zijma");
    dpd_buf4_init(&ZIjMa, CC_MISC, 0, 0, 10, 0, 10, 0, "ZIjMa");
    dpd_buf4_init(&ZIjmA, CC_MISC, 0, 0, 10, 0, 10, 0, "ZIjmA");

    dpd_buf4_init(&tauIJAB, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
    dpd_buf4_init(&tauijab, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauijab");
    dpd_buf4_init(&tauIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
    dpd_buf4_init(&tauIjbA, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjbA");

    dpd_buf4_init(&F_anti, CC_FINTS, 0, 10, 7, 10, 5, 1, "F <ia|bc>");
    dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");

    dpd_contract444(&tauIJAB, &F_anti, &ZIJMA, 0, 0, 1, 0);
    dpd_contract444(&tauijab, &F_anti, &Zijma, 0, 0, 1, 0);
    dpd_contract444(&tauIjAb, &F, &ZIjMa, 0, 0, 1, 0);
    dpd_contract444(&tauIjbA, &F, &ZIjmA, 0, 0, 1, 0);

    dpd_buf4_close(&tauIJAB); 
    dpd_buf4_close(&tauijab); 
    dpd_buf4_close(&tauIjAb);
    dpd_buf4_close(&tauIjbA);

    dpd_buf4_close(&F_anti); 
    dpd_buf4_close(&F);

    dpd_buf4_sort(&ZIJMA, CC_MISC, pqsr, 2, 11, "ZIJAM");
    dpd_buf4_sort(&Zijma, CC_MISC, pqsr, 2, 11, "Zijam");
    dpd_buf4_sort(&ZIjmA, CC_MISC, pqsr, 0, 11, "ZIjAm");

    dpd_buf4_close(&ZIJMA);  
    dpd_buf4_close(&Zijma);
    dpd_buf4_close(&ZIjMa);  
    dpd_buf4_close(&ZIjmA);
  }
  else if(params.ref == 2) { /*** UHF ***/

    dpd_buf4_init(&ZIJMA, CC_MISC, 0, 2, 20, 2, 20, 0, "ZIJMA");
    dpd_buf4_init(&Zijma, CC_MISC, 0, 12, 30, 12, 30, 0, "Zijma");
    dpd_buf4_init(&ZIjMa, CC_MISC, 0, 22, 24, 22, 24, 0, "ZIjMa");
    dpd_buf4_init(&ZIjAm, CC_MISC, 0, 22, 26, 22, 26, 0, "ZIjAm");

    dpd_buf4_init(&tauIJAB, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
    dpd_buf4_init(&tauijab, CC_TAMPS, 0, 12, 17, 12, 17, 0, "tauijab");
    dpd_buf4_init(&tauIjAb, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tauIjAb");

    dpd_buf4_init(&F, CC_FINTS, 0, 20, 7, 20, 5, 1, "F <IA|BC>");
    dpd_contract444(&tauIJAB, &F, &ZIJMA, 0, 0, 1, 0);
    dpd_buf4_close(&F);

    dpd_buf4_init(&F, CC_FINTS, 0, 30, 17, 30, 15, 1, "F <ia|bc>");
    dpd_contract444(&tauijab, &F, &Zijma, 0, 0, 1, 0);
    dpd_buf4_close(&F);

    dpd_buf4_init(&F, CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
    dpd_contract444(&tauIjAb, &F, &ZIjMa, 0, 0, 1, 0);
    dpd_buf4_close(&F);

    dpd_buf4_init(&F, CC_FINTS, 0, 28, 26, 28, 26, 0, "F <Ab|Ci>");
    dpd_contract444(&tauIjAb, &F, &ZIjAm, 0, 1, 1, 0);
    dpd_buf4_close(&F);

    dpd_buf4_close(&tauIJAB); 
    dpd_buf4_close(&tauijab); 
    dpd_buf4_close(&tauIjAb);

    dpd_buf4_sort(&ZIJMA, CC_MISC, pqsr, 2, 21, "ZIJAM");
    dpd_buf4_sort(&Zijma, CC_MISC, pqsr, 12, 31, "Zijam");

    dpd_buf4_close(&ZIJMA);  
    dpd_buf4_close(&Zijma);
    dpd_buf4_close(&ZIjMa);  
    dpd_buf4_close(&ZIjAm);

  }

#ifdef TIME_CCENERGY
  timer_off("Z");
#endif
}

}} // namespace psi::ccenergy
