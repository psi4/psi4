/*! \file
    \ingroup CCENERGY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccenergy {

void Wmnij_build(void)
{
  dpdbuf4 A_anti, A;
  dpdbuf4 WMNIJ, Wmnij, WMnIj, W;
  dpdfile2 tIA, tia;
  dpdbuf4 Eijka, Eijka_anti, Eaijk, Eaijk_anti;
  dpdbuf4 D_anti, D, tauIJAB, tauijab, tauIjAb;

  if(params.ref == 0) { /** RHF **/
    dpd_buf4_init(&A, CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
    dpd_buf4_copy(&A, CC_HBAR, "WMnIj");
    dpd_buf4_close(&A);
  }
  else if(params.ref == 1) { /** ROHF **/
    dpd_buf4_init(&A_anti, CC_AINTS, 0, 2, 2, 0, 0, 1, "A <ij|kl>");
    dpd_buf4_copy(&A_anti, CC_HBAR, "WMNIJ");
    dpd_buf4_copy(&A_anti, CC_HBAR, "Wmnij");
    dpd_buf4_close(&A_anti);

    dpd_buf4_init(&A, CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
    dpd_buf4_copy(&A, CC_HBAR, "WMnIj");
    dpd_buf4_close(&A);
  }
  else if(params.ref == 2) { /*** UHF ***/
    dpd_buf4_init(&A, CC_AINTS, 0, 2, 2, 0, 0, 1, "A <IJ|KL>");
    dpd_buf4_copy(&A, CC_HBAR, "WMNIJ");
    dpd_buf4_close(&A);

    dpd_buf4_init(&A, CC_AINTS, 0, 12, 12, 10, 10, 1, "A <ij|kl>");
    dpd_buf4_copy(&A, CC_HBAR, "Wmnij");
    dpd_buf4_close(&A);

    dpd_buf4_init(&A, CC_AINTS, 0, 22, 22, 22, 22, 0, "A <Ij|Kl>");
    dpd_buf4_copy(&A, CC_HBAR, "WMnIj");
    dpd_buf4_close(&A);
  }

  if(params.ref == 0) { /** RHF **/
    dpd_buf4_init(&WMnIj, CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj");
    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");

    dpd_buf4_init(&Eaijk, CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    dpd_contract244(&tIA, &Eaijk, &WMnIj, 1, 0, 1, 1, 1);
    dpd_buf4_close(&Eaijk);

    dpd_buf4_init(&Eijka, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    dpd_contract424(&Eijka, &tIA, &WMnIj, 3, 1, 0, 1, 1);
    dpd_buf4_close(&Eijka);

    dpd_file2_close(&tIA);
    dpd_buf4_close(&WMnIj);
  }
  else if(params.ref == 1) { /** ROHF **/  
    dpd_buf4_init(&WMNIJ, CC_HBAR, 0, 2, 0, 2, 2, 0, "WMNIJ");
    dpd_buf4_init(&Wmnij, CC_HBAR, 0, 2, 0, 2, 2, 0, "Wmnij");
    dpd_buf4_init(&WMnIj, CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj");

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");

    dpd_buf4_init(&Eijka_anti, CC_EINTS, 0, 2, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
    dpd_buf4_init(&Eijka, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    dpd_buf4_init(&Eaijk_anti, CC_EINTS, 0, 11, 2, 11, 0, 1, "E <ai|jk>");
    dpd_buf4_init(&Eaijk, CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");

    dpd_buf4_init(&W, CC_TMP0, 0, 2, 0, 2, 0, 0, "W (MN,IJ)");
    dpd_contract424(&Eijka_anti, &tIA, &W, 3, 1, 0, 1, 0);
    dpd_contract244(&tIA, &Eaijk_anti, &W, 1, 0, 1, 1, 1);
    dpd_buf4_axpy(&W, &WMNIJ, 1);
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP0, 0, 2, 0, 2, 0, 0, "W (MN,IJ)");
    dpd_contract424(&Eijka_anti, &tia, &W, 3, 1, 0, 1, 0);
    dpd_contract244(&tia, &Eaijk_anti, &W, 1, 0, 1, 1, 1);
    dpd_buf4_axpy(&W, &Wmnij, 1);
    dpd_buf4_close(&W);

    dpd_contract424(&Eijka, &tia, &WMnIj, 3, 1, 0, 1, 1);
    dpd_contract244(&tIA, &Eaijk, &WMnIj, 1, 0, 1, 1, 1);

    dpd_buf4_close(&Eijka_anti);
    dpd_buf4_close(&Eijka);
    dpd_buf4_close(&Eaijk_anti);
    dpd_buf4_close(&Eaijk);

    dpd_file2_close(&tIA);
    dpd_file2_close(&tia);

    dpd_buf4_close(&WMNIJ);
    dpd_buf4_close(&Wmnij);
    dpd_buf4_close(&WMnIj);
  }
  else if(params.ref == 2) { /*** UHF ***/

    dpd_buf4_init(&WMNIJ, CC_HBAR, 0, 2, 0, 2, 2, 0, "WMNIJ");
    dpd_buf4_init(&Wmnij, CC_HBAR, 0, 12, 10, 12, 12, 0, "Wmnij");
    dpd_buf4_init(&WMnIj, CC_HBAR, 0, 22, 22, 22, 22, 0, "WMnIj");

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 2, 3, "tia");

    dpd_buf4_init(&Eijka, CC_EINTS, 0, 2, 20, 2, 20, 0, "E <IJ||KA> (I>J,KA)");
    dpd_buf4_init(&Eaijk, CC_EINTS, 0, 21, 2, 21, 0, 1, "E <AI|JK>");
    dpd_buf4_init(&W, CC_TMP0, 0, 2, 0, 2, 0, 0, "W (MN,IJ)");
    dpd_contract424(&Eijka, &tIA, &W, 3, 1, 0, 1, 0);
    dpd_contract244(&tIA, &Eaijk, &W, 1, 0, 1, 1, 1);
    dpd_buf4_axpy(&W, &WMNIJ, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Eijka);
    dpd_buf4_close(&Eaijk);

    dpd_buf4_init(&Eijka, CC_EINTS, 0, 12, 30, 12, 30, 0, "E <ij||ka> (i>j,ka)");
    dpd_buf4_init(&Eaijk, CC_EINTS, 0, 31, 12, 31, 10, 1, "E <ai|jk>");
    dpd_buf4_init(&W, CC_TMP0, 0, 12, 10, 12, 10, 0, "W (mn,ij)");
    dpd_contract424(&Eijka, &tia, &W, 3, 1, 0, 1, 0);
    dpd_contract244(&tia, &Eaijk, &W, 1, 0, 1, 1, 1);
    dpd_buf4_axpy(&W, &Wmnij, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Eijka);
    dpd_buf4_close(&Eaijk);

    dpd_buf4_init(&Eijka, CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
    dpd_buf4_init(&Eaijk, CC_EINTS, 0, 26, 22, 26, 22, 0, "E <Ai|Jk>");
    dpd_contract424(&Eijka, &tia, &WMnIj, 3, 1, 0, 1, 1);
    dpd_contract244(&tIA, &Eaijk, &WMnIj, 1, 0, 1, 1, 1);
    dpd_buf4_close(&Eijka);
    dpd_buf4_close(&Eaijk);

    dpd_file2_close(&tIA);
    dpd_file2_close(&tia);

    dpd_buf4_close(&WMNIJ);
    dpd_buf4_close(&Wmnij);
    dpd_buf4_close(&WMnIj);

  }

  if(params.ref == 0) { /** RHF **/
    dpd_buf4_init(&WMnIj, CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_buf4_init(&tauIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
    dpd_contract444(&D, &tauIjAb, &WMnIj, 0, 0, 1, 1);
    dpd_buf4_close(&tauIjAb);
    dpd_buf4_close(&D);
    dpd_buf4_close(&WMnIj);
  }
  else if(params.ref == 1) { /** ROHF **/
    dpd_buf4_init(&WMNIJ, CC_HBAR, 0, 2, 2, 2, 2, 0, "WMNIJ");
    dpd_buf4_init(&Wmnij, CC_HBAR, 0, 2, 2, 2, 2, 0, "Wmnij");
    dpd_buf4_init(&WMnIj, CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj");

    dpd_buf4_init(&D_anti, CC_DINTS, 0, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");

    dpd_buf4_init(&tauIJAB, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
    dpd_buf4_init(&tauijab, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauijab");
    dpd_buf4_init(&tauIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");

    dpd_contract444(&D_anti, &tauIJAB, &WMNIJ, 0, 0, 1, 1);
    dpd_contract444(&D_anti, &tauijab, &Wmnij, 0, 0, 1, 1);
    dpd_contract444(&D, &tauIjAb, &WMnIj, 0, 0, 1, 1);

    dpd_buf4_close(&tauIJAB);
    dpd_buf4_close(&tauijab);
    dpd_buf4_close(&tauIjAb);

    dpd_buf4_close(&D_anti);
    dpd_buf4_close(&D);

    dpd_buf4_close(&WMNIJ);
    dpd_buf4_close(&Wmnij);
    dpd_buf4_close(&WMnIj);
  }
  else if(params.ref == 2) { /*** UHF ***/
    dpd_buf4_init(&WMNIJ, CC_HBAR, 0, 2, 2, 2, 2, 0, "WMNIJ");
    dpd_buf4_init(&Wmnij, CC_HBAR, 0, 12, 12, 12, 12, 0, "Wmnij");
    dpd_buf4_init(&WMnIj, CC_HBAR, 0, 22, 22, 22, 22, 0, "WMnIj");

    dpd_buf4_init(&D, CC_DINTS, 0, 2, 7, 2, 7, 0, "D <IJ||AB> (I>J,A>B)");
    dpd_buf4_init(&tauIJAB, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
    dpd_contract444(&D, &tauIJAB, &WMNIJ, 0, 0, 1, 1);
    dpd_buf4_close(&tauIJAB);
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, CC_DINTS, 0, 12, 17, 12, 17, 0, "D <ij||ab> (i>j,a>b)");
    dpd_buf4_init(&tauijab, CC_TAMPS, 0, 12, 17, 12, 17, 0, "tauijab");
    dpd_contract444(&D, &tauijab, &Wmnij, 0, 0, 1, 1);
    dpd_buf4_close(&tauijab);
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    dpd_buf4_init(&tauIjAb, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tauIjAb");
    dpd_contract444(&D, &tauIjAb, &WMnIj, 0, 0, 1, 1);
    dpd_buf4_close(&tauIjAb);
    dpd_buf4_close(&D);

    dpd_buf4_close(&WMNIJ);
    dpd_buf4_close(&Wmnij);
    dpd_buf4_close(&WMnIj);
  }

}
}} // namespace psi::ccenergy
