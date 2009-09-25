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

void WmnijT2(void)
{
  dpdbuf4 newtIJAB, newtijab, newtIjAb;
  dpdbuf4 WMNIJ, Wmnij, WMnIj;
  dpdbuf4 tauIJAB, tauijab, tauIjAb;

  if(params.ref == 0) { /** RHF **/
    dpd_buf4_init(&newtIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    dpd_buf4_init(&WMnIj, CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj");
    dpd_buf4_init(&tauIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
    dpd_contract444(&WMnIj, &tauIjAb, &newtIjAb, 1, 1, 1, 1);
    dpd_buf4_close(&tauIjAb);
    dpd_buf4_close(&WMnIj);
    dpd_buf4_close(&newtIjAb);
  }
  else if(params.ref == 1) { /** ROHF **/
    dpd_buf4_init(&newtIJAB, CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tIJAB");
    dpd_buf4_init(&WMNIJ, CC_HBAR, 0, 2, 2, 2, 2, 0, "WMNIJ");
    dpd_buf4_init(&tauIJAB, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
    dpd_contract444(&WMNIJ, &tauIJAB, &newtIJAB, 1, 1, 1, 1);
    dpd_buf4_close(&tauIJAB);
    dpd_buf4_close(&WMNIJ);
    dpd_buf4_close(&newtIJAB);

    dpd_buf4_init(&newtijab, CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tijab");
    dpd_buf4_init(&Wmnij, CC_HBAR, 0, 2, 2, 2, 2, 0, "Wmnij");
    dpd_buf4_init(&tauijab, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauijab");
    dpd_contract444(&Wmnij, &tauijab, &newtijab, 1, 1, 1, 1);
    dpd_buf4_close(&tauijab);
    dpd_buf4_close(&Wmnij);
    dpd_buf4_close(&newtijab);

    dpd_buf4_init(&newtIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    dpd_buf4_init(&WMnIj, CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj");
    dpd_buf4_init(&tauIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
    dpd_contract444(&WMnIj, &tauIjAb, &newtIjAb, 1, 1, 1, 1);
    dpd_buf4_close(&tauIjAb);
    dpd_buf4_close(&WMnIj);
    dpd_buf4_close(&newtIjAb);
  }
  else if(params.ref == 2) { /*** UHF ***/

    dpd_buf4_init(&newtIJAB, CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tIJAB");
    dpd_buf4_init(&WMNIJ, CC_HBAR, 0, 2, 2, 2, 2, 0, "WMNIJ");
    dpd_buf4_init(&tauIJAB, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
    dpd_contract444(&WMNIJ, &tauIJAB, &newtIJAB, 1, 1, 1, 1);
    dpd_buf4_close(&tauIJAB);
    dpd_buf4_close(&WMNIJ);
    dpd_buf4_close(&newtIJAB);

    dpd_buf4_init(&newtijab, CC_TAMPS, 0, 12, 17, 12, 17, 0, "New tijab");
    dpd_buf4_init(&Wmnij, CC_HBAR, 0, 12, 12, 12, 12, 0, "Wmnij");
    dpd_buf4_init(&tauijab, CC_TAMPS, 0, 12, 17, 12, 17, 0, "tauijab");
    dpd_contract444(&Wmnij, &tauijab, &newtijab, 1, 1, 1, 1);
    dpd_buf4_close(&tauijab);
    dpd_buf4_close(&Wmnij);
    dpd_buf4_close(&newtijab);

    dpd_buf4_init(&newtIjAb, CC_TAMPS, 0, 22, 28, 22, 28, 0, "New tIjAb");
    dpd_buf4_init(&WMnIj, CC_HBAR, 0, 22, 22, 22, 22, 0, "WMnIj");
    dpd_buf4_init(&tauIjAb, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tauIjAb");
    dpd_contract444(&WMnIj, &tauIjAb, &newtIjAb, 1, 1, 1, 1);
    dpd_buf4_close(&tauIjAb);
    dpd_buf4_close(&WMnIj);
    dpd_buf4_close(&newtIjAb);

  }
}
}} // namespace psi::ccenergy
