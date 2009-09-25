/*! \file
    \ingroup CCLAMBDA
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cclambda {

void WijmnL2(int L_irr)
{
  dpdbuf4 Lijab, LIJAB, LIjAb;
  dpdbuf4 newLijab, newLIJAB, newLIjAb;
  dpdbuf4 WMNIJ, Wmnij, WMnIj;

  /* RHS += Lmnab*Wijmn */
  if(params.ref == 0) { /** RHF **/
    dpd_buf4_init(&LIjAb, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    dpd_buf4_init(&newLIjAb, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
    dpd_buf4_init(&WMnIj, CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj");
    dpd_contract444(&WMnIj, &LIjAb, &newLIjAb, 0, 1, 1.0, 1.0);
    dpd_buf4_close(&WMnIj);
    dpd_buf4_close(&LIjAb);
    dpd_buf4_close(&newLIjAb);
  }
  else if(params.ref == 1) { /** ROHF **/

    dpd_buf4_init(&LIJAB, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
    dpd_buf4_init(&newLIJAB, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB");
    dpd_buf4_init(&WMNIJ, CC_HBAR, 0, 2, 2, 2, 2, 0, "WMNIJ");
    dpd_contract444(&WMNIJ, &LIJAB, &newLIJAB, 0, 1, 1.0, 1.0);
    dpd_buf4_close(&WMNIJ);
    dpd_buf4_close(&LIJAB);
    dpd_buf4_close(&newLIJAB);

    dpd_buf4_init(&Lijab, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "Lijab");
    dpd_buf4_init(&newLijab, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New Lijab");
    dpd_buf4_init(&Wmnij, CC_HBAR, 0, 2, 2, 2, 2, 0, "Wmnij");
    dpd_contract444(&Wmnij, &Lijab, &newLijab, 0, 1, 1.0, 1.0);
    dpd_buf4_close(&Wmnij);
    dpd_buf4_close(&Lijab);
    dpd_buf4_close(&newLijab);

    dpd_buf4_init(&LIjAb, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    dpd_buf4_init(&newLIjAb, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
    dpd_buf4_init(&WMnIj, CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj");
    dpd_contract444(&WMnIj, &LIjAb, &newLIjAb, 0, 1, 1.0, 1.0);
    dpd_buf4_close(&WMnIj);
    dpd_buf4_close(&LIjAb);
    dpd_buf4_close(&newLIjAb);
  }
  else if(params.ref == 2) { /** UHF **/

    dpd_buf4_init(&LIJAB, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
    dpd_buf4_init(&newLIJAB, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB");
    dpd_buf4_init(&WMNIJ, CC_HBAR, 0, 2, 2, 2, 2, 0, "WMNIJ");
    dpd_contract444(&WMNIJ, &LIJAB, &newLIJAB, 0, 1, 1, 1);
    dpd_buf4_close(&WMNIJ);
    dpd_buf4_close(&LIJAB);
    dpd_buf4_close(&newLIJAB);

    dpd_buf4_init(&Lijab, CC_LAMBDA, L_irr, 12, 17, 12, 17, 0, "Lijab");
    dpd_buf4_init(&newLijab, CC_LAMBDA, L_irr, 12, 17, 12, 17, 0, "New Lijab");
    dpd_buf4_init(&Wmnij, CC_HBAR, 0, 12, 12, 12, 12, 0, "Wmnij");
    dpd_contract444(&Wmnij, &Lijab, &newLijab, 0, 1, 1, 1);
    dpd_buf4_close(&Wmnij);
    dpd_buf4_close(&Lijab);
    dpd_buf4_close(&newLijab);

    dpd_buf4_init(&LIjAb, CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "LIjAb");
    dpd_buf4_init(&newLIjAb, CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "New LIjAb");
    dpd_buf4_init(&WMnIj, CC_HBAR, 0, 22, 22, 22, 22, 0, "WMnIj");
    dpd_contract444(&WMnIj, &LIjAb, &newLIjAb, 0, 1, 1, 1);
    dpd_buf4_close(&WMnIj);
    dpd_buf4_close(&LIjAb);
    dpd_buf4_close(&newLIjAb);
  }
}


}} // namespace psi::cclambda
