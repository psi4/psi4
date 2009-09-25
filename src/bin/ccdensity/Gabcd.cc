/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

void Gabcd(void)
{
  dpdbuf4 G, L, T;
  int G_irr;
  G_irr = params.G_irr;

  if(params.ref == 0) { /** RHF **/
    dpd_buf4_init(&G, CC_GAMMA, G_irr, 5, 5, 5, 5, 0, "GAbCd");
    dpd_buf4_init(&L, CC_GLG, G_irr, 0, 5, 0, 5, 0, "LIjAb");
    dpd_buf4_init(&T, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
    dpd_contract444(&L, &T, &G, 1, 1, 1.0, 0.0);
    dpd_buf4_close(&T);
    dpd_buf4_close(&L);
    if (params.ground)
      dpd_buf4_symm(&G);
    dpd_buf4_close(&G);
  }
  else if(params.ref == 1) { /** RHF/ROHF **/

    dpd_buf4_init(&G, CC_GAMMA, G_irr, 7, 7, 7, 7, 0, "GABCD");
    dpd_buf4_init(&L, CC_GLG, G_irr, 2, 7, 2, 7, 0, "LIJAB");
    dpd_buf4_init(&T, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
    dpd_contract444(&L, &T, &G, 1, 1, 1.0, 0.0);
    dpd_buf4_close(&T);
    dpd_buf4_close(&L);
    if (params.ground)
      dpd_buf4_symm(&G);
    dpd_buf4_close(&G);

    dpd_buf4_init(&G, CC_GAMMA, G_irr, 7, 7, 7, 7, 0, "Gabcd");
    dpd_buf4_init(&L, CC_GLG, G_irr, 2, 7, 2, 7, 0, "Lijab");
    dpd_buf4_init(&T, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauijab");
    dpd_contract444(&L, &T, &G, 1, 1, 1.0, 0.0);
    dpd_buf4_close(&T);
    dpd_buf4_close(&L);
    if (params.ground)
      dpd_buf4_symm(&G);
    dpd_buf4_close(&G);

    dpd_buf4_init(&G, CC_GAMMA, G_irr, 5, 5, 5, 5, 0, "GAbCd");
    dpd_buf4_init(&L, CC_GLG, G_irr, 0, 5, 0, 5, 0, "LIjAb");
    dpd_buf4_init(&T, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
    dpd_contract444(&L, &T, &G, 1, 1, 1.0, 0.0);
    dpd_buf4_close(&T);
    dpd_buf4_close(&L);
    if (params.ground)
      dpd_buf4_symm(&G);
    dpd_buf4_close(&G);
  }
  else if(params.ref == 2) { /** UHF **/
    dpd_buf4_init(&G, CC_GAMMA, G_irr, 7, 7, 7, 7, 0, "GABCD");
    dpd_buf4_init(&L, CC_GLG, G_irr, 2, 7, 2, 7, 0, "LIJAB");
    dpd_buf4_init(&T, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
    dpd_contract444(&L, &T, &G, 1, 1, 1.0, 0.0);
    dpd_buf4_close(&T);
    dpd_buf4_close(&L);
    dpd_buf4_symm(&G);
    dpd_buf4_close(&G);

    dpd_buf4_init(&G, CC_GAMMA, G_irr, 17, 17, 17, 17, 0, "Gabcd");
    dpd_buf4_init(&L, CC_GLG, G_irr, 12, 17, 12, 17, 0, "Lijab");
    dpd_buf4_init(&T, CC_TAMPS, 0, 12, 17, 12, 17, 0, "tauijab");
    dpd_contract444(&L, &T, &G, 1, 1, 1.0, 0.0);
    dpd_buf4_close(&T);
    dpd_buf4_close(&L);
    dpd_buf4_symm(&G);
    dpd_buf4_close(&G);

    dpd_buf4_init(&G, CC_GAMMA, G_irr, 28, 28, 28, 28, 0, "GAbCd");
    dpd_buf4_init(&L, CC_GLG, G_irr, 22, 28, 22, 28, 0, "LIjAb");
    dpd_buf4_init(&T, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tauIjAb");
    dpd_contract444(&L, &T, &G, 1, 1, 1.0, 0.0);
    dpd_buf4_close(&T);
    dpd_buf4_close(&L);
    dpd_buf4_symm(&G);
    dpd_buf4_close(&G);
  }
}

}} // namespace psi::ccdensity
