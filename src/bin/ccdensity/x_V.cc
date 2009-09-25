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

/* compute intermediates for excited state pdm that look like
   ground state intermediates, only use GL */

void V_build_x(void)
{
  int L_irr;
  dpdbuf4 V, L, T;
  L_irr = params.L_irr;

  if(params.ref == 0 || params.ref == 1) { /** RHF/ROHF **/

    dpd_buf4_init(&V, EOM_TMP, L_irr, 10, 10, 10, 10, 0, "VIAJB");
    dpd_buf4_init(&T, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
    dpd_buf4_init(&L, CC_GL, L_irr, 10, 10, 10, 10, 0, "LIAJB");
    dpd_contract444(&T, &L, &V, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&L);
    dpd_buf4_close(&T);
    dpd_buf4_init(&T, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    dpd_buf4_init(&L, CC_GL, L_irr, 10, 10, 10, 10, 0, "LIAjb");
    dpd_contract444(&T, &L, &V, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&L);
    dpd_buf4_close(&T);
    dpd_buf4_close(&V);

    dpd_buf4_init(&V, EOM_TMP, L_irr, 10, 10, 10, 10, 0, "Viajb");
    dpd_buf4_init(&T, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
    dpd_buf4_init(&L, CC_GL, L_irr, 10, 10, 10, 10, 0, "Liajb");
    dpd_contract444(&T, &L, &V, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&L);
    dpd_buf4_close(&T);
    dpd_buf4_init(&T, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    dpd_buf4_init(&L, CC_GL, L_irr, 10, 10, 10, 10, 0, "LIAjb");
    dpd_contract444(&T, &L, &V, 1, 1, 1.0, 1.0);
    dpd_buf4_close(&L);
    dpd_buf4_close(&T);
    dpd_buf4_close(&V);

    dpd_buf4_init(&V, EOM_TMP, L_irr, 10, 10, 10, 10, 0, "VIAjb");
    dpd_buf4_init(&T, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    dpd_buf4_init(&L, CC_GL, L_irr, 10, 10, 10, 10, 0, "Liajb");
    dpd_contract444(&T, &L, &V, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&L);
    dpd_buf4_close(&T);
    dpd_buf4_init(&T, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
    dpd_buf4_init(&L, CC_GL, L_irr, 10, 10, 10, 10, 0, "LIAjb");
    dpd_contract444(&T, &L, &V, 0, 1, 1.0, 1.0);
    dpd_buf4_close(&L);
    dpd_buf4_close(&T);
    dpd_buf4_close(&V);

    dpd_buf4_init(&V, EOM_TMP, L_irr, 10, 10, 10, 10, 0, "ViaJB");
    dpd_buf4_init(&T, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiaJB");
    dpd_buf4_init(&L, CC_GL, L_irr, 10, 10, 10, 10, 0, "LIAJB");
    dpd_contract444(&T, &L, &V, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&L);
    dpd_buf4_close(&T);
    dpd_buf4_init(&T, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
    dpd_buf4_init(&L, CC_GL, L_irr, 10, 10, 10, 10, 0, "LIAjb");
    dpd_contract444(&T, &L, &V, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&L);
    dpd_buf4_close(&T);
    dpd_buf4_close(&V);

    dpd_buf4_init(&V, EOM_TMP, L_irr, 10, 10, 10, 10, 0, "ViAjB");
    dpd_buf4_init(&L, CC_GL, L_irr, 10, 10, 10, 10, 0, "LIbjA");
    dpd_buf4_init(&T, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tjAIb");
    dpd_contract444(&T, &L, &V, 0, 1, 1.0, 0.0);
    dpd_buf4_close(&T);
    dpd_buf4_close(&L);
    dpd_buf4_close(&V);

    dpd_buf4_init(&V, EOM_TMP, L_irr, 10, 10, 10, 10, 0, "VIaJb");
    dpd_buf4_init(&L, CC_GL, L_irr, 10, 10, 10, 10, 0, "LjAIb");
    dpd_buf4_init(&T, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIbjA");
    dpd_contract444(&T, &L, &V, 0, 1, 1.0, 0.0);
    dpd_buf4_close(&T);
    dpd_buf4_close(&L);
    dpd_buf4_close(&V);
  }
  else if(params.ref == 2) { /** UHF **/

    dpd_buf4_init(&V, EOM_TMP, L_irr, 20, 20, 20, 20, 0, "VIAJB");
    dpd_buf4_init(&T, CC_TAMPS, 0, 20, 20, 20, 20, 0, "tIAJB");
    dpd_buf4_init(&L, CC_GL, L_irr, 20, 20, 20, 20, 0, "LIAJB");
    dpd_contract444(&T, &L, &V, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&L);
    dpd_buf4_close(&T);
    dpd_buf4_init(&T, CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
    dpd_buf4_init(&L, CC_GL, L_irr, 20, 30, 20, 30, 0, "LIAjb");
    dpd_contract444(&T, &L, &V, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&L);
    dpd_buf4_close(&T);
    dpd_buf4_close(&V);

    dpd_buf4_init(&V, EOM_TMP, L_irr, 30, 30, 30, 30, 0, "Viajb");
    dpd_buf4_init(&T, CC_TAMPS, 0, 30, 30, 30, 30, 0, "tiajb");
    dpd_buf4_init(&L, CC_GL, L_irr, 30, 30, 30, 30, 0, "Liajb");
    dpd_contract444(&T, &L, &V, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&L);
    dpd_buf4_close(&T);
    dpd_buf4_init(&T, CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
    dpd_buf4_init(&L, CC_GL, L_irr, 20, 30, 20, 30, 0, "LIAjb");
    dpd_contract444(&T, &L, &V, 1, 1, 1.0, 1.0);
    dpd_buf4_close(&L);
    dpd_buf4_close(&T);
    dpd_buf4_close(&V);

    dpd_buf4_init(&V, EOM_TMP, L_irr, 20, 30, 20, 30, 0, "VIAjb");
    dpd_buf4_init(&T, CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
    dpd_buf4_init(&L, CC_GL, L_irr, 30, 30, 30, 30, 0, "Liajb");
    dpd_contract444(&T, &L, &V, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&L);
    dpd_buf4_close(&T);
    dpd_buf4_init(&T, CC_TAMPS, 0, 20, 20, 20, 20, 0, "tIAJB");
    dpd_buf4_init(&L, CC_GL, L_irr, 20, 30, 20, 30, 0, "LIAjb");
    dpd_contract444(&T, &L, &V, 0, 1, 1.0, 1.0);
    dpd_buf4_close(&L);
    dpd_buf4_close(&T);
    dpd_buf4_close(&V);

    dpd_buf4_init(&V, EOM_TMP, L_irr, 30, 20, 30, 20, 0, "ViaJB");
    dpd_buf4_init(&T, CC_TAMPS, 0, 30, 20, 30, 20, 0, "tiaJB");
    dpd_buf4_init(&L, CC_GL, L_irr, 20, 20, 20, 20, 0, "LIAJB");
    dpd_contract444(&T, &L, &V, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&L);
    dpd_buf4_close(&T);
    dpd_buf4_init(&T, CC_TAMPS, 0, 30, 30, 30, 30, 0, "tiajb");
    dpd_buf4_init(&L, CC_GL, L_irr, 20, 30, 20, 30, 0, "LIAjb");
    dpd_contract444(&T, &L, &V, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&L);
    dpd_buf4_close(&T);
    dpd_buf4_close(&V);

    dpd_buf4_init(&V, EOM_TMP, L_irr, 27, 27, 27, 27, 0, "ViAjB");
    dpd_buf4_init(&L, CC_GL, L_irr, 24, 27, 24, 27, 0, "LIbjA");
    dpd_buf4_init(&T, CC_TAMPS, 0, 27, 24, 27, 24, 0, "tjAIb");
    dpd_contract444(&T, &L, &V, 0, 1, 1.0, 0.0);
    dpd_buf4_close(&T);
    dpd_buf4_close(&L);
    dpd_buf4_close(&V);

    dpd_buf4_init(&V, EOM_TMP, L_irr, 24, 24, 24, 24, 0, "VIaJb");
    dpd_buf4_init(&L, CC_GL, L_irr, 27, 24, 27, 24, 0, "LjAIb");
    dpd_buf4_init(&T, CC_TAMPS, 0, 24, 27, 24, 27, 0, "tIbjA");
    dpd_contract444(&T, &L, &V, 0, 1, 1.0, 0.0);
    dpd_buf4_close(&T);
    dpd_buf4_close(&L);
    dpd_buf4_close(&V);
  }
}

}} // namespace psi::ccdensity
