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

void zero_onepdm(struct RHO_Params rho_params)
{
  dpdfile2 D;
  int G_irr;
  G_irr = params.G_irr;

  if ( (params.ref == 0) || (params.ref == 1) ) {
    dpd_file2_init(&D, CC_OEI, G_irr, 0, 0, rho_params.DIJ_lbl);
    dpd_file2_scm(&D, 0.0);
    dpd_file2_close(&D);
    dpd_file2_init(&D, CC_OEI, G_irr, 0, 0, rho_params.Dij_lbl);
    dpd_file2_scm(&D, 0.0);
    dpd_file2_close(&D);
  
    dpd_file2_init(&D, CC_OEI, G_irr, 1, 1, rho_params.DAB_lbl);
    dpd_file2_scm(&D, 0.0);
    dpd_file2_close(&D);
    dpd_file2_init(&D, CC_OEI, G_irr, 1, 1, rho_params.Dab_lbl);
    dpd_file2_scm(&D, 0.0);
    dpd_file2_close(&D);
  
    dpd_file2_init(&D, CC_OEI, G_irr, 0, 1, rho_params.DIA_lbl);
    dpd_file2_scm(&D, 0.0);
    dpd_file2_close(&D);
    dpd_file2_init(&D, CC_OEI, G_irr, 0, 1, rho_params.Dia_lbl);
    dpd_file2_scm(&D, 0.0);
    dpd_file2_close(&D);
  
    dpd_file2_init(&D, CC_OEI, G_irr, 0, 1, rho_params.DAI_lbl);
    dpd_file2_scm(&D, 0.0);
    dpd_file2_close(&D);
    dpd_file2_init(&D, CC_OEI, G_irr, 0, 1, rho_params.Dai_lbl);
    dpd_file2_scm(&D, 0.0);
    dpd_file2_close(&D);
  }
  else if (params.ref == 2) {
    dpd_file2_init(&D, CC_OEI, G_irr, 0, 0, rho_params.DIJ_lbl);
    dpd_file2_scm(&D, 0.0);
    dpd_file2_close(&D);
    dpd_file2_init(&D, CC_OEI, G_irr, 2, 2, rho_params.Dij_lbl);
    dpd_file2_scm(&D, 0.0);
    dpd_file2_close(&D);
  
    dpd_file2_init(&D, CC_OEI, G_irr, 1, 1, rho_params.DAB_lbl);
    dpd_file2_scm(&D, 0.0);
    dpd_file2_close(&D);
    dpd_file2_init(&D, CC_OEI, G_irr, 3, 3, rho_params.Dab_lbl);
    dpd_file2_scm(&D, 0.0);
    dpd_file2_close(&D);
  
    dpd_file2_init(&D, CC_OEI, G_irr, 0, 1, rho_params.DIA_lbl);
    dpd_file2_scm(&D, 0.0);
    dpd_file2_close(&D);
    dpd_file2_init(&D, CC_OEI, G_irr, 2, 3, rho_params.Dia_lbl);
    dpd_file2_scm(&D, 0.0);
    dpd_file2_close(&D);
  
    dpd_file2_init(&D, CC_OEI, G_irr, 0, 1, rho_params.DAI_lbl);
    dpd_file2_scm(&D, 0.0);
    dpd_file2_close(&D);
    dpd_file2_init(&D, CC_OEI, G_irr, 2, 3, rho_params.Dai_lbl);
    dpd_file2_scm(&D, 0.0);
    dpd_file2_close(&D);
  }
}

void zero_twopdm(void) 
{
  dpdbuf4 G;
  int G_irr;
  G_irr = params.G_irr;

  if ( (params.ref == 0) || (params.ref == 1) ) {
    dpd_buf4_init(&G, CC_GAMMA, G_irr, 2, 2, 2, 2, 0, "GIJKL");
    dpd_buf4_scm(&G, 0.0);
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, CC_GAMMA, G_irr, 2, 2, 2, 2, 0, "Gijkl");
    dpd_buf4_scm(&G, 0.0);
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, CC_GAMMA, G_irr, 0, 0, 0, 0, 0, "GIjKl");
    dpd_buf4_scm(&G, 0.0);
    dpd_buf4_close(&G);

    dpd_buf4_init(&G, CC_GAMMA, G_irr, 7, 7, 7, 7, 0, "GABCD");
    dpd_buf4_scm(&G, 0.0);
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, CC_GAMMA, G_irr, 7, 7, 7, 7, 0, "Gabcd");
    dpd_buf4_scm(&G, 0.0);
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, CC_GAMMA, G_irr, 5, 5, 5, 5, 0, "GAbCd");
    dpd_buf4_scm(&G, 0.0);
    dpd_buf4_close(&G);

    dpd_buf4_init(&G, CC_GAMMA, G_irr, 2, 10, 2, 10, 0, "GIJKA");
    dpd_buf4_scm(&G, 0.0);
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, CC_GAMMA, G_irr, 2, 10, 2, 10, 0, "Gijka");
    dpd_buf4_scm(&G, 0.0);
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, CC_GAMMA, G_irr, 0, 10, 0, 10, 0, "GIjKa");
    dpd_buf4_scm(&G, 0.0);
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, CC_GAMMA, G_irr, 0, 10, 0, 10, 0, "GiJkA");
    dpd_buf4_scm(&G, 0.0);
    dpd_buf4_close(&G);

    dpd_buf4_init(&G, CC_GAMMA, G_irr, 10, 10, 10, 10, 0, "GIBJA");
    dpd_buf4_scm(&G, 0.0);
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, CC_GAMMA, G_irr, 10, 10, 10, 10, 0, "Gibja");
    dpd_buf4_scm(&G, 0.0);
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, CC_GAMMA, G_irr, 10, 10, 10, 10, 0, "GIbJa");
    dpd_buf4_scm(&G, 0.0);
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, CC_GAMMA, G_irr, 10, 10, 10, 10, 0, "GiBjA");
    dpd_buf4_scm(&G, 0.0);
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, CC_GAMMA, G_irr, 10, 10, 10, 10, 0, "GIbjA");
    dpd_buf4_scm(&G, 0.0);
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, CC_GAMMA, G_irr, 10, 10, 10, 10, 0, "GiBJa");
    dpd_buf4_scm(&G, 0.0);
    dpd_buf4_close(&G);

    dpd_buf4_init(&G, CC_GAMMA, G_irr, 11, 7, 11, 7, 0, "GCIAB");
    dpd_buf4_scm(&G, 0.0);
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, CC_GAMMA, G_irr, 11, 7, 11, 7, 0, "Gciab");
    dpd_buf4_scm(&G, 0.0);
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, CC_GAMMA, G_irr, 11, 5, 11, 5, 0, "GCiAb");
    dpd_buf4_scm(&G, 0.0);
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, CC_GAMMA, G_irr, 11, 5, 11, 5, 0, "GcIaB");
    dpd_buf4_scm(&G, 0.0);
    dpd_buf4_close(&G);

    dpd_buf4_init(&G, CC_GAMMA, G_irr, 2, 7, 2, 7, 0, "GIJAB");
    dpd_buf4_scm(&G, 0.0);
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, CC_GAMMA, G_irr, 2, 7, 2, 7, 0, "Gijab");
    dpd_buf4_scm(&G, 0.0);
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, CC_GAMMA, G_irr, 0, 5, 0, 5, 0, "GIjAb");
    dpd_buf4_scm(&G, 0.0);
    dpd_buf4_close(&G);
  }
  else if (params.ref == 2) {
    dpd_buf4_init(&G, CC_GAMMA, G_irr, 2, 2, 2, 2, 0, "GIJKL");
    dpd_buf4_scm(&G, 0.0);
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, CC_GAMMA, G_irr, 12, 12, 12, 12, 0, "Gijkl");
    dpd_buf4_scm(&G, 0.0);
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, CC_GAMMA, G_irr, 22, 22, 22, 22, 0, "GIjKl");
    dpd_buf4_scm(&G, 0.0);
    dpd_buf4_close(&G);

    dpd_buf4_init(&G, CC_GAMMA, G_irr, 7, 7, 7, 7, 0, "GABCD");
    dpd_buf4_scm(&G, 0.0);
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, CC_GAMMA, G_irr, 17, 17, 17, 17, 0, "Gabcd");
    dpd_buf4_scm(&G, 0.0);
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, CC_GAMMA, G_irr, 28, 28, 28, 28, 0, "GAbCd");
    dpd_buf4_scm(&G, 0.0);
    dpd_buf4_close(&G);

    dpd_buf4_init(&G, CC_GAMMA, G_irr, 2, 20, 2, 20, 0, "GIJKA");
    dpd_buf4_scm(&G, 0.0);
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, CC_GAMMA, G_irr, 12, 30, 12, 30, 0, "Gijka");
    dpd_buf4_scm(&G, 0.0);
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, CC_GAMMA, G_irr, 22, 24, 22, 24, 0, "GIjKa");
    dpd_buf4_scm(&G, 0.0);
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, CC_GAMMA, G_irr, 23, 27, 23, 27, 0, "GiJkA");
    dpd_buf4_scm(&G, 0.0);
    dpd_buf4_close(&G);

    dpd_buf4_init(&G, CC_GAMMA, G_irr, 20, 20, 20, 20, 0, "GIBJA");
    dpd_buf4_scm(&G, 0.0);
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, CC_GAMMA, G_irr, 30, 30, 30, 30, 0, "Gibja");
    dpd_buf4_scm(&G, 0.0);
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, CC_GAMMA, G_irr, 24, 24, 24, 24, 0, "GIbJa");
    dpd_buf4_scm(&G, 0.0);
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, CC_GAMMA, G_irr, 27, 27, 27, 27, 0, "GiBjA");
    dpd_buf4_scm(&G, 0.0);
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, CC_GAMMA, G_irr, 24, 27, 24, 27, 0, "GIbjA");
    dpd_buf4_scm(&G, 0.0);
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, CC_GAMMA, G_irr, 27, 24, 27, 24, 0, "GiBJa");
    dpd_buf4_scm(&G, 0.0);
    dpd_buf4_close(&G);

    dpd_buf4_init(&G, CC_GAMMA, G_irr, 21, 7, 21, 7, 0, "GCIAB");
    dpd_buf4_scm(&G, 0.0);
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, CC_GAMMA, G_irr, 31, 17, 31, 17, 0, "Gciab");
    dpd_buf4_scm(&G, 0.0);
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, CC_GAMMA, G_irr, 26, 28, 26, 28, 0, "GCiAb");
    dpd_buf4_scm(&G, 0.0);
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, CC_GAMMA, G_irr, 25, 29, 25, 29, 0, "GcIaB");
    dpd_buf4_scm(&G, 0.0);
    dpd_buf4_close(&G);

    dpd_buf4_init(&G, CC_GAMMA, G_irr, 2, 7, 2, 7, 0, "GIJAB");
    dpd_buf4_scm(&G, 0.0);
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, CC_GAMMA, G_irr, 12, 17, 12, 17, 0, "Gijab");
    dpd_buf4_scm(&G, 0.0);
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, CC_GAMMA, G_irr, 24, 28, 24, 28, 0, "GIjAb");
    dpd_buf4_scm(&G, 0.0);
    dpd_buf4_close(&G);
  }
}

}} // namespace psi::ccdensity
