/*! \file
    \ingroup CCEOM
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cceom {

/*
 copied from cclamba to make consistent copies of R and L
 changes should be made to sort_amps in both programs as
 cceom_density is written
*/

void sort_amps(void)
{
  dpdbuf4 R2;
  int R_irr;

  /* calculate irrep of R, the irrep for root of interest. */
  R_irr = eom_params.prop_sym^moinfo.sym;

  if(params.ref == 0 || params.ref == 1) { /* RHF/ROHF */
    /* Build R2iJaB list */
    dpd_buf4_init(&R2, CC_RAMPS, R_irr, 0, 5, 0, 5, 0, "RIjAb");
    dpd_buf4_sort(&R2, CC_RAMPS, qpsr, 0, 5, "RiJaB");
    dpd_buf4_close(&R2);

    /* Build R2IAJB list */
    dpd_buf4_init(&R2, CC_RAMPS, R_irr, 0, 5, 2, 7, 0, "RIJAB");
    dpd_buf4_sort(&R2, CC_RAMPS, prqs, 10, 10, "RIAJB");
    dpd_buf4_close(&R2);

    /* Build R2iajb list */
    dpd_buf4_init(&R2, CC_RAMPS, R_irr, 0, 5, 2, 7, 0, "Rijab");
    dpd_buf4_sort(&R2, CC_RAMPS, prqs, 10, 10, "Riajb");
    dpd_buf4_close(&R2);

    /* Build R2IAjb list */
    dpd_buf4_init(&R2, CC_RAMPS, R_irr, 0, 5, 0, 5, 0, "RIjAb");
    dpd_buf4_sort(&R2, CC_RAMPS, prqs, 10, 10, "RIAjb");
    dpd_buf4_close(&R2);

    /* Build R2iaJB list */
    dpd_buf4_init(&R2, CC_RAMPS, R_irr, 0, 5, 0, 5, 0, "RiJaB");
    dpd_buf4_sort(&R2, CC_RAMPS, prqs, 10, 10, "RiaJB");
    dpd_buf4_close(&R2);

    /* Build R2IbjA and R2 jAIb list */
    dpd_buf4_init(&R2, CC_RAMPS, R_irr, 10, 10, 10, 10, 0, "RIAjb");
    dpd_buf4_sort(&R2, CC_RAMPS, psrq, 10, 10, "RIbjA");
    dpd_buf4_sort(&R2, CC_RAMPS, rqps, 10, 10, "RjAIb");
    dpd_buf4_close(&R2);
  }
  else if(params.ref == 2) { /* UHF */

    dpd_buf4_init(&R2, CC_RAMPS, R_irr, 22, 28, 22, 28, 0, "RIjAb");
    dpd_buf4_sort(&R2, CC_RAMPS, qpsr, 23, 29, "RiJaB");
    dpd_buf4_close(&R2);

    dpd_buf4_init(&R2, CC_RAMPS, R_irr, 0, 5, 2, 7, 0, "RIJAB");
    dpd_buf4_sort(&R2, CC_RAMPS, prqs, 20, 20, "RIAJB");
    dpd_buf4_close(&R2);

    dpd_buf4_init(&R2, CC_RAMPS, R_irr, 10, 15, 12, 17, 0, "Rijab");
    dpd_buf4_sort(&R2, CC_RAMPS, prqs, 30, 30, "Riajb");
    dpd_buf4_close(&R2);

    dpd_buf4_init(&R2, CC_RAMPS, R_irr, 22, 28, 22, 28, 0, "RIjAb");
    dpd_buf4_sort(&R2, CC_RAMPS, prqs, 20, 30, "RIAjb");
    dpd_buf4_close(&R2);

    dpd_buf4_init(&R2, CC_RAMPS, R_irr, 23, 29, 23, 29, 0, "RiJaB");
    dpd_buf4_sort(&R2, CC_RAMPS, prqs, 30, 20, "RiaJB");
    dpd_buf4_close(&R2);

    dpd_buf4_init(&R2, CC_RAMPS, R_irr, 20, 30, 20, 30, 0, "RIAjb");
    dpd_buf4_sort(&R2, CC_RAMPS, psrq, 24, 27, "RIbjA");
    dpd_buf4_sort(&R2, CC_RAMPS, rqps, 27, 24, "RjAIb");
    dpd_buf4_close(&R2);
  }

}


}} // namespace psi::cceom
