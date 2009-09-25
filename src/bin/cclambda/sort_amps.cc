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

void sort_amps(int L_irr)
{
  dpdbuf4 L2;

  if(params.ref == 0) {
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    dpd_buf4_scmcopy(&L2, CC_LAMBDA, "2 LIjAb - LIjBa", 2);
    dpd_buf4_sort_axpy(&L2, CC_LAMBDA, pqsr, 0, 5, "2 LIjAb - LIjBa", -1);
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    dpd_buf4_sort(&L2, CC_LAMBDA, prqs, 10, 10, "LIAjb");
    dpd_buf4_sort(&L2, CC_LAMBDA, psqr, 10, 10, "LIbjA");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 10, 10, 10, 10, 0, "LIAjb");
    dpd_buf4_scmcopy(&L2, CC_LAMBDA, "2 LIAjb - LIbjA", 2);
    dpd_buf4_sort_axpy(&L2, CC_LAMBDA, psrq, 10, 10, "2 LIAjb - LIbjA", -1);
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    dpd_buf4_sort(&L2, CC_LAMBDA, qpsr, 0, 5, "LiJaB");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LiJaB");
    dpd_buf4_sort(&L2, CC_LAMBDA, prqs, 10, 10, "LiaJB");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 10, 10, 10, 10, 0, "LIAjb");
    dpd_buf4_sort(&L2, CC_LAMBDA, rqps, 10, 10, "LjAIb");
    dpd_buf4_close(&L2);
  }
  
  if(params.ref == 1) { /** ROHF **/
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    dpd_buf4_sort(&L2, CC_LAMBDA, prqs, 10, 10, "LIAjb");
    dpd_buf4_sort(&L2, CC_LAMBDA, psqr, 10, 10, "LIbjA");
    dpd_buf4_sort(&L2, CC_LAMBDA, qpsr, 0, 5, "LiJaB");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LiJaB");
    dpd_buf4_sort(&L2, CC_LAMBDA, prqs, 10, 10, "LiaJB");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 10, 10, 10, 10, 0, "LIAjb");
    dpd_buf4_sort(&L2, CC_LAMBDA, rqps, 10, 10, "LjAIb");
    dpd_buf4_close(&L2);
    
    /* Build L2IAJB List */
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 0, 5, 2, 7, 0, "LIJAB");
    dpd_buf4_sort(&L2, CC_LAMBDA, prqs, 10, 10, "LIAJB");
    dpd_buf4_close(&L2);
    /* Build L2iajb List */
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 0, 5, 2, 7, 0, "Lijab");
    dpd_buf4_sort(&L2, CC_LAMBDA, prqs, 10, 10, "Liajb");
    dpd_buf4_close(&L2);
  }
  else if(params.ref == 2) { /** UHF **/

    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "LIjAb");
    dpd_buf4_sort(&L2, CC_LAMBDA, qpsr, 23, 29, "LiJaB");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 0, 5, 2, 7, 0, "LIJAB");
    dpd_buf4_sort(&L2, CC_LAMBDA, prqs, 20, 20, "LIAJB");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 10, 15, 12, 17, 0, "Lijab");
    dpd_buf4_sort(&L2, CC_LAMBDA, prqs, 30, 30, "Liajb");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "LIjAb");
    dpd_buf4_sort(&L2, CC_LAMBDA, prqs, 20, 30, "LIAjb");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 23, 29, 23, 29, 0, "LiJaB");
    dpd_buf4_sort(&L2, CC_LAMBDA, prqs, 30, 20, "LiaJB");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 20, 30, 20, 30, 0, "LIAjb");
    dpd_buf4_sort(&L2, CC_LAMBDA, psrq, 24, 27, "LIbjA");
    dpd_buf4_sort(&L2, CC_LAMBDA, rqps, 27, 24, "LjAIb");
    dpd_buf4_close(&L2);
  }

}


}} // namespace psi::cclambda
