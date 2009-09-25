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

void Lsave(int L_irr)
{
  dpdfile2 L1;
  dpdbuf4 L2;

  if(params.ref == 0 || params.ref == 1) { /** ROHF **/

    dpd_file2_init(&L1, CC_LAMBDA, L_irr, 0, 1, "New LIA");
    dpd_file2_copy(&L1, CC_LAMBDA, "LIA");
    dpd_file2_close(&L1);

    dpd_file2_init(&L1, CC_LAMBDA, L_irr, 0, 1, "New Lia");
    dpd_file2_copy(&L1, CC_LAMBDA, "Lia");
    dpd_file2_close(&L1);

    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB");
    dpd_buf4_copy(&L2, CC_LAMBDA, "LIJAB");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New Lijab");
    dpd_buf4_copy(&L2, CC_LAMBDA, "Lijab");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
    dpd_buf4_copy(&L2, CC_LAMBDA, "LIjAb");
    dpd_buf4_close(&L2);
  }
  else if(params.ref == 2) { /** UHF **/
    dpd_file2_init(&L1, CC_LAMBDA, L_irr, 0, 1, "New LIA");
    dpd_file2_copy(&L1, CC_LAMBDA, "LIA");
    dpd_file2_close(&L1);

    dpd_file2_init(&L1, CC_LAMBDA, L_irr, 2, 3, "New Lia");
    dpd_file2_copy(&L1, CC_LAMBDA, "Lia");
    dpd_file2_close(&L1);

    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB");
    dpd_buf4_copy(&L2, CC_LAMBDA, "LIJAB");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 12, 17, 12, 17, 0, "New Lijab");
    dpd_buf4_copy(&L2, CC_LAMBDA, "Lijab");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "New LIjAb");
    dpd_buf4_copy(&L2, CC_LAMBDA, "LIjAb");
    dpd_buf4_close(&L2);

  }
}


}} // namespace psi::cclambda
