/*! \file
    \ingroup CCENERGY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <libdpd/dpd.h>
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccenergy {

void FaetT2(void)
{
  dpdfile2 FAEt, Faet;
  dpdbuf4 newtIJAB, newtijab, newtIjAb;
  dpdbuf4 tIJAB, tijab, tIjAb;
  dpdbuf4 t2, Z;

  if(params.ref == 0) { /** RHF **/
    dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_file2_init(&FAEt, CC_OEI, 0, 1, 1, "FAEt");
    dpd_buf4_init(&Z, CC_TMP0, 0, 0, 5, 0, 5, 0, "Zijab");
    dpd_contract424(&tIjAb, &FAEt, &Z, 3, 1, 0, 1, 0);
    dpd_file2_close(&FAEt);
    dpd_buf4_close(&tIjAb);
    dpd_buf4_init(&newtIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    dpd_buf4_axpy(&Z, &newtIjAb, 1);
    dpd_buf4_close(&newtIjAb);
    dpd_buf4_sort_axpy(&Z, CC_TAMPS, qpsr, 0, 5, "New tIjAb", 1);
    dpd_buf4_close(&Z);
  }
  else if(params.ref == 1) { /** ROHF **/
    dpd_buf4_init(&newtIJAB, CC_TAMPS, 0, 2, 5, 2, 7, 0, "New tIJAB");
    dpd_buf4_init(&newtijab, CC_TAMPS, 0, 2, 5, 2, 7, 0, "New tijab");
    dpd_buf4_init(&newtIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");

    dpd_buf4_init(&tIJAB, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    dpd_buf4_init(&tijab, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tijab");
    dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");

    dpd_file2_init(&FAEt, CC_OEI, 0, 1, 1, "FAEt");
    dpd_file2_init(&Faet, CC_OEI, 0, 1, 1, "Faet");

    dpd_buf4_init(&t2, CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,AB)");
    dpd_contract424(&tIJAB, &FAEt, &t2, 3, 1, 0, 1, 0);
    dpd_contract244(&FAEt, &tIJAB, &t2, 1, 2, 1, 1, 1);
    dpd_buf4_axpy(&t2, &newtIJAB, 1);
    dpd_buf4_close(&t2);

    dpd_buf4_init(&t2, CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,AB)");
    dpd_contract424(&tijab, &Faet, &t2, 3, 1, 0, 1, 0);
    dpd_contract244(&Faet, &tijab, &t2, 1, 2, 1, 1, 1);
    dpd_buf4_axpy(&t2, &newtijab, 1);
    dpd_buf4_close(&t2);

    dpd_contract424(&tIjAb, &Faet, &newtIjAb, 3, 1, 0, 1, 1);
    dpd_contract244(&FAEt, &tIjAb, &newtIjAb, 1, 2, 1, 1, 1);

    dpd_file2_close(&FAEt);  
    dpd_file2_close(&Faet);

    dpd_buf4_close(&tIJAB);
    dpd_buf4_close(&tijab);
    dpd_buf4_close(&tIjAb);
    dpd_buf4_close(&newtIJAB);
    dpd_buf4_close(&newtijab);
    dpd_buf4_close(&newtIjAb);
  }
  else if(params.ref == 2) { /*** UHF ***/

    dpd_buf4_init(&newtIJAB, CC_TAMPS, 0, 2, 5, 2, 7, 0, "New tIJAB");
    dpd_buf4_init(&newtijab, CC_TAMPS, 0, 12, 15, 12, 17, 0, "New tijab");
    dpd_buf4_init(&newtIjAb, CC_TAMPS, 0, 22, 28, 22, 28, 0, "New tIjAb");

    dpd_buf4_init(&tIJAB, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    dpd_buf4_init(&tijab, CC_TAMPS, 0, 12, 15, 12, 17, 0, "tijab");
    dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");

    dpd_file2_init(&FAEt, CC_OEI, 0, 1, 1, "FAEt");
    dpd_file2_init(&Faet, CC_OEI, 0, 3, 3, "Faet");

    dpd_buf4_init(&t2, CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,AB)");
    dpd_contract424(&tIJAB, &FAEt, &t2, 3, 1, 0, 1, 0);
    dpd_contract244(&FAEt, &tIJAB, &t2, 1, 2, 1, 1, 1);
    dpd_buf4_axpy(&t2, &newtIJAB, 1);
    dpd_buf4_close(&t2);

    dpd_buf4_init(&t2, CC_TMP0, 0, 12, 15, 12, 15, 0, "T (i>j,ab)");
    dpd_contract424(&tijab, &Faet, &t2, 3, 1, 0, 1, 0);
    dpd_contract244(&Faet, &tijab, &t2, 1, 2, 1, 1, 1);
    dpd_buf4_axpy(&t2, &newtijab, 1);
    dpd_buf4_close(&t2);

    dpd_contract424(&tIjAb, &Faet, &newtIjAb, 3, 1, 0, 1, 1);
    dpd_contract244(&FAEt, &tIjAb, &newtIjAb, 1, 2, 1, 1, 1);

    dpd_file2_close(&FAEt);  
    dpd_file2_close(&Faet);

    dpd_buf4_close(&tIJAB);
    dpd_buf4_close(&tijab);
    dpd_buf4_close(&tIjAb);

    dpd_buf4_close(&newtIJAB);
    dpd_buf4_close(&newtijab);
    dpd_buf4_close(&newtIjAb);

  }
}
}} // namespace psi::ccenergy
