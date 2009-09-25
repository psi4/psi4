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

void cc2_WmbijT2(void) {

  dpdfile2 t1, tia, tIA;
  dpdbuf4 Z, W;
  dpdbuf4 t2, t2a, t2b, tIJAB, tijab, tIjAb;

  if(params.ref == 0) { /** RHF **/

    dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");

    dpd_buf4_init(&Z, CC_TMP0, 0, 5, 0, 5, 0, 0, "CC2 ZAbIj");
    dpd_buf4_init(&W, CC2_HET1, 0, 10, 0, 10, 0, 0, "CC2 WMbIj");
    dpd_contract244(&t1, &W, &Z, 0, 0, 0, -1, 0);
    dpd_buf4_close(&W);

    dpd_buf4_sort_axpy(&Z, CC_TAMPS, rspq, 0, 5, "New tIjAb", 1);
    dpd_buf4_sort_axpy(&Z, CC_TAMPS, srqp, 0, 5, "New tIjAb", 1);
    dpd_buf4_close(&Z);

    dpd_file2_close(&t1);
  }
  else if(params.ref == 1) { /** ROHF **/  
    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");

    /*** AA ***/
    dpd_buf4_init(&W, CC2_HET1, 0, 10, 2, 10, 2, 0, "CC2 WMBIJ (MB,I>J)");
    dpd_buf4_init(&t2, CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,AB)");
    dpd_contract244(&tIA, &W, &t2, 0, 0, 1, -1, 0);
    dpd_buf4_sort(&t2, CC_TMP0, pqsr, 2, 5, "T (I>J,BA)");
    dpd_buf4_close(&t2);
    dpd_buf4_init(&t2a, CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,AB)");
    dpd_buf4_init(&t2b, CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,BA)");
    dpd_buf4_axpy(&t2b, &t2a, -1);
    dpd_buf4_close(&t2b);
    dpd_buf4_init(&tIJAB, CC_TAMPS, 0, 2, 5, 2, 7, 0, "New tIJAB");
    dpd_buf4_axpy(&t2a, &tIJAB, 1);
    dpd_buf4_close(&tIJAB);
    dpd_buf4_close(&t2a);
    dpd_buf4_close(&W);

    /*** BB ***/
    dpd_buf4_init(&W, CC2_HET1, 0, 10, 2, 10, 2, 0, "CC2 Wmbij (mb,i>j)");
    dpd_buf4_init(&t2, CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,AB)");
    dpd_contract244(&tia, &W, &t2, 0, 0, 1, -1, 0);
    dpd_buf4_sort(&t2, CC_TMP0, pqsr, 2, 5, "T (I>J,BA)");
    dpd_buf4_close(&t2);
    dpd_buf4_init(&t2a, CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,AB)");
    dpd_buf4_init(&t2b, CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,BA)");
    dpd_buf4_axpy(&t2b, &t2a, -1);
    dpd_buf4_close(&t2b);
    dpd_buf4_init(&tijab, CC_TAMPS, 0, 2, 5, 2, 7, 0, "New tijab");
    dpd_buf4_axpy(&t2a, &tijab, 1);
    dpd_buf4_close(&tijab);
    dpd_buf4_close(&t2a);
    dpd_buf4_close(&W);

    /*** AB ***/
    dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    dpd_buf4_init(&W, CC2_HET1, 0, 10, 0, 10, 0, 0, "CC2 WMbIj");
    dpd_contract244(&tIA, &W, &tIjAb, 0, 0, 1, -1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&tIjAb);
    dpd_buf4_init(&Z, CC_TMP0, 0, 0, 5, 0, 5, 0, "CC2 ZjIbA");
    dpd_buf4_init(&W, CC2_HET1, 0, 10, 0, 10, 0, 0, "CC2 WmBiJ (mB,iJ)");
    dpd_contract244(&tia, &W, &Z, 0, 0, 1, -1, 0);
    dpd_buf4_close(&W);
    dpd_buf4_sort_axpy(&Z, CC_TAMPS, qpsr, 0, 5, "New tIjAb", 1);
    dpd_buf4_close(&Z);

    dpd_file2_close(&tIA); 
    dpd_file2_close(&tia);
  }
  else if(params.ref == 2) { /*** UHF ***/

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 2, 3, "tia");

    /*** AA ***/
    dpd_buf4_init(&W, CC2_HET1, 0, 20, 2, 20, 2, 0, "CC2 WMBIJ (MB,I>J)");
    dpd_buf4_init(&t2, CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,AB)");
    dpd_contract244(&tIA, &W, &t2, 0, 0, 1, -1, 0);
    dpd_buf4_sort(&t2, CC_TMP0, pqsr, 2, 5, "T (I>J,BA)");
    dpd_buf4_close(&t2);
    dpd_buf4_init(&t2a, CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,AB)");
    dpd_buf4_init(&t2b, CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,BA)");
    dpd_buf4_axpy(&t2b, &t2a, -1);
    dpd_buf4_close(&t2b);
    dpd_buf4_init(&tIJAB, CC_TAMPS, 0, 2, 5, 2, 7, 0, "New tIJAB");
    dpd_buf4_axpy(&t2a, &tIJAB, 1);
    dpd_buf4_close(&t2a);
    dpd_buf4_close(&W);
    dpd_buf4_close(&tIJAB);

    /*** BB ***/
    dpd_buf4_init(&W, CC2_HET1, 0, 30, 12, 30, 12, 0, "CC2 Wmbij (mb,i>j)");
    dpd_buf4_init(&t2, CC_TMP0, 0, 12, 15, 12, 15, 0, "T (i>j,ab)");
    dpd_contract244(&tia, &W, &t2, 0, 0, 1, -1, 0);
    dpd_buf4_sort(&t2, CC_TMP0, pqsr, 12, 15, "T (i>j,ba)");
    dpd_buf4_close(&t2);
    dpd_buf4_init(&t2a, CC_TMP0, 0, 12, 15, 12, 15, 0, "T (i>j,ab)");
    dpd_buf4_init(&t2b, CC_TMP0, 0, 12, 15, 12, 15, 0, "T (i>j,ba)");
    dpd_buf4_axpy(&t2b, &t2a, -1);
    dpd_buf4_close(&t2b);
    dpd_buf4_init(&tijab, CC_TAMPS, 0, 12, 15, 12, 17, 0, "New tijab");
    dpd_buf4_axpy(&t2a, &tijab, 1);
    dpd_buf4_close(&t2a);
    dpd_buf4_close(&W);
    dpd_buf4_close(&tijab);

    /*** AB ***/
    dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 22, 28, 22, 28, 0, "New tIjAb");
    dpd_buf4_init(&W, CC2_HET1, 0, 24, 22, 24, 22, 0, "CC2 WMbIj");
    dpd_contract244(&tIA, &W, &tIjAb, 0, 0, 1, -1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&tIjAb);
    dpd_buf4_init(&Z, CC_TMP0, 0, 23, 29, 23, 29, 0, "CC2 ZjIbA");
    dpd_buf4_init(&W, CC2_HET1, 0, 27, 23, 27, 23, 0, "CC2 WmBiJ (mB,iJ)");
    dpd_contract244(&tia, &W, &Z, 0, 0, 1, -1, 0);
    dpd_buf4_close(&W);
    dpd_buf4_sort_axpy(&Z, CC_TAMPS, qpsr, 22, 28, "New tIjAb", 1);
    dpd_buf4_close(&Z);

    dpd_file2_close(&tIA);
    dpd_file2_close(&tia);

  }

}
}} // namespace psi::ccenergy
