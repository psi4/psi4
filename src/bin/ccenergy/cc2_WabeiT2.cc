/*! \file
    \ingroup CCENERGY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <libdpd/dpd.h>
#include <libpsio/psio.h>
#include <libqt/qt.h>
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccenergy {

void cc2_WabeiT2(void) {

  int rowx, colx, rowz, colz, ab;
  int GX, GZ, Ge, Gi, Gj, hxbuf, hzbuf;
  dpdfile2 t1, tia, tIA;
  dpdbuf4 Z, W, X;
  dpdbuf4 t2, t2a, t2b, tIJAB, tijab, tIjAb;

  if(params.ref == 0) { /** RHF **/

    dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");

    dpd_buf4_init(&Z, CC_TMP0, 0, 5, 0, 5, 0, 0, "CC2 ZAbIj");
    dpd_buf4_init(&X, CC2_HET1, 0, 5, 11, 5, 11, 0, "CC2 WAbEi");
    dpd_contract244(&t1, &X, &Z, 1, 2, 1, 1, 0);

    dpd_buf4_sort_axpy(&Z, CC_TAMPS, rspq, 0, 5, "New tIjAb", 1);
    dpd_buf4_sort_axpy(&Z, CC_TAMPS, srqp, 0, 5, "New tIjAb", 1);
    dpd_buf4_close(&Z);

    dpd_file2_close(&t1);
  }

  else if(params.ref == 1) { /** ROHF **/

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");

    /*** AA ***/
    dpd_buf4_init(&W, CC2_HET1, 0, 11, 7, 11, 7, 0, "CC2 WABEI (EI,A>B)");
    dpd_buf4_init(&t2, CC_TMP0, 0, 0, 7, 0, 7, 0, "T (IJ,A>B)");
    dpd_contract244(&tIA, &W, &t2, 1, 0, 0, 1, 0);
    dpd_buf4_sort(&t2, CC_TMP0, qprs, 0, 7, "T (JI,A>B)");
    dpd_buf4_close(&t2);
    dpd_buf4_init(&t2a, CC_TMP0, 0, 0, 7, 0, 7, 0, "T (IJ,A>B)");
    dpd_buf4_init(&t2b, CC_TMP0, 0, 0, 7, 0, 7, 0, "T (JI,A>B)");
    dpd_buf4_axpy(&t2b, &t2a, -1);
    dpd_buf4_init(&tIJAB, CC_TAMPS, 0, 0, 7, 2, 7, 0, "New tIJAB");
    dpd_buf4_axpy(&t2a, &tIJAB, 1);
    dpd_buf4_close(&tIJAB);
    dpd_buf4_close(&t2b);
    dpd_buf4_close(&t2a);
    dpd_buf4_close(&W);

    /*** BB ***/
    dpd_buf4_init(&W, CC2_HET1, 0, 11, 7, 11, 7, 0, "CC2 Wabei (ei,a>b)");
    dpd_buf4_init(&t2, CC_TMP0, 0, 0, 7, 0, 7, 0, "T (IJ,A>B)");
    dpd_contract244(&tia, &W, &t2, 1, 0, 0, 1, 0);
    dpd_buf4_sort(&t2, CC_TMP0, qprs, 0, 7, "T (JI,A>B)");
    dpd_buf4_close(&t2);
    dpd_buf4_init(&t2a, CC_TMP0, 0, 0, 7, 0, 7, 0, "T (IJ,A>B)");
    dpd_buf4_init(&t2b, CC_TMP0, 0, 0, 7, 0, 7, 0, "T (JI,A>B)");
    dpd_buf4_axpy(&t2b, &t2a, -1);
    dpd_buf4_init(&tijab, CC_TAMPS, 0, 0, 7, 2, 7, 0, "New tijab");
    dpd_buf4_axpy(&t2a, &tijab, 1);
    dpd_buf4_close(&tijab);
    dpd_buf4_close(&t2b);
    dpd_buf4_close(&t2a);
    dpd_buf4_close(&W);

    /*** AB ***/
    dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    dpd_buf4_init(&W, CC2_HET1, 0, 11, 5, 11, 5, 0, "CC2 WAbEi (Ei,Ab)");
    dpd_contract244(&tIA, &W, &tIjAb, 1, 0, 0, 1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&tIjAb);
    dpd_buf4_init(&Z, CC_TMP0, 0, 0, 5, 0, 5, 0, "CC2 ZjIbA");
    dpd_buf4_init(&W, CC2_HET1, 0, 11, 5, 11, 5, 0, "CC2 WaBeI (eI,aB)");
    dpd_contract244(&tia, &W, &Z, 1, 0, 0, 1, 0);
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
    dpd_buf4_init(&W, CC2_HET1, 0, 21, 7, 21, 7, 0, "CC2 WABEI (EI,A>B)");
    dpd_buf4_init(&t2, CC_TMP0, 0, 0, 7, 0, 7, 0, "T (IJ,A>B)");
    dpd_contract244(&tIA, &W, &t2, 1, 0, 0, 1, 0);
    dpd_buf4_sort(&t2, CC_TMP0, qprs, 0, 7, "T (JI,A>B)");
    dpd_buf4_close(&t2);
    dpd_buf4_init(&t2a, CC_TMP0, 0, 0, 7, 0, 7, 0, "T (IJ,A>B)");
    dpd_buf4_init(&t2b, CC_TMP0, 0, 0, 7, 0, 7, 0, "T (JI,A>B)");
    dpd_buf4_axpy(&t2b, &t2a, -1);
    dpd_buf4_init(&tIJAB, CC_TAMPS, 0, 0, 7, 2, 7, 0, "New tIJAB");
    dpd_buf4_axpy(&t2a, &tIJAB, 1);
    dpd_buf4_close(&t2b);
    dpd_buf4_close(&t2a);
    dpd_buf4_close(&W);
    dpd_buf4_close(&tIJAB);

    /*** BB ***/
    dpd_buf4_init(&W, CC2_HET1, 0, 31, 17, 31, 17, 0, "CC2 Wabei (ei,a>b)");
    dpd_buf4_init(&t2, CC_TMP0, 0, 10, 17, 10, 17, 0, "T (ij,a>b)");
    dpd_contract244(&tia, &W, &t2, 1, 0, 0, 1, 0);
    dpd_buf4_sort(&t2, CC_TMP0, qprs, 10, 17, "T (ji,a>b)");
    dpd_buf4_close(&t2);
    dpd_buf4_init(&t2a, CC_TMP0, 0, 10, 17, 10, 17, 0, "T (ij,a>b)");
    dpd_buf4_init(&t2b, CC_TMP0, 0, 10, 17, 10, 17, 0, "T (ji,a>b)");
    dpd_buf4_axpy(&t2b, &t2a, -1);
    dpd_buf4_init(&tijab, CC_TAMPS, 0, 10, 17, 12, 17, 0, "New tijab");
    dpd_buf4_axpy(&t2a, &tijab, 1);
    dpd_buf4_close(&t2b);
    dpd_buf4_close(&t2a);
    dpd_buf4_close(&W);
    dpd_buf4_close(&tijab);

    /*** AB ***/
    dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 22, 28, 22, 28, 0, "New tIjAb");
    dpd_buf4_init(&W, CC2_HET1, 0, 26, 28, 26, 28, 0, "CC2 WAbEi (Ei,Ab)");
    dpd_contract244(&tIA, &W, &tIjAb, 1, 0, 0, 1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&tIjAb);
    dpd_buf4_init(&Z, CC_TMP0, 0, 23, 29, 23, 29, 0, "CC2 ZjIbA");
    dpd_buf4_init(&W, CC2_HET1, 0, 25, 29, 25, 29, 0, "CC2 WaBeI (eI,aB)");
    dpd_contract244(&tia, &W, &Z, 1, 0, 0, 1, 0);
    dpd_buf4_close(&W);
    dpd_buf4_sort_axpy(&Z, CC_TAMPS, qpsr, 22, 28, "New tIjAb", 1);
    dpd_buf4_close(&Z);

    dpd_file2_close(&tIA); 
    dpd_file2_close(&tia);

  }

}
}} // namespace psi::ccenergy
