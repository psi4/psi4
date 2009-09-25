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

void cc2_faeT2(void) {

  dpdfile2 fme, fME, Fae, FAE, fAE, tIA, tia;
  dpdbuf4 tIjAb, tIJAB, tijab, t2;
  dpdbuf4 newtIjAb, newtIJAB, newtijab;
  dpdbuf4 Zijab;

  if(params.ref == 0) { /** RHF **/
    dpd_file2_init(&fAE, CC_OEI, 0, 1, 1, "fAB");

    dpd_buf4_init(&Zijab, CC_TMP0, 0, 0, 5, 0, 5, 0, "CC2 ZIjAb");
    dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_contract424(&tIjAb, &fAE, &Zijab, 3, 1, 0, 1, 0);
    dpd_buf4_close(&tIjAb);
    dpd_buf4_init(&newtIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    dpd_buf4_axpy(&Zijab, &newtIjAb, 1);
    dpd_buf4_close(&newtIjAb);
    dpd_buf4_sort_axpy(&Zijab, CC_TAMPS, qpsr, 0, 5, "New tIjAb", 1);
    dpd_buf4_close(&Zijab);

    dpd_file2_close(&fAE);
  }
  else if(params.ref == 1) { /** ROHF **/

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&fME, CC_OEI, 0, 0, 1, "fIA");
    dpd_file2_init(&FAE, CC2_HET1, 0, 1, 1, "CC2 FAE");
    dpd_contract222(&tIA, &fME, &FAE, 1, 1, -1, 0);
    dpd_file2_close(&FAE);  
    dpd_file2_close(&fME);  
    dpd_file2_close(&tIA);

    dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");
    dpd_file2_init(&fme, CC_OEI, 0, 0, 1, "fia");
    dpd_file2_init(&Fae, CC2_HET1, 0, 1, 1, "CC2 Fae");
    dpd_contract222(&tia, &fme, &Fae, 1, 1, -1, 0);
    dpd_file2_close(&Fae);
    dpd_file2_close(&fme);
    dpd_file2_close(&tia);

    /** F -> tijab **/
    dpd_file2_init(&FAE, CC2_HET1, 0, 1, 1, "CC2 FAE");
    dpd_file2_init(&Fae, CC2_HET1, 0, 1, 1, "CC2 Fae");

    /*** AA ***/
    dpd_buf4_init(&tIJAB, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    dpd_buf4_init(&t2, CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,AB)");
    dpd_contract424(&tIJAB, &FAE, &t2, 3, 1, 0, 1, 0);
    dpd_contract244(&FAE, &tIJAB, &t2, 1, 2, 1, 1, 1);
    dpd_buf4_init(&newtIJAB, CC_TAMPS, 0, 2, 5, 2, 7, 0, "New tIJAB");
    dpd_buf4_axpy(&t2, &newtIJAB, 1);
    dpd_buf4_close(&newtIJAB);
    dpd_buf4_close(&t2);
    dpd_buf4_close(&tIJAB);

    /*** BB ***/
    dpd_buf4_init(&tijab, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tijab");
    dpd_buf4_init(&t2, CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,AB)");
    dpd_contract424(&tijab, &Fae, &t2, 3, 1, 0, 1, 0);
    dpd_contract244(&Fae, &tijab, &t2, 1, 2, 1, 1, 1);
    dpd_buf4_init(&newtijab, CC_TAMPS, 0, 2, 5, 2, 7, 0, "New tijab");
    dpd_buf4_axpy(&t2, &newtijab, 1);
    dpd_buf4_close(&newtijab);
    dpd_buf4_close(&t2);
    dpd_buf4_close(&tijab);

    /*** AB ***/
    dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_buf4_init(&newtIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    dpd_contract424(&tIjAb, &Fae, &newtIjAb, 3, 1, 0, 1, 1);
    dpd_contract244(&FAE, &tIjAb, &newtIjAb, 1, 2, 1, 1, 1);
    dpd_buf4_close(&newtIjAb);
    dpd_buf4_close(&tIjAb);

    dpd_file2_close(&FAE);  
    dpd_file2_close(&Fae);
  }
  else if(params.ref == 2) { /** UHF **/

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&fME, CC_OEI, 0, 0, 1, "fIA");
    dpd_file2_init(&FAE, CC2_HET1, 0, 1, 1, "CC2 FAE");
    dpd_contract222(&tIA, &fME, &FAE, 1, 1, -1, 0);
    dpd_file2_close(&FAE);  
    dpd_file2_close(&fME);  
    dpd_file2_close(&tIA);

    dpd_file2_init(&tia, CC_OEI, 0, 2, 3, "tia");
    dpd_file2_init(&fme, CC_OEI, 0, 2, 3, "fia");
    dpd_file2_init(&Fae, CC2_HET1, 0, 3, 3, "CC2 Fae");
    dpd_contract222(&tia, &fme, &Fae, 1, 1, -1, 0);
    dpd_file2_close(&Fae);
    dpd_file2_close(&fme);
    dpd_file2_close(&tia);

    /** F -> tijab **/
    dpd_file2_init(&FAE, CC2_HET1, 0, 1, 1, "CC2 FAE");
    dpd_file2_init(&Fae, CC2_HET1, 0, 3, 3, "CC2 Fae");

    /*** AA ***/
    dpd_buf4_init(&tIJAB, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    dpd_buf4_init(&t2, CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,AB)");
    dpd_contract424(&tIJAB, &FAE, &t2, 3, 1, 0, 1, 0);
    dpd_contract244(&FAE, &tIJAB, &t2, 1, 2, 1, 1, 1);
    dpd_buf4_close(&tIJAB);
    dpd_buf4_init(&newtIJAB, CC_TAMPS, 0, 2, 5, 2, 7, 0, "New tIJAB");
    dpd_buf4_axpy(&t2, &newtIJAB, 1);
    dpd_buf4_close(&newtIJAB);
    dpd_buf4_close(&t2);

    /*** BB ***/
    dpd_buf4_init(&tijab, CC_TAMPS, 0, 12, 15, 12, 17, 0, "tijab");
    dpd_buf4_init(&t2, CC_TMP0, 0, 12, 15, 12, 15, 0, "T (i>j,ab)");
    dpd_contract424(&tijab, &Fae, &t2, 3, 1, 0, 1, 0);
    dpd_contract244(&Fae, &tijab, &t2, 1, 2, 1, 1, 1);
    dpd_buf4_close(&tijab);
    dpd_buf4_init(&newtijab, CC_TAMPS, 0, 12, 15, 12, 17, 0, "New tijab");
    dpd_buf4_axpy(&t2, &newtijab, 1);
    dpd_buf4_close(&newtijab);
    dpd_buf4_close(&t2);

    /*** AB ***/
    dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_buf4_init(&newtIjAb, CC_TAMPS, 0, 22, 28, 22, 28, 0, "New tIjAb");
    dpd_contract424(&tIjAb, &Fae, &newtIjAb, 3, 1, 0, 1, 1);
    dpd_contract244(&FAE, &tIjAb, &newtIjAb, 1, 2, 1, 1, 1);
    dpd_buf4_close(&newtIjAb);
    dpd_buf4_close(&tIjAb);

    dpd_file2_close(&FAE);  
    dpd_file2_close(&Fae);
  }
}
}} // namespace psi::ccenergy
