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

void ZT2(void)
{
  dpdbuf4 ZIJMA, ZIJAM, Zijma, Zijam, ZIjMa, ZIjAm, Z;
  dpdbuf4 newtIJAB, newtijab, newtIjAb, T2;
  dpdfile2 tIA, tia, T1;
  dpdbuf4 t2, X;

  if(params.ref == 0) { /** RHF **/
    dpd_buf4_init(&X, CC_TMP0, 0, 5, 0, 5, 0, 0, "X(Ab,Ij)");
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_buf4_init(&Z, CC_MISC, 0, 10, 0, 10, 0, 0, "ZMbIj");
    dpd_contract244(&T1, &Z, &X, 0, 0, 0, -1, 0);
    dpd_buf4_close(&Z); 
    dpd_file2_close(&T1); 
    dpd_buf4_sort_axpy(&X, CC_TAMPS, rspq, 0, 5, "New tIjAb", 1);
    dpd_buf4_sort_axpy(&X, CC_TAMPS, srqp, 0, 5, "New tIjAb", 1);
    dpd_buf4_close(&X);
  }
  else if(params.ref == 1) { /** ROHF **/
    dpd_buf4_init(&ZIJMA, CC_MISC, 0, 2, 10, 2, 10, 0, "ZIJMA");
    dpd_buf4_init(&ZIJAM, CC_MISC, 0, 2, 11, 2, 11, 0, "ZIJAM");
    dpd_buf4_init(&Zijma, CC_MISC, 0, 2, 10, 2, 10, 0, "Zijma");
    dpd_buf4_init(&Zijam, CC_MISC, 0, 2, 11, 2, 11, 0, "Zijam");
    dpd_buf4_init(&ZIjMa, CC_MISC, 0, 0, 10, 0, 10, 0, "ZIjMa");
    dpd_buf4_init(&ZIjAm, CC_MISC, 0, 0, 11, 0, 11, 0, "ZIjAm");

    dpd_buf4_init(&newtIJAB, CC_TAMPS, 0, 2, 5, 2, 7, 0, "New tIJAB");
    dpd_buf4_init(&newtijab, CC_TAMPS, 0, 2, 5, 2, 7, 0, "New tijab");
    dpd_buf4_init(&newtIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");

    dpd_buf4_init(&t2, CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,AB)");
    dpd_contract424(&ZIJAM, &tIA, &t2, 3, 0, 0, 1, 0);
    dpd_contract244(&tIA, &ZIJMA, &t2, 0, 2, 1, -1, 1);
    dpd_buf4_axpy(&t2, &newtIJAB, 1);
    dpd_buf4_close(&t2);

    dpd_buf4_init(&t2, CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,AB)");
    dpd_contract424(&Zijam, &tia, &t2, 3, 0, 0, 1, 0);
    dpd_contract244(&tia, &Zijma, &t2, 0, 2, 1, -1, 1);
    dpd_buf4_axpy(&t2, &newtijab, 1);
    dpd_buf4_close(&t2);

    dpd_contract424(&ZIjAm, &tia, &newtIjAb, 3, 0, 0, -1, 1);
    dpd_contract244(&tIA, &ZIjMa, &newtIjAb, 0, 2, 1, -1, 1);

    dpd_file2_close(&tIA); 
    dpd_file2_close(&tia);

    dpd_buf4_close(&newtIJAB); 
    dpd_buf4_close(&newtijab); 
    dpd_buf4_close(&newtIjAb); 

    dpd_buf4_close(&ZIJMA); 
    dpd_buf4_close(&ZIJAM); 
    dpd_buf4_close(&Zijma);
    dpd_buf4_close(&Zijam); 
    dpd_buf4_close(&ZIjMa); 
    dpd_buf4_close(&ZIjAm);
  }
  else if(params.ref == 2) { /*** UHF ***/

    dpd_buf4_init(&ZIJMA, CC_MISC, 0, 2, 20, 2, 20, 0, "ZIJMA");
    dpd_buf4_init(&ZIJAM, CC_MISC, 0, 2, 21, 2, 21, 0, "ZIJAM");
    dpd_buf4_init(&Zijma, CC_MISC, 0, 12, 30, 12, 30, 0, "Zijma");
    dpd_buf4_init(&Zijam, CC_MISC, 0, 12, 31, 12, 31, 0, "Zijam");
    dpd_buf4_init(&ZIjMa, CC_MISC, 0, 22, 24, 22, 24, 0, "ZIjMa");
    dpd_buf4_init(&ZIjAm, CC_MISC, 0, 22, 26, 22, 26, 0, "ZIjAm");

    dpd_buf4_init(&newtIJAB, CC_TAMPS, 0, 2, 5, 2, 7, 0, "New tIJAB");
    dpd_buf4_init(&newtijab, CC_TAMPS, 0, 12, 15, 12, 17, 0, "New tijab");
    dpd_buf4_init(&newtIjAb, CC_TAMPS, 0, 22, 28, 22, 28, 0, "New tIjAb");

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 2, 3, "tia");

    dpd_buf4_init(&t2, CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,AB)");
    dpd_contract424(&ZIJAM, &tIA, &t2, 3, 0, 0, 1, 0);
    dpd_contract244(&tIA, &ZIJMA, &t2, 0, 2, 1, -1, 1);
    dpd_buf4_axpy(&t2, &newtIJAB, 1);
    dpd_buf4_close(&t2);

    dpd_buf4_init(&t2, CC_TMP0, 0, 12, 15, 12, 15, 0, "T (i>j,ab)");
    dpd_contract424(&Zijam, &tia, &t2, 3, 0, 0, 1, 0);
    dpd_contract244(&tia, &Zijma, &t2, 0, 2, 1, -1, 1);
    dpd_buf4_axpy(&t2, &newtijab, 1);
    dpd_buf4_close(&t2);

    dpd_contract424(&ZIjAm, &tia, &newtIjAb, 3, 0, 0, -1, 1);
    dpd_contract244(&tIA, &ZIjMa, &newtIjAb, 0, 2, 1, -1, 1);

    dpd_file2_close(&tIA); 
    dpd_file2_close(&tia);

    dpd_buf4_close(&newtIJAB); 
    dpd_buf4_close(&newtijab); 
    dpd_buf4_close(&newtIjAb); 

    dpd_buf4_close(&ZIJMA); 
    dpd_buf4_close(&ZIJAM); 
    dpd_buf4_close(&Zijma);
    dpd_buf4_close(&Zijam); 
    dpd_buf4_close(&ZIjMa); 
    dpd_buf4_close(&ZIjAm);

  }
}
}} // namespace psi::ccenergy
