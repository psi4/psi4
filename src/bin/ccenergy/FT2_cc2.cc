/*! \file
    \ingroup CCENERGY
    \brief Enter brief description of file here 
*/
#include <libdpd/dpd.h>
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccenergy {

void FT2_CC2(void)
{
  dpdbuf4 newT2, T2, Z;
  dpdfile2 F;

  if(params.ref == 0) { /* RHF */
    dpd_buf4_init(&newT2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_file2_init(&F, CC_OEI, 0, 1, 1, "fAB");
    dpd_contract424(&T2, &F, &newT2, 3, 1, 0, 1, 1);
    dpd_contract244(&F, &T2, &newT2, 1, 2, 1, 1, 1);
    dpd_file2_close(&F);
    dpd_file2_init(&F, CC_OEI, 0, 0, 0, "fIJ");
    dpd_contract424(&T2, &F, &newT2, 1, 0, 1, -1, 1);
    dpd_contract244(&F, &T2, &newT2, 0, 0, 0, -1, 1);
    dpd_file2_close(&F);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&newT2);
  }
  else if(params.ref == 1) { /* ROHF */
    dpd_buf4_init(&newT2, CC_TAMPS, 0, 2, 5, 2, 7, 0, "New tIJAB");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    dpd_file2_init(&F, CC_OEI, 0, 1, 1, "fAB");
    dpd_buf4_init(&Z, CC_TMP0, 0, 2, 5, 2, 5, 0, "Z(I>J,AB)");
    dpd_contract424(&T2, &F, &Z, 3, 1, 0, 1, 0);
    dpd_contract244(&F, &T2, &Z, 1, 2, 1, 1, 1);
    dpd_buf4_axpy(&Z, &newT2, 1);
    dpd_buf4_close(&Z);
    dpd_file2_close(&F);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&newT2);

    dpd_buf4_init(&newT2, CC_TAMPS, 0, 2, 5, 2, 7, 0, "New tijab");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tijab");
    dpd_file2_init(&F, CC_OEI, 0, 1, 1, "fab");
    dpd_buf4_init(&Z, CC_TMP0, 0, 2, 5, 2, 5, 0, "Z(I>J,AB)");
    dpd_contract424(&T2, &F, &Z, 3, 1, 0, 1, 0);
    dpd_contract244(&F, &T2, &Z, 1, 2, 1, 1, 1);
    dpd_buf4_axpy(&Z, &newT2, 1);
    dpd_buf4_close(&Z);
    dpd_file2_close(&F);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&newT2);

    dpd_buf4_init(&newT2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_file2_init(&F, CC_OEI, 0, 1, 1, "fab");
    dpd_contract424(&T2, &F, &newT2, 3, 1, 0, 1, 1);
    dpd_file2_close(&F);
    dpd_file2_init(&F, CC_OEI, 0, 1, 1, "fAB");
    dpd_contract244(&F, &T2, &newT2, 1, 2, 1, 1, 1);
    dpd_file2_close(&F);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&newT2);
  }
  else if(params.ref == 2) { /* UHF */

    dpd_buf4_init(&newT2, CC_TAMPS, 0, 2, 5, 2, 7, 0, "New tIJAB");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    dpd_file2_init(&F, CC_OEI, 0, 1, 1, "fAB");
    dpd_buf4_init(&Z, CC_TMP0, 0, 2, 5, 2, 5, 0, "Z(I>J,AB)");
    dpd_contract424(&T2, &F, &Z, 3, 1, 0, 1, 0);
    dpd_contract244(&F, &T2, &Z, 1, 2, 1, 1, 1);
    dpd_buf4_axpy(&Z, &newT2, 1);
    dpd_buf4_close(&Z);
    dpd_file2_close(&F);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&newT2);

    dpd_buf4_init(&newT2, CC_TAMPS, 0, 12, 15, 12, 17, 0, "New tijab");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 12, 15, 12, 17, 0, "tijab");
    dpd_file2_init(&F, CC_OEI, 0, 3, 3, "fab");
    dpd_buf4_init(&Z, CC_TMP0, 0, 12, 15, 12, 15, 0, "Z(i>j,ab)");
    dpd_contract424(&T2, &F, &Z, 3, 1, 0, 1, 0);
    dpd_contract244(&F, &T2, &Z, 1, 2, 1, 1, 1);
    dpd_buf4_axpy(&Z, &newT2, 1);
    dpd_buf4_close(&Z);
    dpd_file2_close(&F);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&newT2);

    dpd_buf4_init(&newT2, CC_TAMPS, 0, 22, 28, 22, 28, 0, "New tIjAb");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_file2_init(&F, CC_OEI, 0, 3, 3, "fab");
    dpd_contract424(&T2, &F, &newT2, 3, 1, 0, 1, 1);
    dpd_file2_close(&F);
    dpd_file2_init(&F, CC_OEI, 0, 1, 1, "fAB");
    dpd_contract244(&F, &T2, &newT2, 1, 2, 1, 1, 1);
    dpd_file2_close(&F);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&newT2);

  }
}
}} // namespace psi::ccenergy
