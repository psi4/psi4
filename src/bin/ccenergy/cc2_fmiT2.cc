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

void cc2_fmiT2(void) {

  dpdfile2 fia, fIA, Fmi, FMI, fMI, tIA, tia;
  dpdbuf4 tIjAb, tIJAB, tijab, t2;
  dpdbuf4 newtIjAb, newtIJAB, newtijab;
  dpdbuf4 Zijab;

  if(params.ref == 0) { /** RHF **/
    dpd_file2_init(&fMI, CC_OEI, 0, 0, 0, "fIJ");

    dpd_buf4_init(&Zijab, CC_TMP0, 0, 0, 5, 0, 5, 0, "CC2 ZIjAb");
    dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_contract244(&fMI, &tIjAb, &Zijab, 0, 0, 0, -1, 0);
    dpd_buf4_close(&tIjAb);
    dpd_buf4_init(&newtIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    dpd_buf4_axpy(&Zijab, &newtIjAb, 1);
    dpd_buf4_close(&newtIjAb);
    dpd_buf4_sort_axpy(&Zijab, CC_TAMPS, qpsr, 0, 5, "New tIjAb", 1);
    dpd_buf4_close(&Zijab);

    dpd_file2_close(&fMI);
  }
  else if(params.ref == 1) { /** ROHF **/

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&fIA, CC_OEI, 0, 0, 1, "fIA");
    dpd_file2_init(&FMI, CC2_HET1, 0, 0, 0, "CC2 FMI");
    dpd_contract222(&fIA, &tIA, &FMI, 0, 0, 1, 0);
    dpd_file2_close(&FMI); 
    dpd_file2_close(&fIA); 
    dpd_file2_close(&tIA);

    dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");
    dpd_file2_init(&fia, CC_OEI, 0, 0, 1, "fia");
    dpd_file2_init(&Fmi, CC2_HET1, 0, 0, 0, "CC2 Fmi");
    dpd_contract222(&fia, &tia, &Fmi, 0, 0, 1, 0);
    dpd_file2_close(&Fmi);
    dpd_file2_close(&fia); 
    dpd_file2_close(&tia);

    /** F -> tijab **/
    dpd_file2_init(&FMI, CC2_HET1, 0, 0, 0, "CC2 FMI");
    dpd_file2_init(&Fmi, CC2_HET1, 0, 0, 0, "CC2 Fmi");

    /*** AA ***/
    dpd_buf4_init(&tIJAB, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
    dpd_buf4_init(&t2, CC_TMP0, 0, 0, 7, 0, 7, 0, "T (IJ,A>B)");
    dpd_contract424(&tIJAB, &FMI, &t2, 1, 0, 1, -1, 0);
    dpd_contract244(&FMI, &tIJAB, &t2, 0, 0, 0, -1, 1);
    dpd_buf4_init(&newtIJAB, CC_TAMPS, 0, 0, 7, 2, 7, 0, "New tIJAB");
    dpd_buf4_axpy(&t2, &newtIJAB, 1);
    dpd_buf4_close(&newtIJAB);
    dpd_buf4_close(&t2);
    dpd_buf4_close(&tIJAB);

    /*** BB ***/
    dpd_buf4_init(&tijab, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tijab");
    dpd_buf4_init(&t2, CC_TMP0, 0, 0, 7, 0, 7, 0, "T (IJ,A>B)");
    dpd_contract424(&tijab, &Fmi, &t2, 1, 0, 1, -1, 0);
    dpd_contract244(&Fmi, &tijab, &t2, 0, 0, 0, -1, 1);
    dpd_buf4_init(&newtijab, CC_TAMPS, 0, 0, 7, 2, 7, 0, "New tijab");
    dpd_buf4_axpy(&t2, &newtijab, 1);
    dpd_buf4_close(&newtijab);
    dpd_buf4_close(&t2);
    dpd_buf4_close(&tijab);

    /*** AB ***/
    dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_buf4_init(&newtIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    dpd_contract424(&tIjAb, &Fmi, &newtIjAb, 1, 0, 1, -1, 1);
    dpd_contract244(&FMI, &tIjAb, &newtIjAb, 0, 0, 0, -1, 1);
    dpd_buf4_close(&newtIjAb);
    dpd_buf4_close(&tIjAb);

    dpd_file2_close(&FMI); 
    dpd_file2_close(&Fmi);
  }
  else if(params.ref == 2) { /** UHF **/

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&fIA, CC_OEI, 0, 0, 1, "fIA");
    dpd_file2_init(&FMI, CC2_HET1, 0, 0, 0, "CC2 FMI");
    dpd_contract222(&fIA, &tIA, &FMI, 0, 0, 1, 0);
    dpd_file2_close(&FMI); 
    dpd_file2_close(&fIA); 
    dpd_file2_close(&tIA);

    dpd_file2_init(&tia, CC_OEI, 0, 2, 3, "tia");
    dpd_file2_init(&fia, CC_OEI, 0, 2, 3, "fia");
    dpd_file2_init(&Fmi, CC2_HET1, 0, 2, 2, "CC2 Fmi");
    dpd_contract222(&fia, &tia, &Fmi, 0, 0, 1, 0);
    dpd_file2_close(&Fmi);
    dpd_file2_close(&fia);
    dpd_file2_close(&tia);

    /** tijab <- Fmi **/
    dpd_file2_init(&FMI, CC2_HET1, 0, 0, 0, "CC2 FMI");
    dpd_file2_init(&Fmi, CC2_HET1, 0, 2, 2, "CC2 Fmi");

    /*** AA ***/
    dpd_buf4_init(&tIJAB, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
    dpd_buf4_init(&t2, CC_TMP0, 0, 0, 7, 0, 7, 0, "T (IJ,A>B)");
    dpd_contract424(&tIJAB, &FMI, &t2, 1, 0, 1, -1, 0);
    dpd_contract244(&FMI, &tIJAB, &t2, 0, 0, 0, -1, 1);
    dpd_buf4_close(&tIJAB);
    dpd_buf4_init(&newtIJAB, CC_TAMPS, 0, 0, 7, 2, 7, 0, "New tIJAB");
    dpd_buf4_axpy(&t2, &newtIJAB, 1);
    dpd_buf4_close(&newtIJAB);
    dpd_buf4_close(&t2);

    /*** BB ***/
    dpd_buf4_init(&tijab, CC_TAMPS, 0, 10, 17, 12, 17, 0, "tijab");
    dpd_buf4_init(&t2, CC_TMP0, 0, 10, 17, 10, 17, 0, "T (ij,a>b)");
    dpd_contract424(&tijab, &Fmi, &t2, 1, 0, 1, -1, 0);
    dpd_contract244(&Fmi, &tijab, &t2, 0, 0, 0, -1, 1);
    dpd_buf4_close(&tijab);
    dpd_buf4_init(&newtijab, CC_TAMPS, 0, 10, 17, 12, 17, 0, "New tijab");
    dpd_buf4_axpy(&t2, &newtijab, 1);
    dpd_buf4_close(&newtijab);
    dpd_buf4_close(&t2);

    /*** AB ***/
    dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_buf4_init(&newtIjAb, CC_TAMPS, 0, 22, 28, 22, 28, 0, "New tIjAb");
    dpd_contract424(&tIjAb, &Fmi, &newtIjAb, 1, 0, 1, -1, 1);
    dpd_contract244(&FMI, &tIjAb, &newtIjAb, 0, 0, 0, -1, 1);
    dpd_buf4_close(&newtIjAb);
    dpd_buf4_close(&tIjAb);

    dpd_file2_close(&Fmi);
    dpd_file2_close(&FMI); 
  }
}
}} // namespace psi::ccenergy
