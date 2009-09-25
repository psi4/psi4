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

void cc2_Gai_build(int L_irr) {
  dpdbuf4 tIJAB, tijab, tiJaB, tIjAb, tijAB, tIJab, t2;
  dpdfile2 G, GAI, Gai, L1, LIA, Lia;
  dpdbuf4 LIJAB, Lijab, LIjAb, LiJaB;

    if(params.ref == 0) {
      dpd_file2_init(&G, CC_TMP0, L_irr, 1, 0, "CC2 GAI");

      dpd_file2_init(&L1, CC_LAMBDA, L_irr, 0, 1, "LIA");
      dpd_buf4_init(&t2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "2 tIAjb - tIBja");
      dpd_contract422(&t2, &L1, &G, 0, 1, 1, 0);
      dpd_buf4_close(&t2);
      dpd_file2_close(&L1);

      dpd_file2_close(&G);
    }

    else if(params.ref == 1) { /** ROHF **/

      dpd_file2_init(&G, CC_TMP0, L_irr, 1, 0, "GAI");
      dpd_file2_init(&G, CC_TMP0, L_irr, 4, 3, "Gai");

      /* T2(MJ,AB) * L2(IJ,AB) --> G(M,I) */
      dpd_buf4_init(&tIJAB, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
      dpd_buf4_init(&LIJAB, CC_LAMBDA, L_irr, 0, 7, 2, 7, 0, "LIJAB");
      dpd_contract442(&tIJAB, &LIJAB, &G, 0, 0, 1.0, 0.0);
      dpd_buf4_close(&tIJAB);
      dpd_buf4_close(&LIJAB);  

      /* T2(Mj,Ab) * L2(Ij,Ab) --> G(M,I) */
      dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
      dpd_buf4_init(&LIjAb, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
      dpd_contract442(&tIjAb, &LIjAb, &G, 0, 0, 1.0, 1.0);
      dpd_buf4_close(&tIjAb);
      dpd_buf4_close(&LIjAb);

      /* T2(mj,ab) * L2(ij,ab) --> G(m,i) */
      dpd_buf4_init(&tijab, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tijab");
      dpd_buf4_init(&Lijab, CC_LAMBDA, L_irr, 0, 7, 2, 7, 0, "Lijab");
      dpd_contract442(&tijab, &Lijab, &G, 0, 0, 1.0, 0.0);
      dpd_buf4_close(&tijab);
      dpd_buf4_close(&Lijab); 

      /* T2(mJ,aB) * L2(iJ,aB) --> G(m,i) */
      dpd_buf4_init(&tiJaB, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
      dpd_buf4_init(&LiJaB, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LiJaB");
      dpd_contract442(&tiJaB, &LiJaB, &G, 0, 0, 1.0, 1.0);
      dpd_buf4_close(&tiJaB);
      dpd_buf4_close(&LiJaB);

      dpd_file2_close(&G);
      dpd_file2_close(&G);
    }

    else if(params.ref == 2) { /** UHF **/

      dpd_file2_init(&GAI, CC_TMP0, L_irr, 1, 0, "CC2 GAI");
      dpd_file2_init(&Gai, CC_TMP0, L_irr, 3, 2, "CC2 Gai");

      /** AA **/
      dpd_buf4_init(&tIJAB, CC_TAMPS, 0, 20, 20, 20, 20, 0, "tIAJB");
      dpd_file2_init(&LIA, CC_LAMBDA, L_irr, 0, 1, "LIA");
      dpd_contract422(&tIJAB, &LIA, &GAI, 0, 1, 1, 0);
      dpd_file2_close(&LIA);
      dpd_buf4_close(&tIJAB);

      dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
      dpd_file2_init(&Lia, CC_LAMBDA, L_irr, 2, 3, "Lia");
      dpd_contract422(&tIjAb, &Lia, &GAI, 0, 1, 1, 1);
      dpd_file2_close(&Lia);
      dpd_buf4_close(&tIjAb);

      /** BB **/
      dpd_buf4_init(&tijab, CC_TAMPS, 0, 30, 30, 30, 30, 0, "tiajb");
      dpd_file2_init(&Lia, CC_LAMBDA, L_irr, 2, 3, "Lia");
      dpd_contract422(&tijab, &Lia, &Gai, 0, 1, 1, 0);
      dpd_file2_close(&Lia);
      dpd_buf4_close(&tijab);

      dpd_buf4_init(&tiJaB, CC_TAMPS, 0, 30, 20, 30, 20, 0, "tiaJB");
      dpd_file2_init(&LIA, CC_LAMBDA, L_irr, 0, 1, "LIA");
      dpd_contract422(&tiJaB, &LIA, &Gai, 0, 1, 1, 1);
      dpd_file2_close(&LIA);
      dpd_buf4_close(&tiJaB);

      dpd_file2_close(&Gai);
      dpd_file2_close(&GAI);
    }

  return;
}



}} // namespace psi::cclambda
