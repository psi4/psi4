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

/** The ROHF version of this contraction can be done with fewer contractions.
 **/

void cc2_faeL2(int L_irr)
{
  int h, e;
  dpdbuf4 Lijab, LIJAB, LIjAb;
  dpdbuf4 newLIJAB, newLijab, newLIjAb;
  dpdfile2 fab, fAB, F;
  dpdbuf4 X, X1, X2;
  dpdbuf4 L2, newL2;

  /* RHS += P(ab)*Lijae*Feb */

  if(params.ref == 0) { /** RHF **/

    dpd_buf4_init(&X, CC_TMP0, L_irr, 0, 5, 0, 5, 0, "X(Ij,Ab)");

    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    dpd_file2_init(&F, CC_OEI, 0, 1, 1, "fAB");
    dpd_contract424(&L2, &F, &X, 3, 0, 0, 1, 0);
    dpd_file2_close(&F);
    dpd_buf4_close(&L2);

    dpd_buf4_sort_axpy(&X, CC_LAMBDA, qpsr, 0, 5, "New LIjAb", 1);
    dpd_buf4_init(&newL2, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
    dpd_buf4_axpy(&X, &newL2, 1);
    dpd_buf4_close(&newL2);

    dpd_buf4_close(&X);

  }
  else if(params.ref == 1) { /** ROHF **/

    dpd_file2_init(&fAB, CC_OEI, 0, 1, 1, "fAB");
    dpd_file2_init(&fab, CC_OEI, 0, 1, 1, "fab");

    dpd_buf4_init(&LIJAB, CC_LAMBDA, L_irr, 2, 5, 2, 7, 0, "LIJAB");
    dpd_buf4_init(&X1, CC_TMP1, L_irr, 2, 5, 2, 5, 0, "X(2,5) 1");
    dpd_contract424(&LIJAB, &fAB, &X1, 3, 0, 0, 1.0, 0.0);
    dpd_buf4_init(&X2, CC_TMP1, L_irr, 2, 5, 2, 5, 0, "X(2,5) 2");
    dpd_contract244(&fAB, &LIJAB, &X2, 0, 2, 1, 1.0, 0.0);
    dpd_buf4_close(&LIJAB);
    dpd_buf4_axpy(&X1, &X2, 1.0);
    dpd_buf4_close(&X1);
    dpd_buf4_init(&newLIJAB, CC_LAMBDA, L_irr, 2, 5, 2, 7, 0, "New LIJAB");
    dpd_buf4_axpy(&X2, &newLIJAB, 1.0);
    dpd_buf4_close(&X2);
    dpd_buf4_close(&newLIJAB);

    dpd_buf4_init(&Lijab, CC_LAMBDA, L_irr, 2, 5, 2, 7, 0, "Lijab");
    dpd_buf4_init(&X1, CC_TMP1, L_irr, 2, 5, 2, 5, 0, "X(2,5) 1");
    dpd_contract424(&Lijab, &fab, &X1, 3, 0, 0, 1.0, 0.0);
    dpd_buf4_init(&X2, CC_TMP1, L_irr, 2, 5, 2, 5, 0, "X(2,5) 2");
    dpd_contract244(&fab, &Lijab, &X2, 0, 2, 1, 1.0, 0.0);
    dpd_buf4_close(&Lijab);
    dpd_buf4_axpy(&X1, &X2, 1.0);
    dpd_buf4_close(&X1);
    dpd_buf4_init(&newLijab, CC_LAMBDA, L_irr, 2, 5, 2, 7, 0, "New Lijab");
    dpd_buf4_axpy(&X2, &newLijab, 1.0);
    dpd_buf4_close(&X2);
    dpd_buf4_close(&newLijab);

    dpd_buf4_init(&LIjAb, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    dpd_buf4_init(&newLIjAb, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
    dpd_contract424(&LIjAb, &fab, &newLIjAb, 3, 0, 0, 1.0, 1.0);
    dpd_contract244(&fAB, &LIjAb, &newLIjAb, 0, 2, 1, 1.0, 1.0);
    dpd_buf4_close(&LIjAb);
    dpd_buf4_close(&newLIjAb);

    dpd_file2_close(&fab);
    dpd_file2_close(&fAB);
  }
  else if(params.ref == 2) { /** UHF **/

    dpd_file2_init(&fAB, CC_OEI, 0, 1, 1, "fAB");
    dpd_file2_init(&fab, CC_OEI, 0, 3, 3, "fab");
    dpd_file2_copy(&fAB, CC_OEI, "fAB diag");
    dpd_file2_copy(&fab, CC_OEI, "fab diag");
    dpd_file2_close(&fAB);
    dpd_file2_close(&fab);

    dpd_file2_init(&fAB, CC_OEI, 0, 1, 1, "fAB diag");
    dpd_file2_init(&fab, CC_OEI, 0, 3, 3, "fab diag");
  
    dpd_file2_mat_init(&fAB);
    dpd_file2_mat_rd(&fAB);
    dpd_file2_mat_init(&fab);
    dpd_file2_mat_rd(&fab);
  
    for(h=0; h < moinfo.nirreps; h++) {
      
	for(e=0; e < fAB.params->coltot[h]; e++) 
	  fAB.matrix[h][e][e] = 0;

	for(e=0; e < fab.params->coltot[h]; e++)
	  fab.matrix[h][e][e] = 0;
      
    }

    dpd_file2_mat_wrt(&fAB);
    dpd_file2_mat_close(&fAB);
    dpd_file2_mat_wrt(&fab);
    dpd_file2_mat_close(&fab);

    dpd_file2_close(&fAB);
    dpd_file2_close(&fab);

    dpd_file2_init(&fAB, CC_OEI, 0, 1, 1, "fAB diag");
    dpd_file2_init(&fab, CC_OEI, 0, 3, 3, "fab diag");

    /** X(IJ,AB) = L_IJ^AE F_EB **/
    dpd_buf4_init(&X, CC_TMP1, L_irr, 2, 5, 2, 5, 0, "X(IJ,AB) A");
    dpd_buf4_init(&LIJAB, CC_LAMBDA, L_irr, 2, 5, 2, 7, 0, "LIJAB");
    dpd_contract424(&LIJAB, &fAB, &X, 3, 0, 0, 1, 0);
    dpd_buf4_close(&LIJAB);
    /** X(IJ,AB) --> X'(IJ,BA) **/
    dpd_buf4_sort(&X, CC_TMP1, pqsr, 2, 5, "X'(IJ,BA)");
    dpd_buf4_close(&X);
    /** X(IJ,AB) = X(IJ,AB) - X'(IJ,BA) **/
    dpd_buf4_init(&X1, CC_TMP1, L_irr, 2, 5, 2, 5, 0, "X(IJ,AB) A");
    dpd_buf4_init(&X2, CC_TMP1, L_irr, 2, 5, 2, 5, 0, "X'(IJ,BA)");
    dpd_buf4_axpy(&X2, &X1, -1);
    dpd_buf4_close(&X2);
    dpd_buf4_close(&X1);
    /** L(IJ,AB) <-- X(IJ,AB) **/
    dpd_buf4_init(&X, CC_TMP1, L_irr, 2, 5, 2, 5, 0, "X(IJ,AB) A");
    dpd_buf4_init(&newLIJAB, CC_LAMBDA, L_irr, 2, 5, 2, 7, 0, "New LIJAB");
    dpd_buf4_axpy(&X, &newLIJAB, 1.0);
    dpd_buf4_close(&X);
    dpd_buf4_close(&newLIJAB);

    /** X(ij,ab) = L_ij^ae F_eb **/
    dpd_buf4_init(&X, CC_TMP1, L_irr, 12, 15, 12, 15, 0, "X(ij,ab) A");
    dpd_buf4_init(&LIJAB, CC_LAMBDA, L_irr, 12, 15, 12, 17, 0, "Lijab");
    dpd_contract424(&LIJAB, &fab, &X, 3, 0, 0, 1, 0);
    dpd_buf4_close(&LIJAB);
    /** X(ij,ab) --> X'(ij,ba) **/
    dpd_buf4_sort(&X, CC_TMP1, pqsr, 12, 15, "X'(ij,ba)");
    dpd_buf4_close(&X);
    /** X(ij,ab) = X(ij,ab) - X'(ij,ba) **/
    dpd_buf4_init(&X1, CC_TMP1, L_irr, 12, 15, 12, 15, 0, "X(ij,ab) A");
    dpd_buf4_init(&X2, CC_TMP1, L_irr, 12, 15, 12, 15, 0, "X'(ij,ba)");
    dpd_buf4_axpy(&X2, &X1, -1);
    dpd_buf4_close(&X2);
    dpd_buf4_close(&X1);
    /** L(ij,ab) <-- X(ij,ab) **/
    dpd_buf4_init(&X, CC_TMP1, L_irr, 12, 15, 12, 15, 0, "X(ij,ab) A");
    dpd_buf4_init(&newLIJAB, CC_LAMBDA, L_irr, 12, 15, 12, 17, 0, "New Lijab");
    dpd_buf4_axpy(&X, &newLIJAB, 1.0);
    dpd_buf4_close(&X);
    dpd_buf4_close(&newLIJAB);

    /** L(Ij,Ab) <-- L(Ij,Ae) F(e,b) - F(E,A) L(Ij,Eb) **/
    dpd_buf4_init(&LIjAb, CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "LIjAb");
    dpd_buf4_init(&newLIjAb, CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "New LIjAb");
    dpd_contract424(&LIjAb, &fab, &newLIjAb, 3, 0, 0, 1, 1);
    dpd_contract244(&fAB, &LIjAb, &newLIjAb, 0, 2, 1, 1, 1);
    dpd_buf4_close(&LIjAb);
    dpd_buf4_close(&newLIjAb);

    dpd_file2_close(&fab);
    dpd_file2_close(&fAB);

  }

}

}} // namespace psi::cclambda
