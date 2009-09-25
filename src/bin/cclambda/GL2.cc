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

/* GaeL2(): Computes Gae three-body contributions of HBAR matrix
** elements to the Lambda doubles equations. These are written in
** spin-orbitals as:
**
** L_ij^ab <-- <ij||ae> Gbe - <ij||be> Gae
**
** where Gae = -1/2 t_mn^ef L_mn^af
**
** TDC, July 2002
*/

void GaeL2(int L_irr)
{
  dpdbuf4 L2, newLijab, newLIJAB, newLIjAb, newL2;
  dpdbuf4 D, Z;
  dpdfile2 GAE, Gae, G;
  dpdbuf4 X1, X2;

  /* RHS += P(ab)<ij||ae>Gbe */
  if(params.ref == 0) { /** RHF **/
    dpd_file2_init(&G, CC_LAMBDA, L_irr, 1, 1, "GAE");

    dpd_buf4_init(&Z, CC_TMP0, L_irr, 0, 5, 0, 5, 0, "Z(Ij,Ab)");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_contract424(&D, &G, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&D);
    dpd_buf4_sort_axpy(&Z, CC_LAMBDA, qpsr, 0, 5, "New LIjAb", 1);
    dpd_buf4_init(&newL2, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
    dpd_buf4_axpy(&Z, &newL2, 1);
    dpd_buf4_close(&newL2);
    dpd_buf4_close(&Z);

    dpd_file2_close(&G);
  }
  else if(params.ref == 1) { /** ROHF **/

    dpd_file2_init(&GAE, CC_LAMBDA, L_irr, 1, 1, "GAE");
    dpd_file2_init(&Gae, CC_LAMBDA, L_irr, 1, 1, "Gae");

    dpd_buf4_init(&D, CC_DINTS, 0, 2, 5, 2, 5, 0, "D <ij||ab> (i>j,ab)");
    dpd_buf4_init(&X1, CC_TMP1, L_irr, 2, 5, 2, 5, 0, "X(2,5) 1");
    dpd_contract424(&D, &GAE, &X1, 3, 1, 0, 1.0, 0.0);
    dpd_buf4_init(&X2, CC_TMP1, L_irr, 2, 5, 2, 5, 0, "X(2,5) 2");
    dpd_contract244(&GAE, &D, &X2, 1, 2, 1, 1.0, 0.0);
    dpd_buf4_axpy(&X1, &X2, 1.0);
    dpd_buf4_close(&X1);
    dpd_buf4_init(&newLIJAB, CC_LAMBDA, L_irr, 2, 5, 2, 7, 0, "New LIJAB");
    dpd_buf4_axpy(&X2, &newLIJAB, 1.0);
    dpd_buf4_close(&X2);
    dpd_buf4_close(&newLIJAB);


    dpd_buf4_init(&X1, CC_TMP1, L_irr, 2, 5, 2, 5, 0, "X(2,5) 1");
    dpd_contract424(&D, &Gae, &X1, 3, 1, 0, 1.0, 0.0);
    dpd_buf4_init(&X2, CC_TMP1, L_irr, 2, 5, 2, 5, 0, "X(2,5) 2");
    dpd_contract244(&Gae, &D, &X2, 1, 2, 1, 1.0, 0.0);
    dpd_buf4_close(&D);
    dpd_buf4_axpy(&X1, &X2, 1.0);
    dpd_buf4_close(&X1);
    dpd_buf4_init(&newLijab, CC_LAMBDA, L_irr, 2, 5, 2, 7, 0, "New Lijab");
    dpd_buf4_axpy(&X2, &newLijab, 1.0);
    dpd_buf4_close(&X2);
    dpd_buf4_close(&newLijab);


    dpd_buf4_init(&newLIjAb, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");

    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_contract424(&D, &Gae, &newLIjAb, 3, 1, 0, 1.0, 1.0);
    dpd_contract244(&GAE, &D, &newLIjAb, 1, 2, 1, 1.0, 1.0);
    dpd_buf4_close(&D);

    dpd_buf4_close(&newLIjAb);

    dpd_file2_close(&GAE);
    dpd_file2_close(&Gae);
  }
  else if(params.ref == 2) { /** UHF **/

    dpd_file2_init(&GAE, CC_LAMBDA, L_irr, 1, 1, "GAE");
    dpd_file2_init(&Gae, CC_LAMBDA, L_irr, 3, 3, "Gae");

    /** X(IJ,AB) = <IJ||AE> G(B,E) **/
    dpd_buf4_init(&X1, CC_TMP2, L_irr, 2, 5, 2, 5, 0, "X(IJ,AB)");
    dpd_buf4_init(&D, CC_DINTS, 0, 2, 5, 2, 5, 0, "D <IJ||AB> (I>J,AB)");
    dpd_contract424(&D, &GAE, &X1, 3, 1, 0, 1, 0);
    dpd_buf4_close(&D);
    /** X(IJ,AB) --> X(IJ,BA) **/
    dpd_buf4_sort(&X1, CC_TMP2, pqsr, 2, 5, "X(IJ,BA)");
    /** X(IJ,AB) = X(IJ,AB) - X(IJ,BA) **/
    dpd_buf4_init(&X2, CC_TMP2, L_irr, 2, 5, 2, 5, 0, "X(IJ,BA)");
    dpd_buf4_axpy(&X2, &X1, -1);
    dpd_buf4_close(&X2);
    /** X(IJ,AB) --> New L(IJ,AB) **/
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 2, 5, 2, 7, 0, "New LIJAB");
    dpd_buf4_axpy(&X1, &L2, 1);
    dpd_buf4_close(&X1);
    dpd_buf4_close(&L2);

    /** X(ij,ab) = <ij||ae> G(b,e) **/
    dpd_buf4_init(&X1, CC_TMP2, L_irr, 12, 15, 12, 15, 0, "X(ij,ab)");
    dpd_buf4_init(&D, CC_DINTS, 0, 12, 15, 12, 15, 0, "D <ij||ab> (i>j,ab)");
    dpd_contract424(&D, &Gae, &X1, 3, 1, 0, 1, 0);
    dpd_buf4_close(&D);
    /** X(ij,ab) --> X(ij,ba) **/
    dpd_buf4_sort(&X1, CC_TMP2, pqsr, 12, 15, "X(ij,ba)");
    /** X(ij,ab) = X(ij,ab) - X(ij,ba) **/
    dpd_buf4_init(&X2, CC_TMP2, L_irr, 12, 15, 12, 15, 0, "X(ij,ba)");
    dpd_buf4_axpy(&X2, &X1, -1);
    dpd_buf4_close(&X2);
    /** X(ij,ab) --> New L(ij,ab) **/
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 12, 15, 12, 17, 0, "New Lijab");
    dpd_buf4_axpy(&X1, &L2, 1);
    dpd_buf4_close(&X1);
    dpd_buf4_close(&L2);

    /** New L(Ij,Ab) = <Ij|Ae> G(b,e) + <Ij|Eb> G(A,E) **/
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "New LIjAb");
    dpd_buf4_init(&D, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    dpd_contract424(&D, &Gae, &L2, 3, 1, 0, 1, 1);
    dpd_contract244(&GAE, &D, &L2, 1, 2, 1, 1, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&L2);

    dpd_file2_close(&GAE);
    dpd_file2_close(&Gae);
  }

}

/* GmiL2(): Computes Gmi three-body contributions of HBAR matrix
** elements to the Lambda doubles equations. These are written in
** spin-orbitals as:
**
** L_ij^ab <-- - <im||ab> Gmj + <jm||ab> Gmi
**
** where Gmi = -1/2 t_mn^ef L_in^ef
**
** TDC, July 2002
*/

void GmiL2(int L_irr)
{

  dpdbuf4 L2, newLijab, newLIJAB, newLIjAb, newL2;
  dpdbuf4 D, Z;
  dpdfile2 GMI, Gmi, G;
  dpdbuf4 X1, X2;

  /* RHS -= P(ij) * <im||ab> * Gmj */
  if(params.ref == 0) { /** RHF **/

    dpd_file2_init(&G, CC_LAMBDA, L_irr, 0, 0, "GMI");

    dpd_buf4_init(&Z, CC_TMP0, L_irr, 0, 5, 0, 5, 0, "Z(Ij,Ab)");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_contract244(&G, &D, &Z, 0, 0, 0, -1, 0);
    dpd_buf4_close(&D);
    dpd_buf4_sort_axpy(&Z, CC_LAMBDA, qpsr, 0, 5, "New LIjAb", 1);
    dpd_buf4_init(&newL2, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
    dpd_buf4_axpy(&Z, &newL2, 1);
    dpd_buf4_close(&newL2);
    dpd_buf4_close(&Z);

    dpd_file2_close(&G);
  }
  else if(params.ref == 1) { /** ROHF **/

    dpd_file2_init(&GMI, CC_LAMBDA, L_irr, 0, 0, "GMI");
    dpd_file2_init(&Gmi, CC_LAMBDA, L_irr, 0, 0, "Gmi");

    dpd_buf4_init(&D, CC_DINTS, 0, 0, 7, 0, 7, 0, "D <ij||ab> (ij,a>b)");
    dpd_buf4_init(&X1, CC_TMP2, L_irr, 0, 7, 0, 7, 0, "X(0,7) 1");
    dpd_contract424(&D, &GMI, &X1, 1, 0, 1, -1.0, 0.0);
    dpd_buf4_init(&X2, CC_TMP2, L_irr, 0, 7, 0, 7, 0, "X(0,7) 2");
    dpd_contract244(&GMI, &D, &X2, 0, 0, 0, -1.0, 0.0);
    dpd_buf4_axpy(&X1, &X2, 1.0);
    dpd_buf4_close(&X1);
    dpd_buf4_init(&newLIJAB, CC_LAMBDA, L_irr, 0, 7, 2, 7, 0, "New LIJAB");
    dpd_buf4_axpy(&X2, &newLIJAB, 1.0);
    dpd_buf4_close(&X2);
    dpd_buf4_close(&newLIJAB);


    dpd_buf4_init(&X1, CC_TMP2, L_irr, 0, 7, 0, 7, 0, "X(0,7) 1");
    dpd_contract424(&D, &Gmi, &X1, 1, 0, 1, -1.0, 0.0);
    dpd_buf4_init(&X2, CC_TMP2, L_irr, 0, 7, 0, 7, 0, "X(0,7) 2");
    dpd_contract244(&Gmi, &D, &X2, 0, 0, 0, -1.0, 0.0);
    dpd_buf4_close(&D);
    dpd_buf4_axpy(&X1, &X2, 1.0);
    dpd_buf4_close(&X1);
    dpd_buf4_init(&newLijab, CC_LAMBDA, L_irr, 0, 7, 2, 7, 0, "New Lijab");
    dpd_buf4_axpy(&X2, &newLijab, 1.0);
    dpd_buf4_close(&X2);
    dpd_buf4_close(&newLijab);

    dpd_buf4_init(&newLIjAb, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");

    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_contract424(&D, &Gmi, &newLIjAb, 1, 0, 1, -1.0, 1.0);
    dpd_contract244(&GMI, &D, &newLIjAb, 0, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&D);

    dpd_buf4_close(&newLIjAb);

    dpd_file2_close(&Gmi);
    dpd_file2_close(&GMI);
  }
  else if(params.ref == 2) { /** UHF **/

    dpd_file2_init(&GMI, CC_LAMBDA, L_irr, 0, 0, "GMI");
    dpd_file2_init(&Gmi, CC_LAMBDA, L_irr, 2, 2, "Gmi");

    /** X(IJ,AB) = - G(M,I) <MJ||AB> **/
    dpd_buf4_init(&X1, CC_TMP2, L_irr, 0, 7, 0, 7, 0, "X(IJ,AB) C");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 7, 0, 7, 0, "D <IJ||AB> (IJ,A>B)");
    dpd_contract244(&GMI, &D, &X1, 0, 0, 0, -1, 0);
    dpd_buf4_close(&D);
    /** X(IJ,AB) --> X(JI,AB) **/
    dpd_buf4_sort(&X1, CC_TMP2, qprs, 0, 7, "X(JI,AB)");
    /** X(IJ,AB) = X(IJ,AB) - X(JI,AB) **/
    dpd_buf4_init(&X2, CC_TMP2, L_irr, 0, 7, 0, 7, 0, "X(JI,AB)");
    dpd_buf4_axpy(&X2, &X1, -1);
    dpd_buf4_close(&X2);
    /** X(IJ,AB) --> New L(IJ,AB) **/
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 0, 7, 2, 7, 0, "New LIJAB");
    dpd_buf4_axpy(&X1, &L2, 1);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&X1);

    /** X(ij,ab) = - G(m,i) <mj||ab> **/
    dpd_buf4_init(&X1, CC_TMP2, L_irr, 10, 17, 10, 17, 0, "X(ij,ab) C");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 17, 10, 17, 0, "D <ij||ab> (ij,a>b)");
    dpd_contract244(&Gmi, &D, &X1, 0, 0, 0, -1, 0);
    dpd_buf4_close(&D);
    /** X(ij,ab) --> X(ji,ab) **/
    dpd_buf4_sort(&X1, CC_TMP2, qprs, 10, 17, "X(ji,ab)");
    /** X(ij,ab) = X(ij,ab) - X(ji,ab) **/
    dpd_buf4_init(&X2, CC_TMP2, L_irr, 10, 17, 10, 17, 0, "X(ji,ab)");
    dpd_buf4_axpy(&X2, &X1, -1);
    dpd_buf4_close(&X2);
    /** X(ij,ab) --> New L(ij,ab) **/
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 10, 17, 12, 17, 0, "New Lijab");
    dpd_buf4_axpy(&X1, &L2, 1);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&X1);


    /* New L(Ij,Ab) <-- - <Im|Ab> G(m,j) - G(M,I) <Mj|Ab> **/
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "New LIjAb");
    dpd_buf4_init(&D, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    dpd_contract424(&D, &Gmi, &L2, 1, 0, 1, -1, 1);
    dpd_contract244(&GMI, &D, &L2, 0, 0, 0, -1, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&L2);

    dpd_file2_close(&Gmi);
    dpd_file2_close(&GMI);
  }

}

}} // namespace psi::cclambda
