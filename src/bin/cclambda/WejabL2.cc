/*! \file
    \ingroup CCLAMBDA
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstring>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cclambda {

/* WejabL2(): Computes the contribution of the Wamef HBAR matrix
** elements to the Lambda double de-excitation amplitude equations.
** These contributions are given in spin orbitals as:
**
** L_ij^ab = P(ij) L_i^e Wejab
**
** where Wejab = Wamef and is defined as:
**
** Wamef = <am||ef> - t_n^a <nm||ef>
**
** [cf. Gauss and Stanton, JCP 103, 3561-3577 (1995).]
**
** By spin case of L_ij^ab, these contributions are computed as:
**
** L(IJ,AB) = L(I,E) W(EJ,AB) - L(J,E) W(EI,AB) (only one unique contraction)
** L(ij,ab) = L(i,e) W(ej,ab) - L(j,e) W(ei,ab) (only one unique contraction)
** L(Ij,Ab) = L(I,E) W(Ej,Ab) + L(j,e) W(eI,bA)
**
** TDC, July 2002
**
** NB: The ROHF case needs to be re-written to use (AM,EF) ordering
** for the Wamef matrix elements, as I've done for the UHF case.
*/

void WejabL2(int L_irr)
{
  int GW, GL1, GZ, Gej, Gab, Gi, Ge, Gij, Gj, Ga;
  int e, E, i, I, num_j, num_i, num_e, nlinks;
  dpdbuf4 W, L2;
  dpdfile2 LIA, Lia;
  dpdbuf4 Z, Z1, Z2;
  
  /* RHS += P(ij) Lie * Wejab */
  if(params.ref == 0) { /** RHF **/

    dpd_file2_init(&LIA, CC_LAMBDA, L_irr, 0, 1, "LIA");

    dpd_buf4_init(&Z, CC_TMP0, L_irr, 0, 5, 0, 5, 0, "ZIjAb");
    dpd_buf4_scm(&Z, 0);
    dpd_buf4_init(&W, CC_HBAR, 0, 11, 5, 11, 5,  0, "WAmEf");
    /*       dpd_contract244(&LIA, &W, &Z, 1, 2, 1, 1, 0); */
    /* Out-of-core contract244 */
    GW = W.file.my_irrep;
    GZ = Z.file.my_irrep;
    GL1 = LIA.my_irrep;

    dpd_file2_mat_init(&LIA);
    dpd_file2_mat_rd(&LIA);

    for(Gej=0; Gej < moinfo.nirreps; Gej++) {
      Gab = Gej^GW;
      Gij = Gab^GZ;

      dpd_buf4_mat_irrep_init(&Z, Gij);

      for(Ge=0; Ge < moinfo.nirreps; Ge++) {
	Gi = Ge^GL1;
	Gj = GZ^Gab^Gi;

	num_j = Z.params->qpi[Gj];
	num_i = LIA.params->rowtot[Gi];
	num_e = LIA.params->coltot[Ge];

	dpd_buf4_mat_irrep_init_block(&W, Gej, num_j);

	for(e=0; e < num_e; e++) {

	  E = W.params->poff[Ge] + e;
	  dpd_buf4_mat_irrep_rd_block(&W, Gej, W.row_offset[Gej][E], num_j);

	  for(i=0; i < num_i; i++) {
	    I = Z.params->poff[Gi] + i;

	    nlinks = Z.params->coltot[Gab] * num_j;
	    if(nlinks) {
	      C_DAXPY(nlinks, LIA.matrix[Gi][i][e],
		      &(W.matrix[Gej][0][0]),1,
		      &(Z.matrix[Gij][Z.row_offset[Gij][I]][0]),1);
	    }
	  }
	}
	dpd_buf4_mat_irrep_close_block(&W, Gej, num_j);
      }
      dpd_buf4_mat_irrep_wrt(&Z, Gij);
      dpd_buf4_mat_irrep_close(&Z, Gij);
    }
    dpd_file2_mat_close(&LIA);

    /* End out-of-core contract244 */
    dpd_buf4_close(&W);

    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
    dpd_buf4_axpy(&Z, &L2, 1);
    dpd_buf4_close(&L2);
    dpd_buf4_sort_axpy(&Z, CC_LAMBDA, qpsr, 0, 5, "New LIjAb", 1);
    dpd_buf4_close(&Z);
    dpd_file2_close(&LIA);
  }
  else if(params.ref == 1) { /** ROHF **/
    
    dpd_file2_init(&LIA, CC_LAMBDA, L_irr, 0, 1, "LIA");
    dpd_file2_init(&Lia, CC_LAMBDA, L_irr, 0, 1, "Lia");

    /** Z(IJ,AB) = L(I,E) W(EJ,AB) **/
    dpd_buf4_init(&Z1, CC_TMP1, L_irr, 0, 7, 0, 7, 0, "Z(IJ,A>B)");
    dpd_buf4_init(&W, CC_HBAR, 0, 11, 7, 11, 7, 0, "WAMEF");
    dpd_contract244(&LIA, &W, &Z1, 1, 0, 0, 1, 0);
    dpd_buf4_close(&W);
    /** Z(IJ,AB) --> Z(JI,AB) **/
    dpd_buf4_sort(&Z1, CC_TMP1, qprs, 0, 7, "Z(JI,A>B)");
    /** Z(IJ,AB) = Z(IJ,AB) - Z(JI,AB) **/
    dpd_buf4_init(&Z2, CC_TMP1, L_irr, 0, 7, 0, 7, 0, "Z(JI,A>B)");
    dpd_buf4_axpy(&Z2, &Z1, -1);
    dpd_buf4_close(&Z2);
    /** Z(IJ,AB) --> New L(IJ,AB) **/
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 0, 7, 2, 7, 0, "New LIJAB");
    dpd_buf4_axpy(&Z1, &L2, 1);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&Z1);

    /** Z(ij,ab) = L(i,e) W(ej,ab) **/
    dpd_buf4_init(&Z1, CC_TMP1, L_irr, 0, 7, 0, 7, 0, "Z(ij,a>b)");
    dpd_buf4_init(&W, CC_HBAR, 0, 11, 7, 11, 7, 0, "Wamef");
    dpd_contract244(&Lia, &W, &Z1, 1, 0, 0, 1, 0);
    dpd_buf4_close(&W);
    /** Z(ij,ab) --> Z(ji,ab) **/
    dpd_buf4_sort(&Z1, CC_TMP1, qprs, 0, 7, "Z(ji,a>b)");
    /** Z(ij,ab) = Z(ij,ab) - Z(ji,ab) **/
    dpd_buf4_init(&Z2, CC_TMP1, L_irr, 0, 7, 0, 7, 0, "Z(ji,a>b)");
    dpd_buf4_axpy(&Z2, &Z1, -1);
    dpd_buf4_close(&Z2);
    /** Z(ij,ab) --> New L(ij,ab) **/
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 0, 7, 2, 7, 0, "New Lijab");
    dpd_buf4_axpy(&Z1, &L2, 1);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&Z1);

    /** New L(Ij,Ab) <-- L(I,E) W(Ej,Ab) **/
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
    dpd_buf4_init(&W, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
    dpd_contract244(&LIA, &W, &L2, 1, 0, 0, 1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&L2);

    /** Z(jI,bA) = -L(j,e) W(eI,bA) **/
    dpd_buf4_init(&Z, CC_TMP1, L_irr, 0, 5, 0, 5, 0, "Z(jI,bA)");
    dpd_buf4_init(&W, CC_HBAR, 0, 11, 5, 11, 5, 0, "WaMeF");
    dpd_contract244(&Lia, &W, &Z, 1, 0, 0, 1, 0);
    dpd_buf4_close(&W);
    /** Z(jI,bA) --> New L(Ij,Ab) **/
    dpd_buf4_sort_axpy(&Z, CC_LAMBDA, qpsr, 0, 5, "New LIjAb", 1);
    dpd_buf4_close(&Z);

    dpd_file2_close(&Lia);
    dpd_file2_close(&LIA);

  }
  else if(params.ref == 2) { /** UHF **/

    dpd_file2_init(&LIA, CC_LAMBDA, L_irr, 0, 1, "LIA");
    dpd_file2_init(&Lia, CC_LAMBDA, L_irr, 2, 3, "Lia");

    /** Z(IJ,AB) = L(I,E) W(EJ,AB) **/
    dpd_buf4_init(&Z, CC_TMP1, L_irr, 0, 7, 0, 7, 0, "Z(IJ,AB)");
    dpd_buf4_init(&W, CC_HBAR, 0, 21, 7, 21, 7, 0, "WAMEF");
    dpd_contract244(&LIA, &W, &Z, 1, 0, 0, 1, 0);
    dpd_buf4_close(&W);
    /** Z(IJ,AB) --> Z(JI,AB) **/
    dpd_buf4_sort(&Z, CC_TMP1, qprs, 0, 7, "Z(JI,AB)");
    dpd_buf4_close(&Z);
    /** Z(IJ,AB) = Z(IJ,AB) - Z(JI,AB) **/
    dpd_buf4_init(&Z1, CC_TMP1, L_irr, 0, 7, 0, 7, 0, "Z(IJ,AB)");
    dpd_buf4_init(&Z2, CC_TMP1, L_irr, 0, 7, 0, 7, 0, "Z(JI,AB)");
    dpd_buf4_axpy(&Z2, &Z1, -1);
    dpd_buf4_close(&Z2);
    /** Z(IJ,AB) --> New L(IJ,AB) **/
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 0, 7, 2, 7, 0, "New LIJAB");
    dpd_buf4_axpy(&Z1, &L2, 1);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&Z1);

    /** Z(ij,ab) = L(i,e) W(ej,ab) **/
    dpd_buf4_init(&Z, CC_TMP1, L_irr, 10, 17, 10, 17, 0, "Z(ij,ab)");
    dpd_buf4_init(&W, CC_HBAR, 0, 31, 17, 31, 17, 0, "Wamef");
    dpd_contract244(&Lia, &W, &Z, 1, 0, 0, 1, 0);
    dpd_buf4_close(&W);
    /** Z(ij,ab) --> Z(ji,ab) **/
    dpd_buf4_sort(&Z, CC_TMP1, qprs, 10, 17, "Z(ji,ab)");
    dpd_buf4_close(&Z);
    /** Z(ij,ab) = Z(ij,ab) - Z(ji,ab) **/
    dpd_buf4_init(&Z1, CC_TMP1, L_irr, 10, 17, 10, 17, 0, "Z(ij,ab)");
    dpd_buf4_init(&Z2, CC_TMP1, L_irr, 10, 17, 10, 17, 0, "Z(ji,ab)");
    dpd_buf4_axpy(&Z2, &Z1, -1);
    dpd_buf4_close(&Z2);
    /** Z(ij,ab) --> New L(ij,ab) **/
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 10, 17, 12, 17, 0, "New Lijab");
    dpd_buf4_axpy(&Z1, &L2, 1);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&Z1);


    /** New L(Ij,Ab) <-- L(I,E) W(Ej,Ab) **/
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "New LIjAb");
    dpd_buf4_init(&W, CC_HBAR, 0, 26, 28, 26, 28, 0, "WAmEf");
    dpd_contract244(&LIA, &W, &L2, 1, 0, 0, 1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&L2);

    /** Z(jI,bA) = -L(j,e) W(eI,bA) **/
    dpd_buf4_init(&Z, CC_TMP1, L_irr, 23, 29, 23, 29, 0, "Z(jI,bA)");
    dpd_buf4_init(&W, CC_HBAR, 0, 25, 29, 25, 29, 0, "WaMeF");
    dpd_contract244(&Lia, &W, &Z, 1, 0, 0, 1, 0);
    dpd_buf4_close(&W);
    /** Z(jI,bA) --> New L(Ij,Ab) **/
    dpd_buf4_sort_axpy(&Z, CC_LAMBDA, qpsr, 22, 28, "New LIjAb", 1);
    dpd_buf4_close(&Z);

    dpd_file2_close(&Lia);
    dpd_file2_close(&LIA);
  }
}

}} // namespace psi::cclambda
