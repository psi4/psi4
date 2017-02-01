/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*! \file
    \ingroup CCLAMBDA
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstring>
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
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

    global_dpd_->file2_init(&LIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "LIA");

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, L_irr, 0, 5, 0, 5, 0, "ZIjAb");
    global_dpd_->buf4_scm(&Z, 0);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5,  0, "WAmEf");
    /*       dpd_contract244(&LIA, &W, &Z, 1, 2, 1, 1, 0); */
    /* Out-of-core contract244 */
    GW = W.file.my_irrep;
    GZ = Z.file.my_irrep;
    GL1 = LIA.my_irrep;

    global_dpd_->file2_mat_init(&LIA);
    global_dpd_->file2_mat_rd(&LIA);

    for(Gej=0; Gej < moinfo.nirreps; Gej++) {
      Gab = Gej^GW;
      Gij = Gab^GZ;

      global_dpd_->buf4_mat_irrep_init(&Z, Gij);

      for(Ge=0; Ge < moinfo.nirreps; Ge++) {
	Gi = Ge^GL1;
	Gj = GZ^Gab^Gi;

	num_j = Z.params->qpi[Gj];
	num_i = LIA.params->rowtot[Gi];
	num_e = LIA.params->coltot[Ge];

	global_dpd_->buf4_mat_irrep_init_block(&W, Gej, num_j);

	for(e=0; e < num_e; e++) {

	  E = W.params->poff[Ge] + e;
	  global_dpd_->buf4_mat_irrep_rd_block(&W, Gej, W.row_offset[Gej][E], num_j);

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
	global_dpd_->buf4_mat_irrep_close_block(&W, Gej, num_j);
      }
      global_dpd_->buf4_mat_irrep_wrt(&Z, Gij);
      global_dpd_->buf4_mat_irrep_close(&Z, Gij);
    }
    global_dpd_->file2_mat_close(&LIA);

    /* End out-of-core contract244 */
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
    global_dpd_->buf4_axpy(&Z, &L2, 1);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_LAMBDA, qpsr, 0, 5, "New LIjAb", 1);
    global_dpd_->buf4_close(&Z);
    global_dpd_->file2_close(&LIA);
  }
  else if(params.ref == 1) { /** ROHF **/

    global_dpd_->file2_init(&LIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "LIA");
    global_dpd_->file2_init(&Lia, PSIF_CC_LAMBDA, L_irr, 0, 1, "Lia");

    /** Z(IJ,AB) = L(I,E) W(EJ,AB) **/
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, L_irr, 0, 7, 0, 7, 0, "Z(IJ,A>B)");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 7, 11, 7, 0, "WAMEF");
    global_dpd_->contract244(&LIA, &W, &Z1, 1, 0, 0, 1, 0);
    global_dpd_->buf4_close(&W);
    /** Z(IJ,AB) --> Z(JI,AB) **/
    global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP1, qprs, 0, 7, "Z(JI,A>B)");
    /** Z(IJ,AB) = Z(IJ,AB) - Z(JI,AB) **/
    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP1, L_irr, 0, 7, 0, 7, 0, "Z(JI,A>B)");
    global_dpd_->buf4_axpy(&Z2, &Z1, -1);
    global_dpd_->buf4_close(&Z2);
    /** Z(IJ,AB) --> New L(IJ,AB) **/
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 7, 2, 7, 0, "New LIJAB");
    global_dpd_->buf4_axpy(&Z1, &L2, 1);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&Z1);

    /** Z(ij,ab) = L(i,e) W(ej,ab) **/
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, L_irr, 0, 7, 0, 7, 0, "Z(ij,a>b)");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 7, 11, 7, 0, "Wamef");
    global_dpd_->contract244(&Lia, &W, &Z1, 1, 0, 0, 1, 0);
    global_dpd_->buf4_close(&W);
    /** Z(ij,ab) --> Z(ji,ab) **/
    global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP1, qprs, 0, 7, "Z(ji,a>b)");
    /** Z(ij,ab) = Z(ij,ab) - Z(ji,ab) **/
    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP1, L_irr, 0, 7, 0, 7, 0, "Z(ji,a>b)");
    global_dpd_->buf4_axpy(&Z2, &Z1, -1);
    global_dpd_->buf4_close(&Z2);
    /** Z(ij,ab) --> New L(ij,ab) **/
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 7, 2, 7, 0, "New Lijab");
    global_dpd_->buf4_axpy(&Z1, &L2, 1);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&Z1);

    /** New L(Ij,Ab) <-- L(I,E) W(Ej,Ab) **/
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
    global_dpd_->contract244(&LIA, &W, &L2, 1, 0, 0, 1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&L2);

    /** Z(jI,bA) = -L(j,e) W(eI,bA) **/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, L_irr, 0, 5, 0, 5, 0, "Z(jI,bA)");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WaMeF");
    global_dpd_->contract244(&Lia, &W, &Z, 1, 0, 0, 1, 0);
    global_dpd_->buf4_close(&W);
    /** Z(jI,bA) --> New L(Ij,Ab) **/
    global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_LAMBDA, qpsr, 0, 5, "New LIjAb", 1);
    global_dpd_->buf4_close(&Z);

    global_dpd_->file2_close(&Lia);
    global_dpd_->file2_close(&LIA);

  }
  else if(params.ref == 2) { /** UHF **/

    global_dpd_->file2_init(&LIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "LIA");
    global_dpd_->file2_init(&Lia, PSIF_CC_LAMBDA, L_irr, 2, 3, "Lia");

    /** Z(IJ,AB) = L(I,E) W(EJ,AB) **/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, L_irr, 0, 7, 0, 7, 0, "Z(IJ,AB)");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 21, 7, 21, 7, 0, "WAMEF");
    global_dpd_->contract244(&LIA, &W, &Z, 1, 0, 0, 1, 0);
    global_dpd_->buf4_close(&W);
    /** Z(IJ,AB) --> Z(JI,AB) **/
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP1, qprs, 0, 7, "Z(JI,AB)");
    global_dpd_->buf4_close(&Z);
    /** Z(IJ,AB) = Z(IJ,AB) - Z(JI,AB) **/
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, L_irr, 0, 7, 0, 7, 0, "Z(IJ,AB)");
    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP1, L_irr, 0, 7, 0, 7, 0, "Z(JI,AB)");
    global_dpd_->buf4_axpy(&Z2, &Z1, -1);
    global_dpd_->buf4_close(&Z2);
    /** Z(IJ,AB) --> New L(IJ,AB) **/
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 7, 2, 7, 0, "New LIJAB");
    global_dpd_->buf4_axpy(&Z1, &L2, 1);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&Z1);

    /** Z(ij,ab) = L(i,e) W(ej,ab) **/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, L_irr, 10, 17, 10, 17, 0, "Z(ij,ab)");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 31, 17, 31, 17, 0, "Wamef");
    global_dpd_->contract244(&Lia, &W, &Z, 1, 0, 0, 1, 0);
    global_dpd_->buf4_close(&W);
    /** Z(ij,ab) --> Z(ji,ab) **/
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP1, qprs, 10, 17, "Z(ji,ab)");
    global_dpd_->buf4_close(&Z);
    /** Z(ij,ab) = Z(ij,ab) - Z(ji,ab) **/
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, L_irr, 10, 17, 10, 17, 0, "Z(ij,ab)");
    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP1, L_irr, 10, 17, 10, 17, 0, "Z(ji,ab)");
    global_dpd_->buf4_axpy(&Z2, &Z1, -1);
    global_dpd_->buf4_close(&Z2);
    /** Z(ij,ab) --> New L(ij,ab) **/
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 10, 17, 12, 17, 0, "New Lijab");
    global_dpd_->buf4_axpy(&Z1, &L2, 1);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&Z1);


    /** New L(Ij,Ab) <-- L(I,E) W(Ej,Ab) **/
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "New LIjAb");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 26, 28, 26, 28, 0, "WAmEf");
    global_dpd_->contract244(&LIA, &W, &L2, 1, 0, 0, 1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&L2);

    /** Z(jI,bA) = -L(j,e) W(eI,bA) **/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, L_irr, 23, 29, 23, 29, 0, "Z(jI,bA)");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 25, 29, 25, 29, 0, "WaMeF");
    global_dpd_->contract244(&Lia, &W, &Z, 1, 0, 0, 1, 0);
    global_dpd_->buf4_close(&W);
    /** Z(jI,bA) --> New L(Ij,Ab) **/
    global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_LAMBDA, qpsr, 22, 28, "New LIjAb", 1);
    global_dpd_->buf4_close(&Z);

    global_dpd_->file2_close(&Lia);
    global_dpd_->file2_close(&LIA);
  }
}

}} // namespace psi::cclambda
