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
    \ingroup CCHBAR
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <math.h>
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cchbar {

/* Wamef_build(): Computes all contributions to the Wamef HBAR matrix
** elements, whose spin-orbital definition is:
**
** Wamef = <am||ef> - t_n^a <nm||ef>
**
** (cf. Gauss and Stanton, JCP 103, 3561-3577 (1995).)
**
** The storage and naming convention for each spin case are
** as follows:
**
** Spin Case     Storage      Name
** ----------    ---------    --------
** WAMEF         (MA,E>F)      "WAMEF"
** Wamef         (ma,e>f)      "Wamef"
** WAmEf         (mA,Ef)       "WAmEf"
** WaMeF         (Ma,eF)       "WaMeF"
** -----------------------------------
** TDC, June 2002
**
** RHF Cases:  Note that only the WAmEf spin case is required, and
** we store it AS WRITTEN, (Am,Ef).
** TDC, March 2004
**
** For all spin cases, we now use only the following
** WAMEF         (AM,E>F)      "WAMEF"
** Wamef         (am,e>f)      "Wamef"
** WAmEf         (Am,Ef)       "WAmEf"
** WaMeF         (aM,eF)       "WaMeF"
** RAK, April 2004
**
** For CC3, these are computed by ccenergy in file CC_HET1
** RAK, July 2006
*/

void Wamef_build(void) {
  dpdbuf4 Wamef, WAMEF, WAmEf, WaMeF, W;
  dpdbuf4 F, D;
  dpdfile2 tia, tIA;
  int h, Ga, Gn, Gm, A, a, row, nrows, ncols;

  if(params.ref == 0) {

    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    global_dpd_->buf4_sort(&F, PSIF_CC_HBAR, qpsr, 11, 5, "WAmEf");
    global_dpd_->buf4_close(&F);

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_mat_init(&tIA);
    global_dpd_->file2_mat_rd(&tIA);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");

    /* dpd_contract244(&tIA, &D, &W, 0, 0, 0, -1, 1); */
    /** OOC code below added 05/04/05, -TDC **/
    for(h=0; h < moinfo.nirreps; h++) { /* h = Gam = Gnm = Gef */

      global_dpd_->buf4_mat_irrep_init(&D, h);
      global_dpd_->buf4_mat_irrep_rd(&D, h);

      row = 0;
      for(Ga=0; Ga < moinfo.nirreps; Ga++) {
	Gm = Ga ^ h;
	Gn = Ga; /* T1 is totally symmetric */

	W.matrix[h] = global_dpd_->dpd_block_matrix(moinfo.occpi[Gm], W.params->coltot[h]);

	for(A=0; A < moinfo.virtpi[Ga]; A++) {
	  a = moinfo.vir_off[Ga] + A;

	  global_dpd_->buf4_mat_irrep_rd_block(&W, h, W.row_offset[h][a], moinfo.occpi[Gm]);

	  nrows = moinfo.occpi[Gn];
	  ncols = moinfo.occpi[Gm] * W.params->coltot[h];

	  if(nrows && ncols)
	    C_DGEMV('t',nrows,ncols,-1.0,&(D.matrix[h][row][0]),ncols,&(tIA.matrix[Gn][0][A]),
		    moinfo.virtpi[Ga],1.0, W.matrix[h][0],1);


	  global_dpd_->buf4_mat_irrep_wrt_block(&W, h, W.row_offset[h][a], moinfo.occpi[Gm]);
	}

	row += moinfo.occpi[Gn] * moinfo.occpi[Gm];

	global_dpd_->free_dpd_block(W.matrix[h], moinfo.occpi[Gm], W.params->coltot[h]);
      }

      global_dpd_->buf4_mat_irrep_close(&D, h);
    }

    global_dpd_->buf4_close(&D);
    global_dpd_->file2_mat_close(&tIA);
    global_dpd_->file2_close(&tIA);
    global_dpd_->buf4_close(&W);

  }
  else if(params.ref == 1) { /** ROHF **/

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");

    /* <AM||EF> --> W(AM,E>F) */
    /* <am||ef> --> W(am,e>f) */
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 11, 7, 11, 5, 1, "F <ai|bc>");
    global_dpd_->buf4_copy(&F, PSIF_CC_HBAR, "WAMEF");
    global_dpd_->buf4_copy(&F, PSIF_CC_HBAR, "Wamef");
    global_dpd_->buf4_close(&F);

    /* T(N,A) <NM||EF> --> W(AM,E>F) */
    global_dpd_->buf4_init(&WAMEF, PSIF_CC_HBAR, 0, 11, 7, 11, 7, 0, "WAMEF");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 7, 0, 7, 0, "D <ij||ab> (ij,a>b)");
    global_dpd_->contract244(&tIA, &D, &WAMEF, 0, 0, 0, -1, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&WAMEF);

    /* T(n,a) <nm||ef> --> W(am,e>f) */
    global_dpd_->buf4_init(&Wamef, PSIF_CC_HBAR, 0, 11, 7, 11, 7, 0, "Wamef");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 7, 0, 7, 0, "D <ij||ab> (ij,a>b)");
    global_dpd_->contract244(&tia, &D, &Wamef, 0, 0, 0, -1, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&Wamef);

    /* <Am|Ef> --> W(Am,Ef) */
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 11, 5, 11, 5, 0, "F <ai|bc>");
    global_dpd_->buf4_copy(&F, PSIF_CC_HBAR, "WAmEf");
    global_dpd_->buf4_close(&F);

    /* <aM|eF> --> W(aM,eF) */
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 11, 5, 11, 5, 0, "F <ai|bc>");
    global_dpd_->buf4_copy(&F, PSIF_CC_HBAR, "WaMeF");
    global_dpd_->buf4_close(&F);

    /* T(N,A) <Nm|Ef> --> W(Am,Ef) */
    global_dpd_->buf4_init(&WAmEf, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->contract244(&tIA, &D, &WAmEf, 0, 0, 0, -1, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&WAmEf);

    /* T(n,a) <nM|eF> --> W(aM,eF) */
    global_dpd_->buf4_init(&WaMeF, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WaMeF");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->contract244(&tia, &D, &WaMeF, 0, 0, 0, -1, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&WaMeF);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);

  } /** ROHF **/
  else if(params.ref == 2) { /** UHF **/

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");

    /* <AM||EF> --> W(AM,E>F) */
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 21, 7, 21, 5, 1, "F <AI|BC>");
    global_dpd_->buf4_copy(&F, PSIF_CC_HBAR, "WAMEF");
    global_dpd_->buf4_close(&F);

    /* <am||ef> --> W(am,e>f) */
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 31, 17, 31, 15, 1, "F <ai|bc>");
    global_dpd_->buf4_copy(&F, PSIF_CC_HBAR, "Wamef");
    global_dpd_->buf4_close(&F);

    /* T(N,A) <NM||EF> --> W(AM,E>F) */
    global_dpd_->buf4_init(&WAMEF, PSIF_CC_HBAR, 0, 21, 7, 21, 7, 0, "WAMEF");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 7, 0, 7, 0, "D <IJ||AB> (IJ,A>B)");
    global_dpd_->contract244(&tIA, &D, &WAMEF, 0, 0, 0, -1, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&WAMEF);

    /* T(n,a) <nm||ef> --> W(am,e>f) */
    global_dpd_->buf4_init(&Wamef, PSIF_CC_HBAR, 0, 31, 17, 31, 17, 0, "Wamef");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 17, 10, 17, 0, "D <ij||ab> (ij,a>b)");
    global_dpd_->contract244(&tia, &D, &Wamef, 0, 0, 0, -1, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&Wamef);

    /* <Am|Ef> --> W(Am,Ef) */
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 26, 28, 26, 28, 0, "F <Ai|Bc>");
    global_dpd_->buf4_copy(&F, PSIF_CC_HBAR, "WAmEf");
    global_dpd_->buf4_close(&F);

    /* <aM|eF> --> W(aM,eF) */
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 25, 29, 25, 29, 0, "F <aI|bC>");
    global_dpd_->buf4_copy(&F, PSIF_CC_HBAR, "WaMeF");
    global_dpd_->buf4_close(&F);

    /* T(N,A) <Nm|Ef> --> W(Am,Ef) */
    global_dpd_->buf4_init(&WAmEf, PSIF_CC_HBAR, 0, 26, 28, 26, 28, 0, "WAmEf");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    global_dpd_->contract244(&tIA, &D, &WAmEf, 0, 0, 0, -1, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&WAmEf);

    /* T(n,a) <nM|eF> --> W(aM,eF) */
    global_dpd_->buf4_init(&WaMeF, PSIF_CC_HBAR, 0, 25, 29, 25, 29, 0, "WaMeF");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");
    global_dpd_->contract244(&tia, &D, &WaMeF, 0, 0, 0, -1, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&WaMeF);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);
  } /** UHF **/

  return;
}

}} // namespace psi::cchbar
