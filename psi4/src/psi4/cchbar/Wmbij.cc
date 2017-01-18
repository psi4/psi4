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
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cchbar {

/** Wmbij_build(): Constructs the Wmbij HBAR intermediate, defined in
 ** spin orbitals as:
 **
 ** Wmbij = <mb||ij> - Fme t_ij^be - t_n^b Wmnij + 1/2 <mb||ef> tau_ij^ef
 **    + P(ij) <mn||ie> t_jn^be + P(ij) t_i^e { <mb||ej> - t_nj^bf <mn||ef> }
 **
 ** [cf. Gauss and Stanton, JCP 103, 3561-3577 (1995)]
 **
 ** For RHF orbitals, there is only one unique spin case: WMbIj.
 **
 ** TDC, March 2004
 */

void Wmbij_build(void)
{
  dpdfile2 Fme, T1;
  dpdbuf4 W, E, T2, Wmnij, I, Tau, Z, Z1, Z2, C, D;

  if(params.ref == 0) { /** RHF **/

    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    global_dpd_->buf4_sort(&E, PSIF_CC_HBAR, rspq, 10, 0, "WMbIj");
    global_dpd_->buf4_close(&E);

  }
  else if(params.ref == 1) { /** ROHF **/

    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 2, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
    /** <MB||IJ> **/
    global_dpd_->buf4_sort(&E, PSIF_CC_HBAR, rspq, 10, 2, "WMBIJ");
    /** <mb||ij> **/
    global_dpd_->buf4_sort(&E, PSIF_CC_HBAR, rspq, 10, 2, "Wmbij");
    global_dpd_->buf4_close(&E);

    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    /** <Mb|Ij> **/
    global_dpd_->buf4_sort(&E, PSIF_CC_HBAR, rspq, 10, 0, "WMbIj");
    /** <mB|iJ> **/
    global_dpd_->buf4_sort(&E, PSIF_CC_HBAR, rspq, 10, 0, "WmBiJ");
    global_dpd_->buf4_close(&E);

  }
  else if(params.ref == 2) { /** UHF **/

    /** <MB||IJ> **/
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 2, 20, 2, 20, 0, "E <IJ||KA> (I>J,KA)");
    global_dpd_->buf4_sort(&E, PSIF_CC_HBAR, rspq, 20, 2, "WMBIJ");
    global_dpd_->buf4_close(&E);

    /** <mb||ij> **/
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 12, 30, 12, 30, 0, "E <ij||ka> (i>j,ka)");
    global_dpd_->buf4_sort(&E, PSIF_CC_HBAR, rspq, 30, 12, "Wmbij");
    global_dpd_->buf4_close(&E);

    /** <Mb|Ij> **/
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
    global_dpd_->buf4_sort(&E, PSIF_CC_HBAR, rspq, 24, 22, "WMbIj");
    global_dpd_->buf4_close(&E);

    /** <mB|iJ> **/
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 23, 27, 23, 27, 0, "E <iJ|kA>");
    global_dpd_->buf4_sort(&E, PSIF_CC_HBAR, rspq, 27, 23, "WmBiJ");
    global_dpd_->buf4_close(&E);

  }

  if(params.ref == 0) { /** RHF **/

    /** F_ME t_Ij^Eb --> W(Mb,Ij) **/
    global_dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 0, 1, "FME");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 0, 10, 0, 0, "WMbIj");
    global_dpd_->contract244(&Fme, &T2, &W, 1, 2, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_close(&Fme);

  }
  else if(params.ref ==1) { /** ROHF **/

    /** F_ME t_IJ^EB --> W(MB,IJ) **/
    global_dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 0, 1, "FME");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 2, 10, 2, 0, "WMBIJ");
    global_dpd_->contract244(&Fme, &T2, &W, 1, 2, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_close(&Fme);

    /** F_me t_ij^eb --> W(mb,ij) **/
    global_dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 0, 1, "Fme");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tijab");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 2, 10, 2, 0, "Wmbij");
    global_dpd_->contract244(&Fme, &T2, &W, 1, 2, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_close(&Fme);

    /** F_ME t_Ij^Eb --> W(Mb,Ij) **/
    global_dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 0, 1, "FME");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 0, 10, 0, 0, "WMbIj");
    global_dpd_->contract244(&Fme, &T2, &W, 1, 2, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_close(&Fme);

    /** F_me t_iJ^eB --> W(mB,iJ) **/
    global_dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 0, 1, "Fme");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 0, 10, 0, 0, "WmBiJ");
    global_dpd_->contract244(&Fme, &T2, &W, 1, 2, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_close(&Fme);

  }
  else if(params.ref == 2) { /** UHF **/

    /** F_ME t_IJ^EB --> W(MB,IJ) **/
    global_dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 0, 1, "FME");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 20, 2, 20, 2, 0, "WMBIJ");
    global_dpd_->contract244(&Fme, &T2, &W, 1, 2, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_close(&Fme);

    /** F_me t_ij^eb --> W(mb,ij) **/
    global_dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 2, 3, "Fme");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 12, 15, 12, 17, 0, "tijab");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 30, 12, 30, 12, 0, "Wmbij");
    global_dpd_->contract244(&Fme, &T2, &W, 1, 2, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_close(&Fme);

    /** F_ME t_Ij^Eb --> W(Mb,Ij) **/
    global_dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 0, 1, "FME");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 24, 22, 24, 22, 0, "WMbIj");
    global_dpd_->contract244(&Fme, &T2, &W, 1, 2, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_close(&Fme);

    /** F_me t_iJ^eB --> W(mB,iJ) **/
    global_dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 2, 3, "Fme");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 27, 23, 27, 23, 0, "WmBiJ");
    global_dpd_->contract244(&Fme, &T2, &W, 1, 2, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_close(&Fme);

  }

  if(params.ref == 0) { /** RHF **/

    /** - t_n^b W_MnIj **/
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->buf4_init(&Wmnij, PSIF_CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 0, 10, 0, 0, "WMbIj");
    global_dpd_->contract424(&Wmnij, &T1, &W, 1, 0, 1, -1.0, 1.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Wmnij);
    global_dpd_->file2_close(&T1);

  }
  else if(params.ref == 1) { /** ROHF **/

    /** - t_N^B W_MNIJ --> W(MB,IJ) **/
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->buf4_init(&Wmnij, PSIF_CC_HBAR, 0, 0, 2, 2, 2, 0, "WMNIJ");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 2, 10, 2, 0, "WMBIJ");
    global_dpd_->contract424(&Wmnij, &T1, &W, 1, 0, 1, -1.0, 1.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Wmnij);
    global_dpd_->file2_close(&T1);

    /** - t_n^b W_mnij **/
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
    global_dpd_->buf4_init(&Wmnij, PSIF_CC_HBAR, 0, 0, 2, 2, 2, 0, "Wmnij");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 2, 10, 2, 0, "Wmbij");
    global_dpd_->contract424(&Wmnij, &T1, &W, 1, 0, 1, -1.0, 1.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Wmnij);
    global_dpd_->file2_close(&T1);

    /** - t_n^b W_MnIj **/
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
    global_dpd_->buf4_init(&Wmnij, PSIF_CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 0, 10, 0, 0, "WMbIj");
    global_dpd_->contract424(&Wmnij, &T1, &W, 1, 0, 1, -1.0, 1.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Wmnij);
    global_dpd_->file2_close(&T1);

    /** - t_N^B W_mNiJ **/
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->buf4_init(&Wmnij, PSIF_CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj");
    global_dpd_->buf4_sort(&Wmnij, PSIF_CC_TMP0, qprs, 0, 0, "WnMIj");
    global_dpd_->buf4_close(&Wmnij);
    global_dpd_->buf4_init(&Wmnij, PSIF_CC_TMP0, 0, 0, 0, 0, 0, 0, "WnMIj");
    global_dpd_->buf4_sort(&Wmnij, PSIF_CC_TMP1, pqsr, 0, 0, "WnMjI");
    global_dpd_->buf4_close(&Wmnij);
    global_dpd_->buf4_init(&Wmnij, PSIF_CC_TMP1, 0, 0, 0, 0, 0, 0, "WnMjI");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 0, 10, 0, 0, "WmBiJ");
    global_dpd_->contract424(&Wmnij, &T1, &W, 1, 0, 1, -1.0, 1.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Wmnij);
    global_dpd_->file2_close(&T1);
  }
  else if(params.ref == 2) { /** UHF **/

    /** - t_N^B W_MNIJ --> W(MB,IJ) **/
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->buf4_init(&Wmnij, PSIF_CC_HBAR, 0, 0, 2, 2, 2, 0, "WMNIJ");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 20, 2, 20, 2, 0, "WMBIJ");
    global_dpd_->contract424(&Wmnij, &T1, &W, 1, 0, 1, -1.0, 1.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Wmnij);
    global_dpd_->file2_close(&T1);

    /** - t_n^b W_mnij **/
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
    global_dpd_->buf4_init(&Wmnij, PSIF_CC_HBAR, 0, 10, 12, 12, 12, 0, "Wmnij");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 30, 12, 30, 12, 0, "Wmbij");
    global_dpd_->contract424(&Wmnij, &T1, &W, 1, 0, 1, -1.0, 1.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Wmnij);
    global_dpd_->file2_close(&T1);

    /** - t_n^b W_MnIj **/
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
    global_dpd_->buf4_init(&Wmnij, PSIF_CC_HBAR, 0, 22, 22, 22, 22, 0, "WMnIj");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 24, 22, 24, 22, 0, "WMbIj");
    global_dpd_->contract424(&Wmnij, &T1, &W, 1, 0, 1, -1.0, 1.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Wmnij);
    global_dpd_->file2_close(&T1);

    /** - t_N^B W_mNiJ **/
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->buf4_init(&Wmnij, PSIF_CC_HBAR, 0, 22, 22, 22, 22, 0, "WMnIj");
    global_dpd_->buf4_sort(&Wmnij, PSIF_CC_TMP0, qpsr, 23, 23, "WmNiJ");
    global_dpd_->buf4_close(&Wmnij);
    global_dpd_->buf4_init(&Wmnij, PSIF_CC_TMP0, 0, 23, 23, 23, 23, 0, "WmNiJ");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 27, 23, 27, 23, 0, "WmBiJ");
    global_dpd_->contract424(&Wmnij, &T1, &W, 1, 0, 1, -1.0, 1.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Wmnij);
    global_dpd_->file2_close(&T1);

  }

  if(params.ref == 0) { /** RHF **/

    /** <Mb|Ef> tau_Ij^Ef **/
    global_dpd_->buf4_init(&I, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    global_dpd_->buf4_init(&Tau, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 0, 10, 0, 0, "WMbIj");
    global_dpd_->contract444(&I, &Tau, &W, 0, 0, 1.0, 1.0); /* should run OCC, if needed */
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Tau);
    global_dpd_->buf4_close(&I);

  }
  else if(params.ref == 1) { /** ROHF **/

    /** <MB||EF> tau_IJ^EF **/
    global_dpd_->buf4_init(&I, PSIF_CC_FINTS, 0, 10, 7, 10, 5, 1, "F <ia|bc>");
    global_dpd_->buf4_init(&Tau, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 2, 10, 2, 0, "WMBIJ");
    global_dpd_->contract444(&I, &Tau, &W, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Tau);
    global_dpd_->buf4_close(&I);

    /* <mb||ef> tau_ij^ef **/
    global_dpd_->buf4_init(&I, PSIF_CC_FINTS, 0, 10, 7, 10, 5, 1, "F <ia|bc>");
    global_dpd_->buf4_init(&Tau, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauijab");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 2, 10, 2, 0, "Wmbij");
    global_dpd_->contract444(&I, &Tau, &W, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Tau);
    global_dpd_->buf4_close(&I);

    /** <Mb|Ef> tau_Ij^Ef **/
    global_dpd_->buf4_init(&I, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    global_dpd_->buf4_init(&Tau, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 0, 10, 0, 0, "WMbIj");
    global_dpd_->contract444(&I, &Tau, &W, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Tau);
    global_dpd_->buf4_close(&I);

    /** <mB|eF> tau_iJ^eF **/
    global_dpd_->buf4_init(&I, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    global_dpd_->buf4_init(&Tau, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauiJaB");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 0, 10, 0, 0, "WmBiJ");
    global_dpd_->contract444(&I, &Tau, &W, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Tau);
    global_dpd_->buf4_close(&I);
  }
  else if(params.ref == 2) { /** UHF **/

    /** <MB||EF> tau_IJ^EF **/
    global_dpd_->buf4_init(&I, PSIF_CC_FINTS, 0, 20, 7, 20, 5, 1, "F <IA|BC>");
    global_dpd_->buf4_init(&Tau, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 20, 2, 20, 2, 0, "WMBIJ");
    global_dpd_->contract444(&I, &Tau, &W, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&Tau);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&W);

    /* <mb||ef> tau_ij^ef **/
    global_dpd_->buf4_init(&I, PSIF_CC_FINTS, 0, 30, 17, 30, 15, 1, "F <ia|bc>");
    global_dpd_->buf4_init(&Tau, PSIF_CC_TAMPS, 0, 12, 17, 12, 17, 0, "tauijab");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 30, 12, 30, 12, 0, "Wmbij");
    global_dpd_->contract444(&I, &Tau, &W, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Tau);
    global_dpd_->buf4_close(&I);

    /** <Mb|Ef> tau_Ij^Ef **/
    global_dpd_->buf4_init(&I, PSIF_CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
    global_dpd_->buf4_init(&Tau, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tauIjAb");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 24, 22, 24, 22, 0, "WMbIj");
    global_dpd_->contract444(&I, &Tau, &W, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Tau);
    global_dpd_->buf4_close(&I);

    /** <mB|eF> tau_iJ^eF **/
    global_dpd_->buf4_init(&I, PSIF_CC_FINTS, 0, 27, 29, 27, 29, 0, "F <iA|bC>");
    global_dpd_->buf4_init(&Tau, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tauiJaB");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 27, 23, 27, 23, 0, "WmBiJ");
    global_dpd_->contract444(&I, &Tau, &W, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Tau);
    global_dpd_->buf4_close(&I);

  }

  /* Sort <ij||ka> integrals for the E*T2 contributions */
  if(params.ref == 0) { /** RHF **/

    global_dpd_->buf4_init(&I, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    global_dpd_->buf4_sort(&I, PSIF_CC_TMP0, prqs, 0, 10, "<Mn|Ie> (MI,ne)");
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E 2<ai|jk> - <ai|kj>");
    global_dpd_->buf4_sort(&I, PSIF_CC_TMP0, sqrp, 0, 10, "2 <Mn|Ie> - <Nm|Ie> (MI,ne)");
    global_dpd_->buf4_close(&I);

  }
  else if(params.ref == 1) { /** ROHF **/

    global_dpd_->buf4_init(&I, PSIF_CC_EINTS, 0, 0, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
    global_dpd_->buf4_sort(&I, PSIF_CC_TMP0, prqs, 0, 10, "I(MI,NE)");
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    global_dpd_->buf4_sort(&I, PSIF_CC_TMP1, prqs, 0, 10, "I(MI,NE)");
    global_dpd_->buf4_close(&I);
  }
  else if(params.ref == 2) { /** UHF **/

    global_dpd_->buf4_init(&I, PSIF_CC_EINTS, 0, 0, 20, 2, 20, 0, "E <IJ||KA> (I>J,KA)");
    global_dpd_->buf4_sort(&I, PSIF_CC_TMP0, prqs, 0, 20, "<MN||IE> (MI,NE)");
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_CC_EINTS, 0, 10, 30, 12, 30, 0, "E <ij||ka> (i>j,ka)");
    global_dpd_->buf4_sort(&I, PSIF_CC_TMP0, prqs, 10, 30, "<mn||ie> (mi,ne)");
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
    global_dpd_->buf4_sort(&I, PSIF_CC_TMP0, prqs, 0, 30, "<Mn|Ie> (MI,ne)");
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_CC_EINTS, 0, 23, 27, 23, 27, 0, "E <iJ|kA>");
    global_dpd_->buf4_sort(&I, PSIF_CC_TMP0, prqs, 10, 20, "<mN|iE> (mi,NE)");
    global_dpd_->buf4_close(&I);

  }

  if(params.ref == 0) { /** RHF **/

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 0, 10, 0, 10, 0, "Z(MI,jb)");

    global_dpd_->buf4_init(&I, PSIF_CC_TMP0, 0, 0, 10, 0, 10, 0, "2 <Mn|Ie> - <Nm|Ie> (MI,ne)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    global_dpd_->contract444(&I, &T2, &Z, 0, 0, 1, 0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_CC_TMP0, 0, 0, 10, 0, 10, 0, "<Mn|Ie> (MI,ne)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIbjA");
    global_dpd_->contract444(&I, &T2, &Z, 0, 1, -1, 1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_HBAR, psqr, 10, 0, "WMbIj", 1);
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 0, 10, 0, 10, 0, "Z(Mj,Ib)");
    global_dpd_->buf4_init(&I, PSIF_CC_EINTS, 0, 0, 11, 0, 11, 0, "E <ij|ak>");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 11, 10, 11, 0, "tIbAj");
    global_dpd_->contract444(&I, &T2, &Z, 0, 0, 1, 0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_HBAR, psrq, 10, 0, "WMbIj", -1);
    global_dpd_->buf4_close(&Z);

  }
  else if(params.ref == 1) { /** ROHF **/

    /** <MN||IE> t_JN^BE **/
    global_dpd_->buf4_init(&I, PSIF_CC_TMP0, 0, 0, 10, 0, 10, 0, "I(MI,NE)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP2, 0, 0, 10, 0, 10, 0, "Z(MI,JB)");
    global_dpd_->contract444(&I, &T2, &Z, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&I);

    /** <Mn||Ie> t_Jn^Be **/
    global_dpd_->buf4_init(&I, PSIF_CC_TMP1, 0, 0, 10, 0, 10, 0, "I(MI,NE)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    global_dpd_->contract444(&I, &T2, &Z, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&I);

    /** <MN||JE> t_IN^BE **/
    global_dpd_->buf4_init(&I, PSIF_CC_TMP0, 0, 0, 10, 0, 10, 0, "I(MI,NE)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP3, 0, 0, 10, 0, 10, 0, "Z(MJ,IB)");
    global_dpd_->contract444(&I, &T2, &Z, 0, 0, -1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&I);

    /** <Mn||Je> t_In^Be **/
    global_dpd_->buf4_init(&I, PSIF_CC_TMP1, 0, 0, 10, 0, 10, 0, "I(MI,NE)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    global_dpd_->contract444(&I, &T2, &Z, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP2, 0, 0, 10, 0, 10, 0, "Z(MI,JB)");
    global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP4, prqs, 0, 10, "Z(MJ,IB)");
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP4, 0, 0, 10, 0, 10, 0, "Z(MJ,IB)");
    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP3, 0, 0, 10, 0, 10, 0, "Z(MJ,IB)");
    global_dpd_->buf4_axpy(&Z1, &Z2, 1.0);
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_close(&Z2);
    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP3, 0, 0, 10, 0, 10, 0, "Z(MJ,IB)");
    global_dpd_->buf4_sort(&Z2, PSIF_CC_TMP4, psrq, 10, 0, "Z(MB,IJ)");
    global_dpd_->buf4_close(&Z2);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP4, 0, 10, 0, 10, 0, 0, "Z(MB,IJ)");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 0, 10, 2, 0, "WMBIJ");
    global_dpd_->buf4_axpy(&Z, &W, 1.0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&W);

    /** <mn||ie> t_jn^be **/
    global_dpd_->buf4_init(&I, PSIF_CC_TMP0, 0, 0, 10, 0, 10, 0, "I(MI,NE)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP2, 0, 0, 10, 0, 10, 0, "Z(mi,jb)");
    global_dpd_->contract444(&I, &T2, &Z, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&I);

    /** <mN||iE> t_jN^bE **/
    global_dpd_->buf4_init(&I, PSIF_CC_TMP1, 0, 0, 10, 0, 10, 0, "I(MI,NE)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiaJB");
    global_dpd_->contract444(&I, &T2, &Z, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&I);

    /** <mn||je> t_in^be **/
    global_dpd_->buf4_init(&I, PSIF_CC_TMP0, 0, 0, 10, 0, 10, 0, "I(MI,NE)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP3, 0, 0, 10, 0, 10, 0, "Z(mj,ib)");
    global_dpd_->contract444(&I, &T2, &Z, 0, 0, -1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&I);

    /** <mN||jE> t_iN^bE **/
    global_dpd_->buf4_init(&I, PSIF_CC_TMP1, 0, 0, 10, 0, 10, 0, "I(MI,NE)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiaJB");
    global_dpd_->contract444(&I, &T2, &Z, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP2, 0, 0, 10, 0, 10, 0, "Z(mi,jb)");
    global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP4, prqs, 0, 10, "Z(mj,ib)");
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP4, 0, 0, 10, 0, 10, 0, "Z(mj,ib)");
    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP3, 0, 0, 10, 0, 10, 0, "Z(mj,ib)");
    global_dpd_->buf4_axpy(&Z1, &Z2, 1.0);
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_close(&Z2);
    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP3, 0, 0, 10, 0, 10, 0, "Z(mj,ib)");
    global_dpd_->buf4_sort(&Z2, PSIF_CC_TMP4, psrq, 10, 0, "Z(mb,ij)");
    global_dpd_->buf4_close(&Z2);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP4, 0, 10, 0, 10, 0, 0, "Z(mb,ij)");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 0, 10, 2, 0, "Wmbij");
    global_dpd_->buf4_axpy(&Z, &W, 1.0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&W);

    /** <MN||IE> t_jN^bE **/
    global_dpd_->buf4_init(&I, PSIF_CC_TMP0, 0, 0, 10, 0, 10, 0, "I(MI,NE)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiaJB");
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP2, 0, 0, 10, 0, 10, 0, "Z(MI,jb)");
    global_dpd_->contract444(&I, &T2, &Z, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP3, psrq, 10, 0, "Z(Mb,jI)");
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP3, 0, 10, 0, 10, 0, 0, "Z(Mb,jI)");
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP4, pqsr, 10, 0, "Z(Mb,Ij)");
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP4, 0, 10, 0, 10, 0, 0, "Z(Mb,Ij)");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 0, 10, 0, 0, "WMbIj");
    global_dpd_->buf4_axpy(&Z, &W, 1.0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&W);

    /** <Mn|Ie> t_jn^be **/
    global_dpd_->buf4_init(&I, PSIF_CC_TMP1, 0, 0, 10, 0, 10, 0, "I(MI,NE)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP2, 0, 0, 10, 0, 10, 0, "Z(MI,jb)");
    global_dpd_->contract444(&I, &T2, &Z, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP3, psrq, 10, 0, "Z(Mb,jI)");
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP3, 0, 10, 0, 10, 0, 0, "Z(Mb,jI)");
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP4, pqsr, 10, 0, "Z(Mb,Ij)");
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP4, 0, 10, 0, 10, 0, 0, "Z(Mb,Ij)");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 0, 10, 0, 0, "WMbIj");
    global_dpd_->buf4_axpy(&Z, &W, 1.0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&W);

    /** <mn||ie> t_Jn^Be **/
    global_dpd_->buf4_init(&I, PSIF_CC_TMP0, 0, 0, 10, 0, 10, 0, "I(MI,NE)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP2, 0, 0, 10, 0, 10, 0, "Z(mi,JB)");
    global_dpd_->contract444(&I, &T2, &Z, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP3, psrq, 10, 0, "Z(mB,Ji)");
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP3, 0, 10, 0, 10, 0, 0, "Z(mB,Ji)");
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP4, pqsr, 10, 0, "Z(mB,iJ)");
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP4, 0, 10, 0, 10, 0, 0, "Z(mB,iJ)");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 0, 10, 0, 0, "WmBiJ");
    global_dpd_->buf4_axpy(&Z, &W, 1.0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&W);

    /** <mN|iE> t_JN^BE **/
    global_dpd_->buf4_init(&I, PSIF_CC_TMP1, 0, 0, 10, 0, 10, 0, "I(MI,NE)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP2, 0, 0, 10, 0, 10, 0, "Z(mi,JB)");
    global_dpd_->contract444(&I, &T2, &Z, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP3, psrq, 10, 0, "Z(mB,Ji)");
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP3, 0, 10, 0, 10, 0, 0, "Z(mB,Ji)");
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP4, pqsr, 10, 0, "Z(mB,iJ)");
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP4, 0, 10, 0, 10, 0, 0, "Z(mB,iJ)");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 0, 10, 0, 0, "WmBiJ");
    global_dpd_->buf4_axpy(&Z, &W, 1.0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&W);
  } /** RHF or ROHF **/
  else if(params.ref == 2) { /** UHF **/

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 0, 20, 0, 20, 0, "Z(MI,JB)");

    /** <MN||IE> t_JN^BE **/
    global_dpd_->buf4_init(&I, PSIF_CC_TMP0, 0, 0, 20, 0, 20, 0, "<MN||IE> (MI,NE)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 20, 20, 20, 20, 0, "tIAJB");
    global_dpd_->contract444(&I, &T2, &Z, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&I);

    /** <Mn|Ie> t_Jn^Be **/
    global_dpd_->buf4_init(&I, PSIF_CC_TMP0, 0, 0, 30, 0, 30, 0, "<Mn|Ie> (MI,ne)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
    global_dpd_->contract444(&I, &T2, &Z, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&I);

    /** Z(MI,JB) --> Z(MB,IJ) **/
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, psqr, 20, 0, "Z(MB,IJ)");
    global_dpd_->buf4_close(&Z);

    /** Z(MB,IJ) = Z(MB,IJ) - Z(MB,JI) **/
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 20, 0, 20, 0, 0, "Z(MB,IJ)");
    global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP0, pqsr, 20, 0, "Z(MB,JI)");
    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 20, 0, 20, 0, 0, "Z(MB,JI)");
    global_dpd_->buf4_axpy(&Z2, &Z1, -1.0);
    global_dpd_->buf4_close(&Z2);

    /** Z(MB,IJ) --> W(MB,IJ) **/
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 20, 0, 20, 2, 0, "WMBIJ");
    global_dpd_->buf4_axpy(&Z1, &W, 1.0);
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_close(&W);


    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 30, 10, 30, 0, "Z(mi,jb)");

    /** <mn||ie> t_jn^be **/
    global_dpd_->buf4_init(&I, PSIF_CC_TMP0, 0, 10, 30, 10, 30, 0, "<mn||ie> (mi,ne)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 30, 30, 30, 30, 0, "tiajb");
    global_dpd_->contract444(&I, &T2, &Z, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&I);

    /** <mN||iE> t_jN^bE **/
    global_dpd_->buf4_init(&I, PSIF_CC_TMP0, 0, 10, 20, 10, 20, 0, "<mN|iE> (mi,NE)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 30, 20, 30, 20, 0, "tiaJB");
    global_dpd_->contract444(&I, &T2, &Z, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&I);

    /** Z(mi,jb) --> Z(mb,ij) **/
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, psqr, 30, 10, "Z(mb,ij)");
    global_dpd_->buf4_close(&Z);

    /** Z(mb,ij) = Z(mb,ij) - Z(mb,ji) **/
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 30, 10, 30, 10, 0, "Z(mb,ij)");
    global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP0, pqsr, 30, 10, "Z(mb,ji)");
    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 30, 10, 30, 10, 0, "Z(mb,ji)");
    global_dpd_->buf4_axpy(&Z2, &Z1, -1.0);
    global_dpd_->buf4_close(&Z2);

    /** Z(mb,ij) --> W(mb,ij) **/
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 30, 10, 30, 12, 0, "Wmbij");
    global_dpd_->buf4_axpy(&Z1, &W, 1.0);
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_close(&W);


    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 0, 30, 0, 30, 0, "Z(MI,jb)");

    /** <MN||IE> t_jN^bE **/
    global_dpd_->buf4_init(&I, PSIF_CC_TMP0, 0, 0, 20, 0, 20, 0, "<MN||IE> (MI,NE)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 30, 20, 30, 20, 0, "tiaJB");
    global_dpd_->contract444(&I, &T2, &Z, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&I);

    /** <Mn|Ie> t_jn^be **/
    global_dpd_->buf4_init(&I, PSIF_CC_TMP0, 0, 0, 30, 0, 30, 0, "<Mn|Ie> (MI,ne)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 30, 30, 30, 30, 0, "tiajb");
    global_dpd_->contract444(&I, &T2, &Z, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&I);

    /** Z(MI,jb) --> Z(Mb,Ij) **/
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, psqr, 24, 22, "Z(Mb,Ij)");
    global_dpd_->buf4_close(&Z);


    /** -<Mn|Ej> t_In^Eb **/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 22, 24, 22, 24, 0, "Z(Mj,Ib)");
    global_dpd_->buf4_init(&I, PSIF_CC_EINTS, 0, 22, 26, 22, 26, 0, "E <Ij|Ak>");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 24, 26, 24, 26, 0, "tIbAj");
    global_dpd_->contract444(&I, &T2, &Z, 0, 0, -1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, psrq, 24, 22, "Z1(Mb,Ij)");
    global_dpd_->buf4_close(&Z);

    /** Z(Mb,Ij) --> W(Mb,Ij) **/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 24, 22, 24, 22, 0, "Z(Mb,Ij)");
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 24, 22, 24, 22, 0, "Z1(Mb,Ij)");
    global_dpd_->buf4_axpy(&Z1, &Z, 1.0);
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 24, 22, 24, 22, 0, "WMbIj");
    global_dpd_->buf4_axpy(&Z, &W, 1.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Z);


    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 20, 10, 20, 0, "Z(mi,JB)");

    /** <mn||ie> t_Jn^Be **/
    global_dpd_->buf4_init(&I, PSIF_CC_TMP0, 0, 10, 30, 10, 30, 0, "<mn||ie> (mi,ne)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
    global_dpd_->contract444(&I, &T2, &Z, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&I);

    /** <mN|iE> t_JN^BE **/
    global_dpd_->buf4_init(&I, PSIF_CC_TMP0, 0, 10, 20, 10, 20, 0, "<mN|iE> (mi,NE)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 20, 20, 20, 20, 0, "tIAJB");
    global_dpd_->contract444(&I, &T2, &Z, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&I);

    /** Z(mi,JB) --> Z(mB,iJ) **/
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, psqr, 27, 23, "Z(mB,iJ)");
    global_dpd_->buf4_close(&Z);

    /** Z(mB,iJ) --> W(mB,iJ) **/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 27, 23, 27, 23, 0, "Z(mB,iJ)");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 27, 23, 27, 23, 0, "WmBiJ");
    global_dpd_->buf4_axpy(&Z, &W, 1.0);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_close(&Z);

    /** -<mN|eJ> t_iN^eB **/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 23, 27, 23, 27, 0, "Z(mJ,iB)");
    global_dpd_->buf4_init(&I, PSIF_CC_EINTS, 0, 23, 25, 23, 25, 0, "E <iJ|aK>");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 27, 25, 27, 25, 0, "tiBaJ");
    global_dpd_->contract444(&I, &T2, &Z, 0, 0, -1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, psrq, 27, 23, "Z(mB,iJ)");
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 27, 23, 27, 23, 0, "Z(mB,iJ)");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 27, 23, 27, 23, 0, "WmBiJ");
    global_dpd_->buf4_axpy(&Z, &W, 1.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Z);

  } /** UHF **/


  if(params.ref == 1) { /** RHF or ROHF **/

    /* Sort <ai||jk> integrals for remaining E*T2 contributions */
    /** THIS IS NOT ACTUALLY NECESSARY!  FIX THESE CONTRACTIONS! (11/14/01) **/
    global_dpd_->buf4_init(&I, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    global_dpd_->buf4_sort(&I, PSIF_CC_TMP0, rspq, 0, 11, "I(Mn,Ej)");
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_init(&I, PSIF_CC_TMP0, 0, 0, 11, 0, 11, 0, "I(Mn,Ej)");
    global_dpd_->buf4_sort(&I, PSIF_CC_TMP1, psrq, 0, 11, "I(Mj,En)");
    global_dpd_->buf4_close(&I);


    /** -<Mn|Ej> t_In^EB **/ /** Can be written as: - <Mj|En> t_In^Eb !! **/
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->buf4_sort(&T2, PSIF_CC_TMP0, psrq, 10, 11, "T2(Ib,En)");
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&I, PSIF_CC_TMP1, 0, 0, 11, 0, 11, 0, "I(Mj,En)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "T2(Ib,En)");
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP2, 0, 0, 10, 0, 10, 0, "Z(Mj,Ib)");
    global_dpd_->contract444(&I, &T2, &Z, 0, 0, -1.0, 0.0);
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, psrq, 10, 0, "Z(Mb,Ij)");
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Z(Mb,Ij)");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 0, 10, 0, 0, "WMbIj");
    global_dpd_->buf4_axpy(&Z, &W, 1.0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&W);


    /** -<mN|eJ> t_iN^eB **/ /** Can be written as: -<mJ|eN> t_iN^Eb !! **/
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
    global_dpd_->buf4_sort(&T2, PSIF_CC_TMP0, psrq, 10, 11, "T2(iB,eN)");
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&I, PSIF_CC_TMP1, 0, 0, 11, 0, 11, 0, "I(Mj,En)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "T2(iB,eN)");
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP2, 0, 0, 10, 0, 10, 0, "Z(mJ,iB)");
    global_dpd_->contract444(&I, &T2, &Z, 0, 0, -1.0, 0.0);
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, psrq, 10, 0, "Z(mB,iJ)");
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Z(mB,iJ)");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 0, 10, 0, 0, "WmBiJ");
    global_dpd_->buf4_axpy(&Z, &W, 1.0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&W);
  }

  /** Prepare intermediates for final term of Wmbij **/

  if(params.ref == 0) { /** RHF **/

    /* Z(ME,jb) = { <Mb|Ej> + t_jN^bF [2 <Mn|Ef> - <Mn|Fe>] - t_jN^Fb <Mn|Ef> } */
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(ME,jb)");

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D 2<ij|ab> - <ij|ba> (ia,jb)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    global_dpd_->contract444(&D, &T2, &Z, 0, 0, 1, 0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIbjA");
    global_dpd_->contract444(&D, &T2, &Z, 0, 0, -1, 1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    global_dpd_->buf4_axpy(&D, &Z, 1);
    global_dpd_->buf4_close(&D);

    /* W(Mb,Ij) <-- Z(ME,jb) t_I^E */
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 0, 10, 0, 10, 0, "Z(MI,jb)");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract424(&Z, &T1, &Z1, 1, 1, 1, 1, 0);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_sort_axpy(&Z1, PSIF_CC_HBAR, psqr, 10, 0, "WMbIj", 1);
    global_dpd_->buf4_close(&Z1);


    /* Z(Me,Ib) = { - <Mb|Ie> + t_In^Fb <Mn|Fe> } */
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(Me,Ib)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIbjA");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
    global_dpd_->contract444(&D, &T2, &Z, 0, 0, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 11, 10, 11, 10, 0, "C <ia|jb> (bi,ja)");
    global_dpd_->buf4_sort_axpy(&C, PSIF_CC_TMP0, qprs, 10, 10, "Z(Me,Ib)", -1);
    global_dpd_->buf4_close(&C);

    /* W(Mb,Ij) <-- Z(Me,Ib) t_j^e */
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 0, 10, 0, 10, 0, "Z(Mj,Ib)");
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(Me,Ib)");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract424(&Z, &T1, &Z1, 1, 1, 1, 1, 0);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_sort_axpy(&Z1, PSIF_CC_HBAR, psrq, 10, 0, "WMbIj", -1);
    global_dpd_->buf4_close(&Z1);
  }
  else if(params.ref == 1) { /** ROHF **/

    /** t_JN^BF <MN||EF> --> Z_MBJE **/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(ME,JB)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
    global_dpd_->contract444(&D, &T2, &Z, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    global_dpd_->contract444(&D, &T2, &Z, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP1, psrq, 10, 10, "Z(MB,JE)");
    global_dpd_->buf4_close(&Z);

    /** t_I^E ( <MB||JE> + Z1_MBJE ) --> Z2_MBIJ **/
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, 0, 10, 10, 10, 10, 0, "Z(MB,JE)");
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
    global_dpd_->buf4_axpy(&C, &Z1, -1.0);
    global_dpd_->buf4_close(&C);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Z1(MB,JI)");
    global_dpd_->contract424(&Z1, &T1, &Z2, 3, 1, 0, 1.0, 0.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_sort(&Z2, PSIF_CC_TMP1, pqsr, 10, 0, "Z2(MB,IJ)");
    global_dpd_->buf4_close(&Z2);

    /** Z1_MBJI(TMP0) - Z2_MBIJ(TMP1) --> W_MBIJ **/
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Z1(MB,JI)");
    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP1, 0, 10, 0, 10, 0, 0, "Z2(MB,IJ)");
    global_dpd_->buf4_axpy(&Z1, &Z2, -1.0);
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 0, 10, 2, 0, "WMBIJ");
    global_dpd_->buf4_axpy(&Z2, &W, 1.0);
    global_dpd_->buf4_close(&Z2);
    global_dpd_->buf4_close(&W);

    /** t_jn^bf <mn||ef> --> Z_mbje **/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(me,jb)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
    global_dpd_->contract444(&D, &T2, &Z, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiaJB");
    global_dpd_->contract444(&D, &T2, &Z, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP1, psrq, 10, 10, "Z(mb,je)");
    global_dpd_->buf4_close(&Z);

    /** t_i^e ( <mb||je> + Z1_mbje ) --> Z2_mbij **/
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, 0, 10, 10, 10, 10, 0, "Z(mb,je)");
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
    global_dpd_->buf4_axpy(&C, &Z1, -1.0);
    global_dpd_->buf4_close(&C);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Z1(mb,ji)");
    global_dpd_->contract424(&Z1, &T1, &Z2, 3, 1, 0, 1.0, 0.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_sort(&Z2, PSIF_CC_TMP1, pqsr, 10, 0, "Z2(mb,ij)");
    global_dpd_->buf4_close(&Z2);

    /** Z1_mbji(TMP0) - Z2_mbij(TMP1) --> W_mbij **/
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Z1(mb,ji)");
    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP1, 0, 10, 0, 10, 0, 0, "Z2(mb,ij)");
    global_dpd_->buf4_axpy(&Z1, &Z2, -1.0);
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 0, 10, 2, 0, "Wmbij");
    global_dpd_->buf4_axpy(&Z2, &W, 1.0);
    global_dpd_->buf4_close(&Z2);
    global_dpd_->buf4_close(&W);


    /** <Mn|Ef> t_jn^bf + <MN||EF> t_jN^bF --> Z1_MEjb(TMP0) **/
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z1(ME,jb)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
    global_dpd_->contract444(&D, &T2, &Z1, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiaJB");
    global_dpd_->contract444(&D, &T2, &Z1, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&Z1);

    /** <Mn|Fe> t_In^Fb --> Z2_MeIb (TMP3) **/
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->buf4_sort(&D, PSIF_CC_TMP1, psrq, 10, 11, "D(Me,Fn)");
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->buf4_sort(&T2, PSIF_CC_TMP2, psrq, 10, 11, "T2(Ib,Fn)");
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP3, 0, 10, 10, 10, 10, 0, "Z(Me,Ib)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP2, 0, 10, 11, 10, 11, 0, "T2(Ib,Fn)");
    global_dpd_->buf4_init(&D, PSIF_CC_TMP1, 0, 10, 11, 10, 11, 0, "D(Me,Fn)");
    global_dpd_->contract444(&D, &T2, &Z2, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&Z2);

    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z1(ME,jb)");
    global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP1, psrq, 10, 10, "Z1(Mb,jE)");
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP3, 0, 10, 10, 10, 10, 0, "Z(Me,Ib)");
    global_dpd_->buf4_sort(&Z2, PSIF_CC_TMP0, psrq, 10, 10, "Z(Mb,Ie)");
    global_dpd_->buf4_close(&Z2);

    /** t_I^E ( <Mj|Eb> + Z(Mb,jE)(TMP1) ) --> Z(Mb,jI)(TMP1) **/
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, 0, 10, 10, 10, 10, 0, "Z1(Mb,jE)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
    global_dpd_->buf4_axpy(&D, &Z1, 1.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP2, 0, 10, 0, 10, 0, 0, "Z(Mb,jI)");
    global_dpd_->contract424(&Z1, &T1, &Z, 3, 1, 0, 1.0, 0.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP1, pqsr, 10, 0, "Z(Mb,Ij)");
    global_dpd_->buf4_close(&Z);

    /** t_j^e ( <Mb|Ie> - Z(Mb,Ie)(TMP0) ) --> Z(Mb,Ij)(TMP2) **/
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(Mb,Ie)");
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    global_dpd_->buf4_axpy(&C, &Z1, -1.0);
    global_dpd_->buf4_close(&C);
    global_dpd_->buf4_scm(&Z1, -1.0);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP2, 0, 10, 0, 10, 0, 0, "Z(Mb,Ij)");
    global_dpd_->contract424(&Z1, &T1, &Z, 3, 1, 0, 1.0, 0.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_close(&Z);

    /** Z(Mb,Ij) (TMP1) + Z(Mb,Ij) (TMP2) --> W_MbIj **/
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, 0, 10, 0, 10, 0, 0, "Z(Mb,Ij)");
    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP2, 0, 10, 0, 10, 0, 0, "Z(Mb,Ij)");
    global_dpd_->buf4_axpy(&Z2, &Z1, 1.0);
    global_dpd_->buf4_close(&Z2);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 0, 10, 0, 0, "WMbIj");
    global_dpd_->buf4_axpy(&Z1, &W, 1.0);
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_close(&W);

    /** t_JN^BF <mN|eF> + t_Jn^Bf <mn||ef> --> Z(mB,Je) (TMP1) **/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(me,JB)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
    global_dpd_->contract444(&D, &T2, &Z, 0 , 0, 1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    global_dpd_->contract444(&D, &T2, &Z, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP1, psrq, 10, 10, "Z(mB,Je)");
    global_dpd_->buf4_close(&Z);

    /** -t_Ni^Bf <mN|fE> --> Z(mB,iE) (TMP0) **/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP2, 0, 10, 10, 10, 10, 0, "Z(mE,iB)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tjAIb");
    global_dpd_->contract444(&D, &T2, &Z, 0, 0, -1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, psrq, 10, 10, "Z(mB,iE)");
    global_dpd_->buf4_close(&Z);

    /** t_i^e ( <mJ|eB> + Z(mB,Je) ) --> Z1(mB,iJ) (TMP1) **/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, 0, 10, 10, 10, 10, 0, "Z(mB,Je)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
    global_dpd_->buf4_axpy(&D, &Z, 1.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP2, 0, 10, 0, 10, 0, 0, "Z(mB,Ji)");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
    global_dpd_->contract424(&Z, &T1, &Z1, 3, 1, 0, 1.0, 0.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP1, pqsr, 10, 0, "Z1(mB,iJ)");
    global_dpd_->buf4_close(&Z1);

    /** t_J^E ( <mB|iE> + Z(mB,iE) ) + Z1(mB,Ij) (TMP1) --> Z2(mB,iJ) (TMP2) **/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(mB,iE)");
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    global_dpd_->buf4_axpy(&C, &Z, 1.0);
    global_dpd_->buf4_close(&C);
    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP2, 0, 10, 0, 10, 0, 0, "Z2(mB,iJ)");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract424(&Z, &T1, &Z2, 3, 1, 0, 1.0, 0.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, 0, 10, 0, 10, 0, 0, "Z1(mB,iJ)");
    global_dpd_->buf4_axpy(&Z1, &Z2, 1.0);
    global_dpd_->buf4_close(&Z1);

    /** Z2(mB,iJ) --> W_mBiJ **/
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 0, 10, 0, 0, "WmBiJ");
    global_dpd_->buf4_axpy(&Z2, &W, 1.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Z2);
  } /** RHF or ROHF **/
  else if(params.ref == 2) { /** UHF **/

    /** t_JN^BF <MN||EF> --> Z_MBJE **/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 20, 20, 20, 20, 0, "Z(ME,JB)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 20, 20, 20, 20, 0, "D <IJ||AB> (IA,JB)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 20, 20, 20, 20, 0, "tIAJB");
    global_dpd_->contract444(&D, &T2, &Z, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 20, 30, 20, 30, 0, "D <Ij|Ab> (IA,jb)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
    global_dpd_->contract444(&D, &T2, &Z, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, psrq, 20, 20, "Z(MB,JE)");
    global_dpd_->buf4_close(&Z);

    /** t_I^E ( <MB||JE> + Z1_MBJE ) --> Z2_MBIJ **/
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 20, 20, 20, 20, 0, "Z(MB,JE)");
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 20, 20, 20, 20, 0, "C <IA||JB>");
    global_dpd_->buf4_axpy(&C, &Z1, -1.0);
    global_dpd_->buf4_close(&C);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 20, 0, 20, 0, 0, "Z1(MB,JI)");
    global_dpd_->contract424(&Z1, &T1, &Z2, 3, 1, 0, 1.0, 0.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_sort(&Z2, PSIF_CC_TMP0, pqsr, 20, 0, "Z2(MB,IJ)");
    global_dpd_->buf4_close(&Z2);

    /** Z1_MBJI - Z2_MBIJ --> W_MBIJ **/
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 20, 0, 20, 0, 0, "Z1(MB,JI)");
    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 20, 0, 20, 0, 0, "Z2(MB,IJ)");
    global_dpd_->buf4_axpy(&Z1, &Z2, -1.0);
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 20, 0, 20, 2, 0, "WMBIJ");
    global_dpd_->buf4_axpy(&Z2, &W, 1.0);
    global_dpd_->buf4_close(&Z2);
    global_dpd_->buf4_close(&W);

    /** t_jn^bf <mn||ef> --> Z_mbje **/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 30, 30, 30, 30, 0, "Z(me,jb)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 30, 30, 30, 30, 0, "D <ij||ab> (ia,jb)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 30, 30, 30, 30, 0, "tiajb");
    global_dpd_->contract444(&D, &T2, &Z, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 30, 20, 30, 20, 0, "D <Ij|Ab> (ia,JB)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 30, 20, 30, 20, 0, "tiaJB");
    global_dpd_->contract444(&D, &T2, &Z, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, psrq, 30, 30, "Z(mb,je)");
    global_dpd_->buf4_close(&Z);

    /** t_i^e ( <mb||je> + Z1_mbje ) --> Z2_mbij **/
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 30, 30, 30, 30, 0, "Z(mb,je)");
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 30, 30, 30, 30, 0, "C <ia||jb>");
    global_dpd_->buf4_axpy(&C, &Z1, -1.0);
    global_dpd_->buf4_close(&C);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 30, 10, 30, 10, 0, "Z1(mb,ji)");
    global_dpd_->contract424(&Z1, &T1, &Z2, 3, 1, 0, 1.0, 0.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_sort(&Z2, PSIF_CC_TMP0, pqsr, 30, 10, "Z2(mb,ij)");
    global_dpd_->buf4_close(&Z2);

    /** Z1_mbji - Z2_mbij --> W_mbij **/
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 30, 10, 30, 10, 0, "Z1(mb,ji)");
    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 30, 10, 30, 10, 0, "Z2(mb,ij)");
    global_dpd_->buf4_axpy(&Z1, &Z2, -1.0);
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 30, 10, 30, 12, 0, "Wmbij");
    global_dpd_->buf4_axpy(&Z2, &W, 1.0);
    global_dpd_->buf4_close(&Z2);
    global_dpd_->buf4_close(&W);


    /** <Mn|Ef> t_jn^bf + <MN||EF> t_jN^bF --> Z1_MEjb **/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 20, 30, 20, 30, 0, "Z(ME,jb)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 20, 30, 20, 30, 0, "D <Ij|Ab> (IA,jb)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 30, 30, 30, 30, 0, "tiajb");
    global_dpd_->contract444(&D, &T2, &Z, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 20, 20, 20, 20, 0, "D <IJ||AB> (IA,JB)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 30, 20, 30, 20, 0, "tiaJB");
    global_dpd_->contract444(&D, &T2, &Z, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, psrq, 24, 27, "Z(Mb,jE)");
    global_dpd_->buf4_close(&Z);


    /** <Mn|Fe> t_In^Fb --> Z2_MeIb **/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 24, 24, 24, 24, 0, "Z(Me,Ib)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 24, 27, 24, 27, 0, "D <Ij|Ab> (Ib,jA)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 24, 27, 24, 27, 0, "tIbjA");
    global_dpd_->contract444(&D, &T2, &Z, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, psrq, 24, 24, "Z(Mb,Ie)");
    global_dpd_->buf4_close(&Z);


    /** t_I^E ( <Mj|Eb> + Z(Mb,jE) ) --> Z(Mb,jI) **/
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 24, 27, 24, 27, 0, "Z(Mb,jE)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 24, 27, 24, 27, 0, "D <Ij|Ab> (Ib,jA)");
    global_dpd_->buf4_axpy(&D, &Z1, 1.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 24, 23, 24, 23, 0, "Z(Mb,jI)");
    global_dpd_->contract424(&Z1, &T1, &Z, 3, 1, 0, 1.0, 0.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, pqsr, 24, 22, "Z1(Mb,Ij)");
    global_dpd_->buf4_close(&Z);

    /** t_j^e ( <Mb|Ie> - Z(Mb,Ie) ) --> Z(Mb,Ij) **/
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 24, 24, 24, 24, 0, "Z(Mb,Ie)");
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 24, 24, 24, 24, 0, "C <Ia|Jb>");
    global_dpd_->buf4_axpy(&C, &Z1, -1.0);
    global_dpd_->buf4_close(&C);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 24, 22, 24, 22, 0, "Z2(Mb,Ij)");
    global_dpd_->contract424(&Z1, &T1, &Z, 3, 1, 0, 1.0, 0.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_close(&Z);

    /** Z(Mb,Ij) (TMP1) + Z(Mb,Ij) (TMP2) --> W_MbIj **/
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 24, 22, 24, 22, 0, "Z1(Mb,Ij)");
    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 24, 22, 24, 22, 0, "Z2(Mb,Ij)");
    global_dpd_->buf4_axpy(&Z2, &Z1, -1.0);
    global_dpd_->buf4_close(&Z2);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 24, 22, 24, 22, 0, "WMbIj");
    global_dpd_->buf4_axpy(&Z1, &W, 1.0);
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_close(&W);



    /** t_JN^BF <mN|eF> + t_Jn^Bf <mn||ef> --> Z(mB,Je) **/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 30, 20, 30, 20, 0, "Z(me,JB)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 30, 20, 30, 20, 0, "D <Ij|Ab> (ia,JB)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 20, 20, 20, 20, 0, "tIAJB");
    global_dpd_->contract444(&D, &T2, &Z, 0 , 0, 1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 30, 30, 30, 30, 0, "D <ij||ab> (ia,jb)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
    global_dpd_->contract444(&D, &T2, &Z, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, psrq, 27, 24, "Z(mB,Je)");
    global_dpd_->buf4_close(&Z);

    /** t_Ni^Bf <mN|fE> --> Z(mB,iE) **/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 27, 27, 27, 27, 0, "Z(mE,iB)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 27, 24, 27, 24, 0, "D <iJ|aB> (iB,Ja)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 27, 24, 27, 24, 0, "tjAIb");
    global_dpd_->contract444(&D, &T2, &Z, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, psrq, 27, 27, "Z(mB,iE)");
    global_dpd_->buf4_close(&Z);

    /** t_i^e ( <mJ|eB> + Z(mB,Je) ) --> Z1(mB,iJ) **/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 27, 24, 27, 24, 0, "Z(mB,Je)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 27, 24, 27, 24, 0, "D <iJ|aB> (iB,Ja)");
    global_dpd_->buf4_axpy(&D, &Z, 1.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 27, 22, 27, 22, 0, "Z(mB,Ji)");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
    global_dpd_->contract424(&Z, &T1, &Z1, 3, 1, 0, 1.0, 0.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP0, pqsr, 27, 23, "Z1(mB,iJ)");
    global_dpd_->buf4_close(&Z1);

    /** -t_J^E ( -<mB|iE> + Z(mB,iE) ) + Z1(mB,Ij) --> Z1(mB,iJ) **/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 27, 27, 27, 27, 0, "Z(mB,iE)");
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 27, 27, 27, 27, 0, "C <iA|jB>");
    global_dpd_->buf4_axpy(&C, &Z, -1.0);
    global_dpd_->buf4_close(&C);
    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 27, 23, 27, 23, 0, "Z2(mB,iJ)");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract424(&Z, &T1, &Z2, 3, 1, 0, 1.0, 0.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 27, 23, 27, 23, 0, "Z1(mB,iJ)");
    global_dpd_->buf4_axpy(&Z2, &Z1, -1.0);
    global_dpd_->buf4_close(&Z2);

    /** Z2(mB,iJ) --> W_mBiJ **/
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 27, 23, 27, 23, 0, "WmBiJ");
    global_dpd_->buf4_axpy(&Z1, &W, 1.0);
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_close(&W);

  }
}

}} // namespace psi::cchbar
