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
    \ingroup CCENERGY
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "Params.h"
#include "ccwave.h"

namespace psi { namespace ccenergy {

void CCEnergyWavefunction::Z_build(void)
{
  dpdbuf4 ZIJMA, Zijma, ZIjMa, ZIjmA, ZIjAm, ZMaIj, ZmAIj, Z;
  dpdbuf4 tauIJAB, tauijab, tauIjAb, tauIjbA, F_anti, F, tau;
  int Gmb, Gij, mb, nrows, ncols;

#ifdef TIME_CCENERGY
  timer_on("Z");
#endif

  if(params_.ref == 0) { /** RHF **/
    /* ZMbIj = <Mb|Ef> * tau(Ij,Ef) */
    /* OOC code added 3/23/05  -TDC */
    global_dpd_->buf4_init(&Z, PSIF_CC_MISC, 0, 10, 0, 10, 0, 0, "ZMbIj");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    global_dpd_->buf4_init(&tau, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
    global_dpd_->contract444(&F, &tau, &Z, 0, 0, 1, 0);
/*     for(Gmb=0; Gmb < moinfo.nirreps; Gmb++) { */
/*       Gij = Gmb;  /\* tau is totally symmetric *\/ */
/*       dpd_buf4_mat_irrep_init(&tau, Gij); */
/*       dpd_buf4_mat_irrep_rd(&tau, Gij); */

/*       dpd_buf4_mat_irrep_init(&Z, Gmb); */

/*       dpd_buf4_mat_irrep_row_init(&F, Gmb); */

/*       for(mb=0; mb < F.params->rowtot[Gmb]; mb++) { */
/* 	dpd_buf4_mat_irrep_row_rd(&F, Gmb, mb); */

/* 	nrows = tau.params->rowtot[Gij]; */
/* 	ncols = tau.params->coltot[Gij]; */
/* 	if(nrows && ncols) */
/* 	  C_DGEMV('n',nrows,ncols,1.0,tau.matrix[Gij][0],ncols,F.matrix[Gmb][0],1, */
/* 		  0.0,Z.matrix[Gmb][mb],1); */
/*       } */

/*       dpd_buf4_mat_irrep_row_close(&F, Gmb); */
/*       dpd_buf4_mat_irrep_wrt(&Z, Gmb); */
/*       dpd_buf4_mat_irrep_close(&Z, Gmb); */
/*       dpd_buf4_mat_irrep_close(&tau, Gij); */
/*     } */
    global_dpd_->buf4_close(&tau);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&Z);
  }
  else if(params_.ref == 1) { /** ROHF **/
    global_dpd_->buf4_init(&ZIJMA, PSIF_CC_MISC, 0, 2, 10, 2, 10, 0, "ZIJMA");
    global_dpd_->buf4_init(&Zijma, PSIF_CC_MISC, 0, 2, 10, 2, 10, 0, "Zijma");
    global_dpd_->buf4_init(&ZIjMa, PSIF_CC_MISC, 0, 0, 10, 0, 10, 0, "ZIjMa");
    global_dpd_->buf4_init(&ZIjmA, PSIF_CC_MISC, 0, 0, 10, 0, 10, 0, "ZIjmA");

    global_dpd_->buf4_init(&tauIJAB, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
    global_dpd_->buf4_init(&tauijab, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauijab");
    global_dpd_->buf4_init(&tauIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
    global_dpd_->buf4_init(&tauIjbA, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjbA");

    global_dpd_->buf4_init(&F_anti, PSIF_CC_FINTS, 0, 10, 7, 10, 5, 1, "F <ia|bc>");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");

    global_dpd_->contract444(&tauIJAB, &F_anti, &ZIJMA, 0, 0, 1, 0);
    global_dpd_->contract444(&tauijab, &F_anti, &Zijma, 0, 0, 1, 0);
    global_dpd_->contract444(&tauIjAb, &F, &ZIjMa, 0, 0, 1, 0);
    global_dpd_->contract444(&tauIjbA, &F, &ZIjmA, 0, 0, 1, 0);

    global_dpd_->buf4_close(&tauIJAB);
    global_dpd_->buf4_close(&tauijab);
    global_dpd_->buf4_close(&tauIjAb);
    global_dpd_->buf4_close(&tauIjbA);

    global_dpd_->buf4_close(&F_anti);
    global_dpd_->buf4_close(&F);

    global_dpd_->buf4_sort(&ZIJMA, PSIF_CC_MISC, pqsr, 2, 11, "ZIJAM");
    global_dpd_->buf4_sort(&Zijma, PSIF_CC_MISC, pqsr, 2, 11, "Zijam");
    global_dpd_->buf4_sort(&ZIjmA, PSIF_CC_MISC, pqsr, 0, 11, "ZIjAm");

    global_dpd_->buf4_close(&ZIJMA);
    global_dpd_->buf4_close(&Zijma);
    global_dpd_->buf4_close(&ZIjMa);
    global_dpd_->buf4_close(&ZIjmA);
  }
  else if(params_.ref == 2) { /*** UHF ***/

    global_dpd_->buf4_init(&ZIJMA, PSIF_CC_MISC, 0, 2, 20, 2, 20, 0, "ZIJMA");
    global_dpd_->buf4_init(&Zijma, PSIF_CC_MISC, 0, 12, 30, 12, 30, 0, "Zijma");
    global_dpd_->buf4_init(&ZIjMa, PSIF_CC_MISC, 0, 22, 24, 22, 24, 0, "ZIjMa");
    global_dpd_->buf4_init(&ZIjAm, PSIF_CC_MISC, 0, 22, 26, 22, 26, 0, "ZIjAm");

    global_dpd_->buf4_init(&tauIJAB, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
    global_dpd_->buf4_init(&tauijab, PSIF_CC_TAMPS, 0, 12, 17, 12, 17, 0, "tauijab");
    global_dpd_->buf4_init(&tauIjAb, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tauIjAb");

    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 20, 7, 20, 5, 1, "F <IA|BC>");
    global_dpd_->contract444(&tauIJAB, &F, &ZIJMA, 0, 0, 1, 0);
    global_dpd_->buf4_close(&F);

    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 30, 17, 30, 15, 1, "F <ia|bc>");
    global_dpd_->contract444(&tauijab, &F, &Zijma, 0, 0, 1, 0);
    global_dpd_->buf4_close(&F);

    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
    global_dpd_->contract444(&tauIjAb, &F, &ZIjMa, 0, 0, 1, 0);
    global_dpd_->buf4_close(&F);

    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 28, 26, 28, 26, 0, "F <Ab|Ci>");
    global_dpd_->contract444(&tauIjAb, &F, &ZIjAm, 0, 1, 1, 0);
    global_dpd_->buf4_close(&F);

    global_dpd_->buf4_close(&tauIJAB);
    global_dpd_->buf4_close(&tauijab);
    global_dpd_->buf4_close(&tauIjAb);

    global_dpd_->buf4_sort(&ZIJMA, PSIF_CC_MISC, pqsr, 2, 21, "ZIJAM");
    global_dpd_->buf4_sort(&Zijma, PSIF_CC_MISC, pqsr, 12, 31, "Zijam");

    global_dpd_->buf4_close(&ZIJMA);
    global_dpd_->buf4_close(&Zijma);
    global_dpd_->buf4_close(&ZIjMa);
    global_dpd_->buf4_close(&ZIjAm);

  }

#ifdef TIME_CCENERGY
  timer_off("Z");
#endif
}

}} // namespace psi::ccenergy
