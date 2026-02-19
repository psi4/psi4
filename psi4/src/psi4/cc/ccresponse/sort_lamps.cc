/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*! \file
    \ingroup ccresponse
    \brief Enter brief description of file here
*/
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"
#include <cmath>

namespace psi {
namespace ccresponse {

void sort_lamps() {
    dpdbuf4 L;

    /* RAK fixing this for new cclambda, assuming A1 ground lambda? */
    global_dpd_->buf4_init(&L, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "LIjAb 0 -1");
    global_dpd_->buf4_scmcopy(&L, PSIF_CC_LAMPS, "2 LIjAb - LIjBa", 2);
    global_dpd_->buf4_sort_axpy(&L, PSIF_CC_LAMPS, pqsr, 0, 5, "2 LIjAb - LIjBa", -1);
    global_dpd_->buf4_close(&L);
}

void sort_lamps_quadratic_resp() {
    dpdbuf4 L;

    //Sorted  L (ib|ja)
    global_dpd_->buf4_init(&L, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    global_dpd_->buf4_sort(&L, PSIF_CC_LAMPS, psqr, 10, 10, "(2 LIjAb - LIjBa) (ib|ja)");
    global_dpd_->buf4_sort(&L, PSIF_CC_LAMPS, prqs, 10, 10, "(2 LIjAb - LIjBa) (ia|jb)");
    global_dpd_->buf4_sort(&L, PSIF_CC_LAMPS, pqsr, 0, 5, "(2 LIjAb - LIjBa) (ij,ba)");
    global_dpd_->buf4_close(&L);
}

void sort_integrals_quadratic_resp() {
    dpdbuf4 L2, D, W, WL;
    dpdfile2 L1;
    double Y2_norm;

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->buf4_sort(&D, PSIF_CC_DINTS, prqs, 10, 10, "D (ia|jb)");
    global_dpd_->buf4_sort(&D, PSIF_CC_DINTS, psqr, 10, 10, "D (ib|ja)");
    global_dpd_->buf4_sort(&D, PSIF_CC_DINTS, pqsr, 0, 5, "D (ij|ba)");
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 11, 0, 11, 0, "2WMnIe - WnMIe (Mn,eI)");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, qprs, 0, 11, "2WMnIe - WnMIe (nM,eI)");   //sort
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 11, 0, 11, 0, "2WMnIe - WnMIe (Mn,eI)");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, psqr, 0, 10, "2WMnIe - WnMIe (MI,nE)");
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 10, 0, 10, 0, "WMnIe");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, psrq, 10, 0, "WMnIe (Me,nI)");  //Can we replace?
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 10, 0, 10, 0, "WMnIe");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, sqrp, 11, 0, "WMnIe (en,IM)");  //Can we replace?
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, qprs, 0, 0, "WMnIj (nM,Ij)");
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
    global_dpd_->buf4_scmcopy(&W, PSIF_CC_HBAR, "WAmEf 2(Am,Ef) - (Am,fE)", 2);
    global_dpd_->buf4_sort_axpy(&W, PSIF_CC_HBAR, pqsr, 11, 5, "WAmEf 2(Am,Ef) - (Am,fE)", -1);
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, pqsr, 11, 5, "WAmEf (am,fe)");
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf 2(Am,Ef) - (Am,fE)");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, qprs, 10, 5, "WAmEf 2(mA,Ef) - (mA,fE)"); //Compute this part out of core
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, pqsr, 11, 5, "WAmEf 2(Am,Ef) - (Am,fE) (am,fe)"); //Compute this part out of core
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf 2(Am,Ef) - (Am,fE)");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, prqs, 5, 10, "WAmEf 2(Am,Ef) - (Am,fE) (AE,mf)"); //Compute this part out of core
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbeJ");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, psrq, 10, 10, "WMebJ");
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 11, 10, 11, 0, "WMbEj");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, prqs, 10, 11, "WMEbj");
    global_dpd_->buf4_close(&W);

    //Build Hvvvv x L2 
    global_dpd_->buf4_init(&WL, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "WefabL2");
    global_dpd_->buf4_scmcopy(&WL, PSIF_CC_LAMPS, "WefabL2 2(ij,ab) - (ij,ba)", 2);
    global_dpd_->buf4_sort_axpy(&WL, PSIF_CC_LAMPS, pqsr, 0, 5, "WefabL2 2(ij,ab) - (ij,ba)", -1);
    global_dpd_->buf4_close(&WL);

    //Sort
    global_dpd_->buf4_init(&WL, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "WefabL2 2(ij,ab) - (ij,ba)");
//    global_dpd_->buf4_sort(&WL, PSIF_CC_LAMPS, pqsr, 0, 5, "WefabL2 2(ji,ab) - (ji,ba)");
    //global_dpd_->buf4_close(&WL);


    //Sort WefabL2 2(ij,ab) - (ij,ba) ->WefabL2 2(ij,ba) - (ij,ab)
    //global_dpd_->buf4_init(&WL, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "WefabL2 2(ij,ab) - (ij,ba)");
    global_dpd_->buf4_sort(&WL, PSIF_CC_LAMPS, pqsr, 0, 5, "WefabL2 2(ij,ba) - (ij,ab)");
    global_dpd_->buf4_close(&WL);   
/*
    //Sort
    global_dpd_->buf4_init(&WL, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "WefabL2 2(ij,ab) - (ij,ba)");
    //global_dpd_->buf4_sort(&WL, PSIF_CC_LAMPS, pqsr, 0, 5, "WefabL2 2(ji,ab) - (ji,ba)");
    //global_dpd_->buf4_close(&WL);
    //global_dpd_->buf4_scmcopy(&WL, PSIF_CC_LAMPS, "WefabL2 (ij,ab) + WefabL2 (ij,ba)", 1);
    global_dpd_->buf4_sort_axpy(&WL, PSIF_CC_HBAR, pqsr, 0, 5, "WefabL2 (ij,ab) + WefabL2 (ij,ba)", 2);
    //global_dpd_->buf4_sort_axpy(&WL, PSIF_CC_HBAR, pqsr, 0, 5, "WefabL2 (ij,ab) + WefabL2 (ij,ba)", 1);
    global_dpd_->buf4_close(&WL);

    //Sort WefabL2 2(ij,ab) - (ij,ba) ->WefabL2 2(ij,ba) - (ij,ab)
    //global_dpd_->buf4_sort(&WL, PSIF_CC_LAMPS, pqsr, 0, 5, "WefabL2 2(ij,ba) - (ij,ab)");
    //global_dpd_->buf4_close(&WL);     
*/
    
    //OPTIMIZE IT!
    global_dpd_->buf4_init(&WL, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "WL(ij,ab) + WL(ji,ba)"); 
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 5, 10, 5, 0, "WAmEf 2(mA,Ef) - (mA,fE)");
    global_dpd_->file2_init(&L1, PSIF_CC_LAMPS, 0, 0, 1, "LIA 0 -1");
    global_dpd_->contract424(&W, &L1, &WL, 1, 1, 1, 1, 0); 
    global_dpd_->buf4_sort_axpy(&WL, PSIF_CC_HBAR, qpsr, 0, 5, "WL(ij,ab) + WL(ji,ba)", 1);
    global_dpd_->buf4_close(&WL);
    global_dpd_->file2_close(&L1);

    // Type-II L2 residual 
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 10, 10, 10, 10, 0, "LHX1Y1 Residual II");
    global_dpd_->buf4_scmcopy(&L2, PSIF_CC_LAMPS, "LHX1Y1 (ia,jb) + (jb,ia) Residual II", 2);	
    global_dpd_->buf4_sort_axpy(&L2, PSIF_CC_LAMPS, rspq, 10, 10, "LHX1Y1 (ia,jb) + (jb,ia) Residual II", 2);	
    global_dpd_->buf4_close(&L2);

}

}  // namespace ccresponse
}  // namespace psi
