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
#include <cstdlib>
#include <string>
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cchbar {

/* function to compute matrix element <Phi_ij^ab|Hbar|0> which equals zero as long as the
T amplitudes were determined by ccenergy with the same H */

void DT2(void), FaetT2(void), FmitT2(void), WmnijT2(void), WmbejT2(void);
void BT2(void), ZT2(void), FT2(void), ET2(void), CT2(void), dijabT2(void);
void BT2_AO(void); void Z_build(void);
void status(const char *, std::string);
void FT2_CC2(void);
void Wmnij_for_Wabij(void);
void Wmbej_for_Wabij(void);
void purge_Wabij(void);

void Wabij_build(void)
{
  dpdbuf4 newtIJAB, newtijab, newtIjAb;
  double dotval;

  Wmbej_for_Wabij();
  Z_build();
  Wmnij_for_Wabij();

  DT2();
  if(params.print & 2) status("<ij||ab> -> T2", "outfile");

  if(params.wfn != "CC2") { /* skip all this is wfn=CC2 */

    FaetT2();
    FmitT2();
    if(params.print & 2) status("F -> T2", "outfile");

    timer_on("WmnijT2");
    WmnijT2();
    if(params.print & 2) status("Wmnij -> T2", "outfile");
    timer_off("WmnijT2");

    timer_on("BT2");
    /* didn't bother to put BT2_AO in here yet */
    BT2();
    if(params.print & 2) status("<ab||cd> -> T2", "outfile");
    timer_off("BT2");

    timer_on("ZT2");
    ZT2();
    if(params.print & 2) status("Z -> T2", "outfile");
    timer_off("ZT2");

    timer_on("FT2");
    FT2();
    if(params.print & 2) status("<ia||bc> -> T2", "outfile");
    timer_off("FT2");

    timer_on("ET2");
    ET2();
    if(params.print & 2) status("<ij||ka> -> T2", "outfile");
    timer_off("ET2");

    timer_on("WmbejT2");
    WmbejT2();
    if(params.print & 2) status("Wmbej -> T2", "outfile");
    timer_off("WmbejT2");

    timer_on("CT2");
    CT2();
    if(params.print & 2) status("<ia||jb> -> T2", "outfile");
    timer_off("CT2");

  }
  else { /* For CC2, just include (FT2)c->T2 */
    /* didn't implement CC2 yet */
  }

  if (params.ref == 1)
    purge_Wabij();

  if(params.ref == 0) { /** RHF **/
    global_dpd_->buf4_init(&newtIjAb, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "WAbIj residual");
    dotval = global_dpd_->buf4_dot_self(&newtIjAb);
    outfile->Printf("Norm squared of <Phi^IJ_AB|Hbar|0>: %20.15lf\n",dotval);
    global_dpd_->buf4_close(&newtIjAb);
  }
  else if (params.ref == 1) {
    global_dpd_->buf4_init(&newtIJAB, PSIF_CC_HBAR, 0, 2, 7, 2, 7, 0, "WABIJ residual");
    global_dpd_->buf4_print(&newtIJAB,"outfile",1);
    global_dpd_->buf4_close(&newtIJAB);
    global_dpd_->buf4_init(&newtijab, PSIF_CC_HBAR, 0, 2, 7, 2, 7, 0, "Wabij residual");
    global_dpd_->buf4_print(&newtijab,"outfile",1);
    global_dpd_->buf4_close(&newtijab);
    global_dpd_->buf4_init(&newtIjAb, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "WAbIj residual");
    global_dpd_->buf4_print(&newtIjAb,"outfile",1);
    global_dpd_->buf4_close(&newtIjAb);
  }
  else if (params.ref == 2) {
    global_dpd_->buf4_init(&newtIJAB, PSIF_CC_HBAR, 0, 2, 7, 2, 7, 0, "WABIJ residual");
    global_dpd_->buf4_print(&newtIJAB,"outfile",1);
    global_dpd_->buf4_close(&newtIJAB);
    global_dpd_->buf4_init(&newtijab, PSIF_CC_HBAR, 0, 12, 17, 12, 17, 0, "Wabij residual");
    global_dpd_->buf4_print(&newtijab,"outfile",1);
    global_dpd_->buf4_close(&newtijab);
    global_dpd_->buf4_init(&newtIjAb, PSIF_CC_HBAR, 0, 22, 28, 22, 28, 0, "WAbIj residual");
    global_dpd_->buf4_print(&newtIjAb,"outfile",1);
    global_dpd_->buf4_close(&newtIjAb);
  }
}

void DT2(void)
{
  dpdbuf4 D;
  if(params.ref == 0) { /*** RHF ***/
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->buf4_copy(&D, PSIF_CC_HBAR, "WAbIj residual");
    global_dpd_->buf4_close(&D);
  }
  else if(params.ref == 1) { /*** ROHF ***/
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
    global_dpd_->buf4_copy(&D, PSIF_CC_HBAR, "WABIJ residual");
    global_dpd_->buf4_copy(&D, PSIF_CC_HBAR, "Wabij residual");
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->buf4_copy(&D, PSIF_CC_HBAR, "WAbIj residual");
    global_dpd_->buf4_close(&D);
  }
  else if(params.ref == 2) { /*** UHF ***/
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 2, 7, 2, 7, 0, "D <IJ||AB> (I>J,A>B)");
    global_dpd_->buf4_copy(&D, PSIF_CC_HBAR, "WABIJ residual");
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 12, 17, 12, 17, 0, "D <ij||ab> (i>j,a>b)");
    global_dpd_->buf4_copy(&D, PSIF_CC_HBAR, "Wabij residual");
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    global_dpd_->buf4_copy(&D, PSIF_CC_HBAR, "WAbIj residual");
    global_dpd_->buf4_close(&D);
  }
}

void FaetT2(void)
{
  dpdfile2 FAEt, Faet;
  dpdbuf4 newtIJAB, newtijab, newtIjAb;
  dpdbuf4 tIJAB, tijab, tIjAb;
  dpdbuf4 t2;

  if(params.ref == 0) { /** RHF **/
    global_dpd_->buf4_init(&newtIjAb, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "WAbIj residual");
    global_dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->file2_init(&FAEt, PSIF_CC_OEI, 0, 1, 1, "FAEt");
    global_dpd_->contract424(&tIjAb, &FAEt, &newtIjAb, 3, 1, 0, 1, 1);
    global_dpd_->contract244(&FAEt, &tIjAb, &newtIjAb, 1, 2, 1, 1, 1);
    global_dpd_->file2_close(&FAEt);
    global_dpd_->buf4_close(&tIjAb);
    global_dpd_->buf4_close(&newtIjAb);
  }
  else if(params.ref == 1) { /** ROHF **/
    global_dpd_->buf4_init(&newtIJAB, PSIF_CC_HBAR, 0, 2, 5, 2, 7, 0, "WABIJ residual");
    global_dpd_->buf4_init(&newtijab, PSIF_CC_HBAR, 0, 2, 5, 2, 7, 0, "Wabij residual");
    global_dpd_->buf4_init(&newtIjAb, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "WAbIj residual");
    global_dpd_->buf4_init(&tIJAB, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    global_dpd_->buf4_init(&tijab, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tijab");
    global_dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->file2_init(&FAEt, PSIF_CC_OEI, 0, 1, 1, "FAEt");
    global_dpd_->file2_init(&Faet, PSIF_CC_OEI, 0, 1, 1, "Faet");
    global_dpd_->buf4_init(&t2, PSIF_CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,AB)");
    global_dpd_->contract424(&tIJAB, &FAEt, &t2, 3, 1, 0, 1, 0);
    global_dpd_->contract244(&FAEt, &tIJAB, &t2, 1, 2, 1, 1, 1);
    global_dpd_->buf4_axpy(&t2, &newtIJAB, 1);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_init(&t2, PSIF_CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,AB)");
    global_dpd_->contract424(&tijab, &Faet, &t2, 3, 1, 0, 1, 0);
    global_dpd_->contract244(&Faet, &tijab, &t2, 1, 2, 1, 1, 1);
    global_dpd_->buf4_axpy(&t2, &newtijab, 1);
    global_dpd_->buf4_close(&t2);
    global_dpd_->contract424(&tIjAb, &Faet, &newtIjAb, 3, 1, 0, 1, 1);
    global_dpd_->contract244(&FAEt, &tIjAb, &newtIjAb, 1, 2, 1, 1, 1);
    global_dpd_->file2_close(&FAEt);
    global_dpd_->file2_close(&Faet);
    global_dpd_->buf4_close(&tIJAB);
    global_dpd_->buf4_close(&tijab);
    global_dpd_->buf4_close(&tIjAb);
    global_dpd_->buf4_close(&newtIJAB);
    global_dpd_->buf4_close(&newtijab);
    global_dpd_->buf4_close(&newtIjAb);
  }
  else if(params.ref == 2) { /*** UHF ***/
    global_dpd_->buf4_init(&newtIJAB, PSIF_CC_HBAR, 0, 2, 5, 2, 7, 0, "WABIJ residual");
    global_dpd_->buf4_init(&newtijab, PSIF_CC_HBAR, 0, 12, 15, 12, 17, 0, "Wabij residual");
    global_dpd_->buf4_init(&newtIjAb, PSIF_CC_HBAR, 0, 22, 28, 22, 28, 0, "WAbIj residual");
    global_dpd_->buf4_init(&tIJAB, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    global_dpd_->buf4_init(&tijab, PSIF_CC_TAMPS, 0, 12, 15, 12, 17, 0, "tijab");
    global_dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    global_dpd_->file2_init(&FAEt, PSIF_CC_OEI, 0, 1, 1, "FAEt");
    global_dpd_->file2_init(&Faet, PSIF_CC_OEI, 0, 3, 3, "Faet");
    global_dpd_->buf4_init(&t2, PSIF_CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,AB)");
    global_dpd_->contract424(&tIJAB, &FAEt, &t2, 3, 1, 0, 1, 0);
    global_dpd_->contract244(&FAEt, &tIJAB, &t2, 1, 2, 1, 1, 1);
    global_dpd_->buf4_axpy(&t2, &newtIJAB, 1);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_init(&t2, PSIF_CC_TMP0, 0, 12, 15, 12, 15, 0, "T (i>j,ab)");
    global_dpd_->contract424(&tijab, &Faet, &t2, 3, 1, 0, 1, 0);
    global_dpd_->contract244(&Faet, &tijab, &t2, 1, 2, 1, 1, 1);
    global_dpd_->buf4_axpy(&t2, &newtijab, 1);
    global_dpd_->buf4_close(&t2);
    global_dpd_->contract424(&tIjAb, &Faet, &newtIjAb, 3, 1, 0, 1, 1);
    global_dpd_->contract244(&FAEt, &tIjAb, &newtIjAb, 1, 2, 1, 1, 1);
    global_dpd_->file2_close(&FAEt);
    global_dpd_->file2_close(&Faet);
    global_dpd_->buf4_close(&tIJAB);
    global_dpd_->buf4_close(&tijab);
    global_dpd_->buf4_close(&tIjAb);
    global_dpd_->buf4_close(&newtIJAB);
    global_dpd_->buf4_close(&newtijab);
    global_dpd_->buf4_close(&newtIjAb);
  }
}

void FmitT2(void)
{
  dpdfile2 FMIt, Fmit;
  dpdbuf4 newtIJAB, newtijab, newtIjAb;
  dpdbuf4 tIJAB, tijab, tIjAb;
  dpdbuf4 t2;

  if(params.ref == 0) { /** RHF **/
    global_dpd_->buf4_init(&newtIjAb, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "WAbIj residual");
    global_dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->file2_init(&FMIt, PSIF_CC_OEI, 0, 0, 0, "FMIt");

    global_dpd_->contract424(&tIjAb, &FMIt, &newtIjAb, 1, 0, 1, -1, 1);
    global_dpd_->contract244(&FMIt, &tIjAb, &newtIjAb, 0, 0, 0, -1, 1);

    global_dpd_->file2_close(&FMIt);
    global_dpd_->buf4_close(&tIjAb);
    global_dpd_->buf4_close(&newtIjAb);
  }
  else if(params.ref == 1) { /** ROHF **/
    global_dpd_->buf4_init(&newtIJAB, PSIF_CC_HBAR, 0, 0, 7, 2, 7, 0, "WABIJ residual");
    global_dpd_->buf4_init(&newtijab, PSIF_CC_HBAR, 0, 0, 7, 2, 7, 0, "Wabij residual");
    global_dpd_->buf4_init(&newtIjAb, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "WAbIj residual");
    global_dpd_->buf4_init(&tIJAB, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
    global_dpd_->buf4_init(&tijab, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tijab");
    global_dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->file2_init(&FMIt, PSIF_CC_OEI, 0, 0, 0, "FMIt");
    global_dpd_->file2_init(&Fmit, PSIF_CC_OEI, 0, 0, 0, "Fmit");
    global_dpd_->buf4_init(&t2, PSIF_CC_TMP0, 0, 0, 7, 0, 7, 0, "T (IJ,A>B)");
    global_dpd_->contract424(&tIJAB, &FMIt, &t2, 1, 0, 1, -1, 0);
    global_dpd_->contract244(&FMIt, &tIJAB, &t2, 0, 0, 0, -1, 1);
    global_dpd_->buf4_axpy(&t2, &newtIJAB, 1);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_init(&t2, PSIF_CC_TMP0, 0, 0, 7, 0, 7, 0, "T (IJ,A>B)");
    global_dpd_->contract424(&tijab, &Fmit, &t2, 1, 0, 1, -1, 0);
    global_dpd_->contract244(&Fmit, &tijab, &t2, 0, 0, 0, -1, 1);
    global_dpd_->buf4_axpy(&t2, &newtijab, 1);
    global_dpd_->buf4_close(&t2);
    global_dpd_->contract424(&tIjAb, &Fmit, &newtIjAb, 1, 0, 1, -1, 1);
    global_dpd_->contract244(&FMIt, &tIjAb, &newtIjAb, 0, 0, 0, -1, 1);
    global_dpd_->file2_close(&FMIt);
    global_dpd_->file2_close(&Fmit);
    global_dpd_->buf4_close(&tIJAB);
    global_dpd_->buf4_close(&tijab);
    global_dpd_->buf4_close(&tIjAb);
    global_dpd_->buf4_close(&newtIJAB);
    global_dpd_->buf4_close(&newtijab);
    global_dpd_->buf4_close(&newtIjAb);
  }
  else if(params.ref == 2) { /*** UHF ***/
    global_dpd_->buf4_init(&newtIJAB, PSIF_CC_HBAR, 0, 0, 7, 2, 7, 0, "WABIJ residual");
    global_dpd_->buf4_init(&newtijab, PSIF_CC_HBAR, 0, 10, 17, 12, 17, 0, "Wabij residual");
    global_dpd_->buf4_init(&newtIjAb, PSIF_CC_HBAR, 0, 22, 28, 22, 28, 0, "WAbIj residual");
    global_dpd_->buf4_init(&tIJAB, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
    global_dpd_->buf4_init(&tijab, PSIF_CC_TAMPS, 0, 10, 17, 12, 17, 0, "tijab");
    global_dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    global_dpd_->file2_init(&FMIt, PSIF_CC_OEI, 0, 0, 0, "FMIt");
    global_dpd_->file2_init(&Fmit, PSIF_CC_OEI, 0, 2, 2, "Fmit");
    global_dpd_->buf4_init(&t2, PSIF_CC_TMP0, 0, 0, 7, 0, 7, 0, "T (IJ,A>B)");
    global_dpd_->contract424(&tIJAB, &FMIt, &t2, 1, 0, 1, -1, 0);
    global_dpd_->contract244(&FMIt, &tIJAB, &t2, 0, 0, 0, -1, 1);
    global_dpd_->buf4_axpy(&t2, &newtIJAB, 1);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_init(&t2, PSIF_CC_TMP0, 0, 10, 17, 10, 17, 0, "T (ij,a>b)");
    global_dpd_->contract424(&tijab, &Fmit, &t2, 1, 0, 1, -1, 0);
    global_dpd_->contract244(&Fmit, &tijab, &t2, 0, 0, 0, -1, 1);
    global_dpd_->buf4_axpy(&t2, &newtijab, 1);
    global_dpd_->buf4_close(&t2);
    global_dpd_->contract424(&tIjAb, &Fmit, &newtIjAb, 1, 0, 1, -1, 1);
    global_dpd_->contract244(&FMIt, &tIjAb, &newtIjAb, 0, 0, 0, -1, 1);
    global_dpd_->file2_close(&FMIt);
    global_dpd_->file2_close(&Fmit);
    global_dpd_->buf4_close(&tIJAB);
    global_dpd_->buf4_close(&tijab);
    global_dpd_->buf4_close(&tIjAb);
    global_dpd_->buf4_close(&newtIJAB);
    global_dpd_->buf4_close(&newtijab);
    global_dpd_->buf4_close(&newtIjAb);
  }
}

void WmnijT2(void)
{
  dpdbuf4 newtIJAB, newtijab, newtIjAb;
  dpdbuf4 WMNIJ, Wmnij, WMnIj;
  dpdbuf4 tauIJAB, tauijab, tauIjAb;

  if(params.ref == 0) { /** RHF **/
    global_dpd_->buf4_init(&newtIjAb, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "WAbIj residual");
    global_dpd_->buf4_init(&WMnIj, PSIF_CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj");
    global_dpd_->buf4_init(&tauIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
    global_dpd_->contract444(&WMnIj, &tauIjAb, &newtIjAb, 1, 1, 1, 1);
    global_dpd_->buf4_close(&tauIjAb);
    global_dpd_->buf4_close(&WMnIj);
    global_dpd_->buf4_close(&newtIjAb);
  }
  else if(params.ref == 1) { /** ROHF **/
    global_dpd_->buf4_init(&newtIJAB, PSIF_CC_HBAR, 0, 2, 7, 2, 7, 0, "WABIJ residual");
    global_dpd_->buf4_init(&WMNIJ, PSIF_CC_HBAR, 0, 2, 2, 2, 2, 0, "WMNIJ");
    global_dpd_->buf4_init(&tauIJAB, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
    global_dpd_->contract444(&WMNIJ, &tauIJAB, &newtIJAB, 1, 1, 1, 1);
    global_dpd_->buf4_close(&tauIJAB);
    global_dpd_->buf4_close(&WMNIJ);
    global_dpd_->buf4_close(&newtIJAB);

    global_dpd_->buf4_init(&newtijab, PSIF_CC_HBAR, 0, 2, 7, 2, 7, 0, "Wabij residual");
    global_dpd_->buf4_init(&Wmnij, PSIF_CC_HBAR, 0, 2, 2, 2, 2, 0, "Wmnij");
    global_dpd_->buf4_init(&tauijab, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauijab");
    global_dpd_->contract444(&Wmnij, &tauijab, &newtijab, 1, 1, 1, 1);
    global_dpd_->buf4_close(&tauijab);
    global_dpd_->buf4_close(&Wmnij);
    global_dpd_->buf4_close(&newtijab);

    global_dpd_->buf4_init(&newtIjAb, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "WAbIj residual");
    global_dpd_->buf4_init(&WMnIj, PSIF_CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj");
    global_dpd_->buf4_init(&tauIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
    global_dpd_->contract444(&WMnIj, &tauIjAb, &newtIjAb, 1, 1, 1, 1);
    global_dpd_->buf4_close(&tauIjAb);
    global_dpd_->buf4_close(&WMnIj);
    global_dpd_->buf4_close(&newtIjAb);
  }
  else if(params.ref == 2) { /*** UHF ***/
    global_dpd_->buf4_init(&newtIJAB, PSIF_CC_HBAR, 0, 2, 7, 2, 7, 0, "WABIJ residual");
    global_dpd_->buf4_init(&WMNIJ, PSIF_CC_HBAR, 0, 2, 2, 2, 2, 0, "WMNIJ");
    global_dpd_->buf4_init(&tauIJAB, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
    global_dpd_->contract444(&WMNIJ, &tauIJAB, &newtIJAB, 1, 1, 1, 1);
    global_dpd_->buf4_close(&tauIJAB);
    global_dpd_->buf4_close(&WMNIJ);
    global_dpd_->buf4_close(&newtIJAB);

    global_dpd_->buf4_init(&newtijab, PSIF_CC_HBAR, 0, 12, 17, 12, 17, 0, "Wabij residual");
    global_dpd_->buf4_init(&Wmnij, PSIF_CC_HBAR, 0, 12, 12, 12, 12, 0, "Wmnij");
    global_dpd_->buf4_init(&tauijab, PSIF_CC_TAMPS, 0, 12, 17, 12, 17, 0, "tauijab");
    global_dpd_->contract444(&Wmnij, &tauijab, &newtijab, 1, 1, 1, 1);
    global_dpd_->buf4_close(&tauijab);
    global_dpd_->buf4_close(&Wmnij);
    global_dpd_->buf4_close(&newtijab);

    global_dpd_->buf4_init(&newtIjAb, PSIF_CC_HBAR, 0, 22, 28, 22, 28, 0, "WAbIj residual");
    global_dpd_->buf4_init(&WMnIj, PSIF_CC_HBAR, 0, 22, 22, 22, 22, 0, "WMnIj");
    global_dpd_->buf4_init(&tauIjAb, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tauIjAb");
    global_dpd_->contract444(&WMnIj, &tauIjAb, &newtIjAb, 1, 1, 1, 1);
    global_dpd_->buf4_close(&tauIjAb);
    global_dpd_->buf4_close(&WMnIj);
    global_dpd_->buf4_close(&newtIjAb);
  }
}

void BT2(void)
{
  int h;
  dpdbuf4 newtIJAB, newtijab, newtIjAb;
  dpdbuf4 B_anti, B;
  dpdbuf4 tauIJAB, tauijab, tauIjAb;
  dpdbuf4 Z1,Z2;

  if(params.ref == 0) { /** RHF **/
    global_dpd_->buf4_init(&newtIjAb, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "WAbIj residual");
    global_dpd_->buf4_init(&tauIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
    global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 5, 5, 5, 5, 0, "B <ab|cd>");
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 5, 0, 5, 0, 0, "Z(Ab,Ij)");
    global_dpd_->contract444(&B, &tauIjAb, &Z1, 0, 0, 1, 0);
    global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP0, rspq, 0, 5, "Z(Ij,Ab)");
    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab)");
    global_dpd_->buf4_axpy(&Z2, &newtIjAb, 1);
    global_dpd_->buf4_close(&Z2);
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_close(&B);
    global_dpd_->buf4_close(&tauIjAb);
    global_dpd_->buf4_close(&newtIjAb);
  }
  else if(params.ref == 1) { /** ROHF **/
    global_dpd_->buf4_init(&newtIJAB, PSIF_CC_HBAR, 0, 2, 7, 2, 7, 0, "WABIJ residual");
    global_dpd_->buf4_init(&newtijab, PSIF_CC_HBAR, 0, 2, 7, 2, 7, 0, "Wabij residual");
    global_dpd_->buf4_init(&newtIjAb, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "WAbIj residual");
    global_dpd_->buf4_init(&tauIJAB, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
    global_dpd_->buf4_init(&tauijab, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauijab");
    global_dpd_->buf4_init(&tauIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
    global_dpd_->buf4_init(&B_anti, PSIF_CC_BINTS, 0, 7, 7, 5, 5, 1, "B <ab|cd>");
    global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 5, 5, 5, 5, 0, "B <ab|cd>");
    /* AA and BB terms */
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 7, 2, 7, 2, 0, "Z(ab,ij)");
    global_dpd_->contract444(&B_anti, &tauIJAB, &Z1, 0, 0, 1, 0);
    global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP0, rspq, 2, 7, "Z(ij,ab)");
    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 2, 7, 2, 7, 0, "Z(ij,ab)");
    global_dpd_->buf4_axpy(&Z2, &newtIJAB, 1);
    global_dpd_->buf4_close(&Z2);
    global_dpd_->contract444(&B_anti, &tauijab, &Z1, 0, 0, 1, 0);
    global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP0, rspq, 2, 7, "Z(ij,ab)");
    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 2, 7, 2, 7, 0, "Z(ij,ab)");
    global_dpd_->buf4_axpy(&Z2, &newtijab, 1);
    global_dpd_->buf4_close(&Z2);
    global_dpd_->buf4_close(&Z1);
    /* AB term */
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 5, 0, 5, 0, 0, "Z(Ab,Ij)");
    global_dpd_->contract444(&B, &tauIjAb, &Z1, 0, 0, 1, 0);
    global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP0, rspq, 0, 5, "Z(Ij,Ab)");
    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab)");
    global_dpd_->buf4_axpy(&Z2, &newtIjAb, 1);
    global_dpd_->buf4_close(&Z2);
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_close(&B_anti);
    global_dpd_->buf4_close(&B);
    global_dpd_->buf4_close(&tauIJAB);
    global_dpd_->buf4_close(&tauijab);
    global_dpd_->buf4_close(&tauIjAb);
    global_dpd_->buf4_close(&newtIJAB);
    global_dpd_->buf4_close(&newtijab);
    global_dpd_->buf4_close(&newtIjAb);
  }
  else if(params.ref == 2) { /*** UHF ***/
    global_dpd_->buf4_init(&newtIJAB, PSIF_CC_HBAR, 0, 2, 7, 2, 7, 0, "WABIJ residual");
    global_dpd_->buf4_init(&newtijab, PSIF_CC_HBAR, 0, 12, 17, 12, 17, 0, "Wabij residual");
    global_dpd_->buf4_init(&newtIjAb, PSIF_CC_HBAR, 0, 22, 28, 22, 28, 0, "WAbIj residual");
    global_dpd_->buf4_init(&tauIJAB, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
    global_dpd_->buf4_init(&tauijab, PSIF_CC_TAMPS, 0, 12, 17, 12, 17, 0, "tauijab");
    global_dpd_->buf4_init(&tauIjAb, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tauIjAb");
    global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 7, 7, 5, 5, 1, "B <AB|CD>");
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 7, 2, 7, 2, 0, "Z(AB,IJ)");
    global_dpd_->contract444(&B, &tauIJAB, &Z1, 0, 0, 1, 0);
    global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP0, rspq, 2, 7, "Z(IJ,AB)");
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 2, 7, 2, 7, 0, "Z(IJ,AB)");
    global_dpd_->buf4_axpy(&Z2, &newtIJAB, 1);
    global_dpd_->buf4_close(&Z2);
    global_dpd_->buf4_close(&B);
    global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 17, 17, 15, 15, 1, "B <ab|cd>");
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 17, 12, 17, 12, 0, "Z(ab,ij)");
    global_dpd_->contract444(&B, &tauijab, &Z1, 0, 0, 1, 0);
    global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP0, rspq, 12, 17, "Z(ij,ab)");
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 12, 17, 12, 17, 0, "Z(ij,ab)");
    global_dpd_->buf4_axpy(&Z2, &newtijab, 1);
    global_dpd_->buf4_close(&Z2);
    global_dpd_->buf4_close(&B);
    global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 28, 28, 28, 28, 0, "B <Ab|Cd>");
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 28, 22, 28, 22, 0, "Z(Ab,Ij)");
    global_dpd_->contract444(&B, &tauIjAb, &Z1, 0, 0, 1, 0);
    global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP0, rspq, 22, 28, "Z(Ij,Ab)");
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 22, 28, 22, 28, 0, "Z(Ij,Ab)");
    global_dpd_->buf4_axpy(&Z2, &newtIjAb, 1);
    global_dpd_->buf4_close(&Z2);
    global_dpd_->buf4_close(&B);
    global_dpd_->buf4_close(&tauIJAB);
    global_dpd_->buf4_close(&tauijab);
    global_dpd_->buf4_close(&tauIjAb);
    global_dpd_->buf4_close(&newtIJAB);
    global_dpd_->buf4_close(&newtijab);
    global_dpd_->buf4_close(&newtIjAb);
  }
}

void ZT2(void)
{
  dpdbuf4 ZIJMA, ZIJAM, Zijma, Zijam, ZIjMa, ZIjAm, Z;
  dpdbuf4 newtIJAB, newtijab, newtIjAb, T2;
  dpdfile2 tIA, tia, T1;
  dpdbuf4 t2, X;

  if(params.ref == 0) { /** RHF **/
    global_dpd_->buf4_init(&X, PSIF_CC_TMP0, 0, 5, 0, 5, 0, 0, "X(Ab,Ij)");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->buf4_init(&Z, PSIF_CC_MISC, 0, 10, 0, 10, 0, 0, "ZMbIj");
    global_dpd_->contract244(&T1, &Z, &X, 0, 0, 0, -1, 0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_sort_axpy(&X, PSIF_CC_HBAR, rspq, 0, 5, "WAbIj residual", 1);
    global_dpd_->buf4_sort_axpy(&X, PSIF_CC_HBAR, srqp, 0, 5, "WAbIj residual", 1);
    global_dpd_->buf4_close(&X);
  }
  else if(params.ref == 1) { /** ROHF **/
    global_dpd_->buf4_init(&ZIJMA, PSIF_CC_MISC, 0, 2, 10, 2, 10, 0, "ZIJMA");
    global_dpd_->buf4_init(&ZIJAM, PSIF_CC_MISC, 0, 2, 11, 2, 11, 0, "ZIJAM");
    global_dpd_->buf4_init(&Zijma, PSIF_CC_MISC, 0, 2, 10, 2, 10, 0, "Zijma");
    global_dpd_->buf4_init(&Zijam, PSIF_CC_MISC, 0, 2, 11, 2, 11, 0, "Zijam");
    global_dpd_->buf4_init(&ZIjMa, PSIF_CC_MISC, 0, 0, 10, 0, 10, 0, "ZIjMa");
    global_dpd_->buf4_init(&ZIjAm, PSIF_CC_MISC, 0, 0, 11, 0, 11, 0, "ZIjAm");
    global_dpd_->buf4_init(&newtIJAB, PSIF_CC_HBAR, 0, 2, 5, 2, 7, 0, "WABIJ residual");
    global_dpd_->buf4_init(&newtijab, PSIF_CC_HBAR, 0, 2, 5, 2, 7, 0, "Wabij residual");
    global_dpd_->buf4_init(&newtIjAb, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "WAbIj residual");
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");
    global_dpd_->buf4_init(&t2, PSIF_CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,AB)");
    global_dpd_->contract424(&ZIJAM, &tIA, &t2, 3, 0, 0, 1, 0);
    global_dpd_->contract244(&tIA, &ZIJMA, &t2, 0, 2, 1, -1, 1);
    global_dpd_->buf4_axpy(&t2, &newtIJAB, 1);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_init(&t2, PSIF_CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,AB)");
    global_dpd_->contract424(&Zijam, &tia, &t2, 3, 0, 0, 1, 0);
    global_dpd_->contract244(&tia, &Zijma, &t2, 0, 2, 1, -1, 1);
    global_dpd_->buf4_axpy(&t2, &newtijab, 1);
    global_dpd_->buf4_close(&t2);
    global_dpd_->contract424(&ZIjAm, &tia, &newtIjAb, 3, 0, 0, -1, 1);
    global_dpd_->contract244(&tIA, &ZIjMa, &newtIjAb, 0, 2, 1, -1, 1);
    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);
    global_dpd_->buf4_close(&newtIJAB);
    global_dpd_->buf4_close(&newtijab);
    global_dpd_->buf4_close(&newtIjAb);
    global_dpd_->buf4_close(&ZIJMA);
    global_dpd_->buf4_close(&ZIJAM);
    global_dpd_->buf4_close(&Zijma);
    global_dpd_->buf4_close(&Zijam);
    global_dpd_->buf4_close(&ZIjMa);
    global_dpd_->buf4_close(&ZIjAm);
  }
  else if(params.ref == 2) { /*** UHF ***/
    global_dpd_->buf4_init(&ZIJMA, PSIF_CC_MISC, 0, 2, 20, 2, 20, 0, "ZIJMA");
    global_dpd_->buf4_init(&ZIJAM, PSIF_CC_MISC, 0, 2, 21, 2, 21, 0, "ZIJAM");
    global_dpd_->buf4_init(&Zijma, PSIF_CC_MISC, 0, 12, 30, 12, 30, 0, "Zijma");
    global_dpd_->buf4_init(&Zijam, PSIF_CC_MISC, 0, 12, 31, 12, 31, 0, "Zijam");
    global_dpd_->buf4_init(&ZIjMa, PSIF_CC_MISC, 0, 22, 24, 22, 24, 0, "ZIjMa");
    global_dpd_->buf4_init(&ZIjAm, PSIF_CC_MISC, 0, 22, 26, 22, 26, 0, "ZIjAm");
    global_dpd_->buf4_init(&newtIJAB, PSIF_CC_HBAR, 0, 2, 5, 2, 7, 0, "WABIJ residual");
    global_dpd_->buf4_init(&newtijab, PSIF_CC_HBAR, 0, 12, 15, 12, 17, 0, "Wabij residual");
    global_dpd_->buf4_init(&newtIjAb, PSIF_CC_HBAR, 0, 22, 28, 22, 28, 0, "WAbIj residual");
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");
    global_dpd_->buf4_init(&t2, PSIF_CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,AB)");
    global_dpd_->contract424(&ZIJAM, &tIA, &t2, 3, 0, 0, 1, 0);
    global_dpd_->contract244(&tIA, &ZIJMA, &t2, 0, 2, 1, -1, 1);
    global_dpd_->buf4_axpy(&t2, &newtIJAB, 1);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_init(&t2, PSIF_CC_TMP0, 0, 12, 15, 12, 15, 0, "T (i>j,ab)");
    global_dpd_->contract424(&Zijam, &tia, &t2, 3, 0, 0, 1, 0);
    global_dpd_->contract244(&tia, &Zijma, &t2, 0, 2, 1, -1, 1);
    global_dpd_->buf4_axpy(&t2, &newtijab, 1);
    global_dpd_->buf4_close(&t2);
    global_dpd_->contract424(&ZIjAm, &tia, &newtIjAb, 3, 0, 0, -1, 1);
    global_dpd_->contract244(&tIA, &ZIjMa, &newtIjAb, 0, 2, 1, -1, 1);
    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);
    global_dpd_->buf4_close(&newtIJAB);
    global_dpd_->buf4_close(&newtijab);
    global_dpd_->buf4_close(&newtIjAb);
    global_dpd_->buf4_close(&ZIJMA);
    global_dpd_->buf4_close(&ZIJAM);
    global_dpd_->buf4_close(&Zijma);
    global_dpd_->buf4_close(&Zijam);
    global_dpd_->buf4_close(&ZIjMa);
    global_dpd_->buf4_close(&ZIjAm);
  }
}

void Z_build(void)
{
  dpdbuf4 ZIJMA, Zijma, ZIjMa, ZIjmA, ZIjAm, ZMaIj, ZmAIj, Z;
  dpdbuf4 tauIJAB, tauijab, tauIjAb, tauIjbA, F_anti, F, tau;
  int Gmb, Gij, mb, nrows, ncols;

  timer_on("Z");
  if(params.ref == 0) { /** RHF **/
    /* ZMbIj = <Mb|Ef> * tau(Ij,Ef) */
    /* OOC code added 3/23/05  -TDC */
    global_dpd_->buf4_init(&Z, PSIF_CC_MISC, 0, 10, 0, 10, 0, 0, "ZMbIj");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    global_dpd_->buf4_init(&tau, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
/*     dpd_contract444(&F, &tau, &Z, 0, 0, 1, 0); */
    for(Gmb=0; Gmb < moinfo.nirreps; Gmb++) {
      Gij = Gmb;  /* tau is totally symmetric */
      global_dpd_->buf4_mat_irrep_init(&tau, Gij);
      global_dpd_->buf4_mat_irrep_rd(&tau, Gij);
      global_dpd_->buf4_mat_irrep_init(&Z, Gmb);
      global_dpd_->buf4_mat_irrep_row_init(&F, Gmb);
      for(mb=0; mb < F.params->rowtot[Gmb]; mb++) {
	global_dpd_->buf4_mat_irrep_row_rd(&F, Gmb, mb);
	nrows = tau.params->rowtot[Gij];
	ncols = tau.params->coltot[Gij];
	if(nrows && ncols)
	  C_DGEMV('n',nrows,ncols,1.0,tau.matrix[Gij][0],ncols,F.matrix[Gmb][0],1,
		  0.0,Z.matrix[Gmb][mb],1);
      }
      global_dpd_->buf4_mat_irrep_row_close(&F, Gmb);
      global_dpd_->buf4_mat_irrep_wrt(&Z, Gmb);
      global_dpd_->buf4_mat_irrep_close(&Z, Gmb);
      global_dpd_->buf4_mat_irrep_close(&tau, Gij);
    }
    global_dpd_->buf4_close(&tau);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&Z);
  }
  else if(params.ref == 1) { /** ROHF **/
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
  else if(params.ref == 2) { /*** UHF ***/
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
  timer_off("Z");
}

void FT2(void)
{
  dpdfile2 tIA, tia, t1;
  dpdbuf4 newtIJAB, newtijab, newtIjAb, t2, t2a, t2b;
  dpdbuf4 F_anti, F;
  dpdbuf4 Z, X;
  int Gie, Gij, Gab, nrows, ncols, nlinks, Gi, Ge, Gj, i, I;

  if(params.ref == 0) { /** RHF **/
    /* t(ij,ab) <-- t(j,e) * <ie|ab> + t(i,e) * <je|ba> */
    /* OOC code added 3/23/05, TDC */
    global_dpd_->buf4_init(&X, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(Ij,Ab)");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    global_dpd_->file2_init(&t1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_mat_init(&t1);
    global_dpd_->file2_mat_rd(&t1);
    for(Gie=0; Gie < moinfo.nirreps; Gie++) {
      Gab = Gie; /* F is totally symmetric */
      Gij = Gab; /* T2 is totally symmetric */
      global_dpd_->buf4_mat_irrep_init(&X, Gij);
      ncols = F.params->coltot[Gie];
      for(Gi=0; Gi < moinfo.nirreps; Gi++) {
	Gj = Ge = Gi^Gie; /* T1 is totally symmetric */
	nlinks = moinfo.virtpi[Ge];
	nrows = moinfo.occpi[Gj];
	global_dpd_->buf4_mat_irrep_init_block(&F, Gie, nlinks);
	for(i=0; i < moinfo.occpi[Gi]; i++) {
	  I = F.params->poff[Gi] + i;
	  global_dpd_->buf4_mat_irrep_rd_block(&F, Gie, F.row_offset[Gie][I], nlinks);
	  if(nrows && ncols && nlinks)
	    C_DGEMM('n','n',nrows,ncols,nlinks,1.0,t1.matrix[Gj][0],nlinks,F.matrix[Gie][0],ncols,
		    0.0,X.matrix[Gij][X.row_offset[Gij][I]],ncols);
	}
	global_dpd_->buf4_mat_irrep_close_block(&F, Gie, nlinks);
      }
      global_dpd_->buf4_mat_irrep_wrt(&X, Gij);
      global_dpd_->buf4_mat_irrep_close(&X, Gij);
    }
    global_dpd_->file2_mat_close(&t1);
    global_dpd_->file2_close(&t1);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_init(&t2, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "WAbIj residual");
    global_dpd_->buf4_axpy(&X, &t2, 1);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_sort_axpy(&X, PSIF_CC_HBAR, qpsr, 0, 5, "WAbIj residual", 1);
    global_dpd_->buf4_close(&X);
  }
  else if(params.ref == 1) { /** ROHF **/
    global_dpd_->buf4_init(&newtIJAB, PSIF_CC_HBAR, 0, 0, 7, 2, 7, 0, "WABIJ residual");
    global_dpd_->buf4_init(&newtijab, PSIF_CC_HBAR, 0, 0, 7, 2, 7, 0, "Wabij residual");
    global_dpd_->buf4_init(&newtIjAb, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "WAbIj residual");
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");
    /*** AA ***/
    global_dpd_->buf4_init(&F_anti, PSIF_CC_FINTS, 0, 10, 7, 10, 5, 1, "F <ia|bc>");
    global_dpd_->buf4_init(&t2, PSIF_CC_TMP0, 0, 0, 7, 0, 7, 0, "T (IJ,A>B)");
    global_dpd_->contract424(&F_anti, &tIA, &t2, 1, 1, 1, 1, 0);
    global_dpd_->buf4_sort(&t2, PSIF_CC_TMP0, qprs, 0, 7, "T (JI,A>B)");
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_init(&t2a, PSIF_CC_TMP0, 0, 0, 7, 0, 7, 0, "T (IJ,A>B)");
    global_dpd_->buf4_init(&t2b, PSIF_CC_TMP0, 0, 0, 7, 0, 7, 0, "T (JI,A>B)");
    global_dpd_->buf4_axpy(&t2b, &t2a, -1);
    global_dpd_->buf4_axpy(&t2a, &newtIJAB, 1);
    global_dpd_->buf4_close(&t2b);
    global_dpd_->buf4_close(&t2a);
    global_dpd_->buf4_close(&F_anti);
    /*** BB ***/
    global_dpd_->buf4_init(&F_anti, PSIF_CC_FINTS, 0, 10, 7, 10, 5, 1, "F <ia|bc>");
    global_dpd_->buf4_init(&t2, PSIF_CC_TMP0, 0, 0, 7, 0, 7, 0, "T (IJ,A>B)");
    global_dpd_->contract424(&F_anti, &tia, &t2, 1, 1, 1, 1, 0);
    global_dpd_->buf4_sort(&t2, PSIF_CC_TMP0, qprs, 0, 7, "T (JI,A>B)");
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_init(&t2a, PSIF_CC_TMP0, 0, 0, 7, 0, 7, 0, "T (IJ,A>B)");
    global_dpd_->buf4_init(&t2b, PSIF_CC_TMP0, 0, 0, 7, 0, 7, 0, "T (JI,A>B)");
    global_dpd_->buf4_axpy(&t2b, &t2a, -1);
    global_dpd_->buf4_axpy(&t2a, &newtijab, 1);
    global_dpd_->buf4_close(&t2b);
    global_dpd_->buf4_close(&t2a);
    global_dpd_->buf4_close(&F_anti);
    /*** AB ***/
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    global_dpd_->contract424(&F, &tia, &newtIjAb, 1, 1, 1, 1, 1);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 11, 5, 11, 5, 0, "F <ai|bc>");
    global_dpd_->contract244(&tIA, &F, &newtIjAb, 1, 0, 0, 1, 1);
    global_dpd_->buf4_close(&F);
    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);
    global_dpd_->buf4_close(&newtIJAB);
    global_dpd_->buf4_close(&newtijab);
    global_dpd_->buf4_close(&newtIjAb);
  }
  else if(params.ref == 2) { /*** UHF ***/
    global_dpd_->buf4_init(&newtIJAB, PSIF_CC_HBAR, 0, 0, 7, 2, 7, 0, "WABIJ residual");
    global_dpd_->buf4_init(&newtijab, PSIF_CC_HBAR, 0, 10, 17, 12, 17, 0, "Wabij residual");
    global_dpd_->buf4_init(&newtIjAb, PSIF_CC_HBAR, 0, 22, 28, 22, 28, 0, "WAbIj residual");
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");
    /*** AA ***/
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 20, 7, 20, 5, 1, "F <IA|BC>");
    global_dpd_->buf4_init(&t2, PSIF_CC_TMP0, 0, 0, 7, 0, 7, 0, "T (IJ,A>B)");
    global_dpd_->contract424(&F, &tIA, &t2, 1, 1, 1, 1, 0);
    global_dpd_->buf4_sort(&t2, PSIF_CC_TMP0, qprs, 0, 7, "T (JI,A>B)");
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_init(&t2a, PSIF_CC_TMP0, 0, 0, 7, 0, 7, 0, "T (IJ,A>B)");
    global_dpd_->buf4_init(&t2b, PSIF_CC_TMP0, 0, 0, 7, 0, 7, 0, "T (JI,A>B)");
    global_dpd_->buf4_axpy(&t2b, &t2a, -1);
    global_dpd_->buf4_axpy(&t2a, &newtIJAB, 1);
    global_dpd_->buf4_close(&t2b);
    global_dpd_->buf4_close(&t2a);
    global_dpd_->buf4_close(&F);
    /*** BB ***/
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 30, 17, 30, 15, 1, "F <ia|bc>");
    global_dpd_->buf4_init(&t2, PSIF_CC_TMP0, 0, 10, 17, 10, 17, 0, "T (ij,a>b)");
    global_dpd_->contract424(&F, &tia, &t2, 1, 1, 1, 1, 0);
    global_dpd_->buf4_sort(&t2, PSIF_CC_TMP0, qprs, 10, 17, "T (ji,a>b)");
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_init(&t2a, PSIF_CC_TMP0, 0, 10, 17, 10, 17, 0, "T (ij,a>b)");
    global_dpd_->buf4_init(&t2b, PSIF_CC_TMP0, 0, 10, 17, 10, 17, 0, "T (ji,a>b)");
    global_dpd_->buf4_axpy(&t2b, &t2a, -1);
    global_dpd_->buf4_axpy(&t2a, &newtijab, 1);
    global_dpd_->buf4_close(&t2b);
    global_dpd_->buf4_close(&t2a);
    global_dpd_->buf4_close(&F);
    /*** AB ***/
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
    global_dpd_->contract424(&F, &tia, &newtIjAb, 1, 1, 1, 1, 1);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 28, 26, 28, 26, 0, "F <Ab|Ci>");
    global_dpd_->contract244(&tIA, &F, &newtIjAb, 1, 2, 0, 1, 1);
    global_dpd_->buf4_close(&F);
    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);
    global_dpd_->buf4_close(&newtIJAB);
    global_dpd_->buf4_close(&newtijab);
    global_dpd_->buf4_close(&newtIjAb);
  }
}

void ET2(void)
{
  dpdfile2 tIA, tia;
  dpdbuf4 newtIJAB, newtijab, newtIjAb;
  dpdbuf4 E, t2, t2a, t2b;

  if(params.ref == 0) { /** RHF **/
    global_dpd_->buf4_init(&newtIjAb, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "WAbIj residual");
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    global_dpd_->contract424(&E, &tIA, &newtIjAb, 1, 0, 0, -1, 1);
    global_dpd_->buf4_close(&E);
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 10, 0, 10, 0, 0, "E <ia|jk>");
    global_dpd_->contract244(&tIA, &E, &newtIjAb, 0, 0, 1, -1, 1);
    global_dpd_->buf4_close(&E);
    global_dpd_->file2_close(&tIA);
    global_dpd_->buf4_close(&newtIjAb);
  }
  else if(params.ref == 1) { /** ROHF **/
    global_dpd_->buf4_init(&newtIJAB, PSIF_CC_HBAR, 0, 2, 5, 2, 7, 0, "WABIJ residual");
    global_dpd_->buf4_init(&newtijab, PSIF_CC_HBAR, 0, 2, 5, 2, 7, 0, "Wabij residual");
    global_dpd_->buf4_init(&newtIjAb, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "WAbIj residual");
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");
    /*** AA ***/
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 11, 2, 11, 0, 1, "E <ai|jk>");
    global_dpd_->buf4_init(&t2, PSIF_CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,AB)");
    global_dpd_->contract424(&E, &tIA, &t2, 1, 0, 0, -1, 0);
    global_dpd_->buf4_sort(&t2, PSIF_CC_TMP0, pqsr, 2, 5, "T (I>J,BA)");
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_init(&t2a, PSIF_CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,AB)");
    global_dpd_->buf4_init(&t2b, PSIF_CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,BA)");
    global_dpd_->buf4_axpy(&t2b, &t2a, -1);
    global_dpd_->buf4_close(&t2b);
    global_dpd_->buf4_axpy(&t2a, &newtIJAB, 1);
    global_dpd_->buf4_close(&t2a);
    global_dpd_->buf4_close(&E);
    /*** BB ***/
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 11, 2, 11, 0, 1, "E <ai|jk>");
    global_dpd_->buf4_init(&t2, PSIF_CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,AB)");
    global_dpd_->contract424(&E, &tia, &t2, 1, 0, 0, -1, 0);
    global_dpd_->buf4_sort(&t2, PSIF_CC_TMP0, pqsr, 2, 5, "T (I>J,BA)");
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_init(&t2a, PSIF_CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,AB)");
    global_dpd_->buf4_init(&t2b, PSIF_CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,BA)");
    global_dpd_->buf4_axpy(&t2b, &t2a, -1);
    global_dpd_->buf4_close(&t2b);
    global_dpd_->buf4_axpy(&t2a, &newtijab, 1);
    global_dpd_->buf4_close(&t2a);
    global_dpd_->buf4_close(&E);
    /*** AB ***/
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    global_dpd_->contract424(&E, &tia, &newtIjAb, 1, 0, 0, -1, 1);
    global_dpd_->buf4_close(&E);
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 10, 0, 10, 0, 0, "E <ia|jk>");
    global_dpd_->contract244(&tIA, &E, &newtIjAb, 0, 0, 1, -1, 1);
    global_dpd_->buf4_close(&E);
    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);
    global_dpd_->buf4_close(&newtIJAB);
    global_dpd_->buf4_close(&newtijab);
    global_dpd_->buf4_close(&newtIjAb);
  }
  else if(params.ref == 2) { /*** UHF ***/
    global_dpd_->buf4_init(&newtIJAB, PSIF_CC_HBAR, 0, 2, 5, 2, 7, 0, "WABIJ residual");
    global_dpd_->buf4_init(&newtijab, PSIF_CC_HBAR, 0, 12, 15, 12, 17, 0, "Wabij residual");
    global_dpd_->buf4_init(&newtIjAb, PSIF_CC_HBAR, 0, 22, 28, 22, 28, 0, "WAbIj residual");
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");
    /*** AA ***/
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 21, 2, 21, 0, 1, "E <AI|JK>");
    global_dpd_->buf4_init(&t2, PSIF_CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,AB)");
    global_dpd_->contract424(&E, &tIA, &t2, 1, 0, 0, -1, 0);
    global_dpd_->buf4_sort(&t2, PSIF_CC_TMP0, pqsr, 2, 5, "T (I>J,BA)");
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_init(&t2a, PSIF_CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,AB)");
    global_dpd_->buf4_init(&t2b, PSIF_CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,BA)");
    global_dpd_->buf4_axpy(&t2b, &t2a, -1);
    global_dpd_->buf4_close(&t2b);
    global_dpd_->buf4_axpy(&t2a, &newtIJAB, 1);
    global_dpd_->buf4_close(&t2a);
    global_dpd_->buf4_close(&E);
    /*** BB ***/
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 31, 12, 31, 10, 1, "E <ai|jk>");
    global_dpd_->buf4_init(&t2, PSIF_CC_TMP0, 0, 12, 15, 12, 15, 0, "T (i>j,ab)");
    global_dpd_->contract424(&E, &tia, &t2, 1, 0, 0, -1, 0);
    global_dpd_->buf4_sort(&t2, PSIF_CC_TMP0, pqsr, 12, 15, "T (i>j,ba)");
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_init(&t2a, PSIF_CC_TMP0, 0, 12, 15, 12, 15, 0, "T (i>j,ab)");
    global_dpd_->buf4_init(&t2b, PSIF_CC_TMP0, 0, 12, 15, 12, 15, 0, "T (i>j,ba)");
    global_dpd_->buf4_axpy(&t2b, &t2a, -1);
    global_dpd_->buf4_close(&t2b);
    global_dpd_->buf4_axpy(&t2a, &newtijab, 1);
    global_dpd_->buf4_close(&t2a);
    global_dpd_->buf4_close(&E);
    /*** AB ***/
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 22, 26, 22, 26, 0, "E <Ij|Ak>");
    global_dpd_->contract424(&E, &tia, &newtIjAb, 3, 0, 0, -1, 1);
    global_dpd_->buf4_close(&E);
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 24, 22, 24, 22, 0, "E <Ia|Jk>");
    global_dpd_->contract244(&tIA, &E, &newtIjAb, 0, 0, 1, -1, 1);
    global_dpd_->buf4_close(&E);
    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);
    global_dpd_->buf4_close(&newtIJAB);
    global_dpd_->buf4_close(&newtijab);
    global_dpd_->buf4_close(&newtIjAb);
  }
}


void WmbejT2(void)
{
  dpdbuf4 T2new, T2, W, T2B, W1, W2, Z;

  if(params.ref == 0) { /** RHF **/
    /*** AB ***/
    /* 2 W(ME,jb) + W(Me,Jb) */
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbeJ");
    global_dpd_->buf4_copy(&W, PSIF_CC_TMP0, "2 W(ME,jb) + W(Me,Jb)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W1, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "2 W(ME,jb) + W(Me,Jb)");
    global_dpd_->buf4_init(&W2, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbEj");
    global_dpd_->buf4_axpy(&W2, &W1, 2);
    global_dpd_->buf4_close(&W2);
    global_dpd_->buf4_close(&W1);
    /* T2(Ib,mE) * W(mE,jA) --> Z(Ib,jA) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z (Ib,jA)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIbjA");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbeJ");
    timer_on("WmbejT2 444");
    global_dpd_->contract444(&T2, &W, &T2new, 0, 1, 1, 0);
    timer_off("WmbejT2 444");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T2);
    /* T2(Ib,jA) --> T2(IA,jb) (part III) */
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, psrq, 10, 10, "T2 (IA,jb) 3");
    global_dpd_->buf4_close(&T2new);
    /* 1/2 [ (2 T2(IA,me) - T2(IE,ma)) * (2 W(ME,jb) + W(Me,Jb)] --> T2(IA,jb) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "T2 (IA,jb) 1");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "2 tIAjb - tIBja");
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "2 W(ME,jb) + W(Me,Jb)");
    timer_on("WmbejT2 444");
    global_dpd_->contract444(&T2, &W, &T2new, 0, 1, 0.5, 0);
    timer_off("WmbejT2 444");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T2);
    /* 1/2 Z(Ib,jA) + T2(IA,jb) --> T2(IA,jb) (Part I) */
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z (Ib,jA)");
    global_dpd_->buf4_axpy(&Z, &T2new, 0.5);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&T2new);
    /* T2(IA,jb) (I) + T2(IA,jb) (III) --> T2(IA,jb) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "T2 (IA,jb) 1");
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "T2 (IA,jb) 3");
    global_dpd_->buf4_axpy(&T2, &T2new, 1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, prqs, 0, 5, "T2 (Ij,Ab) (1+3)");
    global_dpd_->buf4_close(&T2new);
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "T2 (Ij,Ab) (1+3)");
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, qpsr, 0, 5, "T2 (Ij,Ab) (2+4)");
    global_dpd_->buf4_close(&T2new);
    /* T2(Ij,Ab) <--- I + II + III + IV */
    global_dpd_->buf4_init(&T2new, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "WAbIj residual");
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "T2 (Ij,Ab) (1+3)");
    global_dpd_->buf4_axpy(&T2, &T2new, 1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "T2 (Ij,Ab) (2+4)");
    global_dpd_->buf4_axpy(&T2, &T2new, 1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&T2new);
  }
  else if(params.ref == 1) { /** ROHF **/
    /*** AA ***/
    /* T2(IA,ME) * W(ME,JB) --> T2(IA,JB) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "T2 (IA,JB)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMBEJ");
    timer_on("WmbejT2 444");
    global_dpd_->contract444(&T2, &W, &T2new, 0, 1, 1, 0);
    timer_off("WmbejT2 444");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T2);
    /* T2(IA,me) * W(me,JB) --> T2(IA,JB) */
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WmBeJ");
    timer_on("WmbejT2 444");
    global_dpd_->contract444(&T2, &W, &T2new, 0, 1, 1, 1);
    timer_off("WmbejT2 444");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T2);
    /* T2(IA,JB) --> T2(IJ,AB) */
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, prqs, 0, 5, "X(0,5) 1");
    global_dpd_->buf4_close(&T2new);
    /* P(IJ) P(AB) T2(IA,JB) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 1");
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, qprs, 0, 5, "X(0,5) 2");
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, pqsr, 0, 5, "X(0,5) 3");
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, qpsr, 0, 5, "X(0,5) 4");
    /* T2(IA,JB) - T2(JA,IB) - T2(IB,JA) + T2(JB,IA) --> T2(IA,JB) */
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 2");
    global_dpd_->buf4_axpy(&T2, &T2new, -1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 3");
    global_dpd_->buf4_axpy(&T2, &T2new, -1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 4");
    global_dpd_->buf4_axpy(&T2, &T2new, +1);
    global_dpd_->buf4_close(&T2);
    /* T2(IJ,AB) --> New T2(IJ,AB) */
    global_dpd_->buf4_init(&T2, PSIF_CC_HBAR, 0, 0, 5, 2, 7, 0, "WABIJ residual");
    global_dpd_->buf4_axpy(&T2new, &T2, 1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&T2new);
    /*** BB ***/
    /* T2(ia,me) * W(me,jb) --> T2(ia,jb) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "T2 (ia,jb)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "Wmbej");
    timer_on("WmbejT2 444");
    global_dpd_->contract444(&T2, &W, &T2new, 0, 1, 1, 0);
    timer_off("WmbejT2 444");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T2);
    /* T2(ia,ME) * W(ME,jb) --> T2(ia,jb) */
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiaJB");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbEj");
    timer_on("WmbejT2 444");
    global_dpd_->contract444(&T2, &W, &T2new, 0, 1, 1, 1);
    timer_off("WmbejT2 444");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T2);
    /* T2(ia,jb) --> T2(ij,ab) */
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, prqs, 0, 5, "X(0,5) 1");
    global_dpd_->buf4_close(&T2new);
    /* P(ij) P(ab) T2(ia,jb) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 1");
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, qprs, 0, 5, "X(0,5) 2");
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, pqsr, 0, 5, "X(0,5) 3");
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, qpsr, 0, 5, "X(0,5) 4");
    /* T2(ij,ab) - T2(ji,ab) - T2(ij,ba) + T2(ji,ba) --> T2(ij,ab) */
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 2");
    global_dpd_->buf4_axpy(&T2, &T2new, -1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 3");
    global_dpd_->buf4_axpy(&T2, &T2new, -1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 4");
    global_dpd_->buf4_axpy(&T2, &T2new, +1);
    global_dpd_->buf4_close(&T2);
    /* T2(ij,ab) --> New T2(ij,ab) */
    global_dpd_->buf4_init(&T2, PSIF_CC_HBAR, 0, 0, 5, 2, 7, 0, "Wabij residual");
    global_dpd_->buf4_axpy(&T2new, &T2, 1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&T2new);

    /*** AB ***/
    /* T2(IA,ME) * W(ME,jb) --> T2(IA,jb) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "T2 (IA,jb)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbEj");
    timer_on("WmbejT2 444");
    global_dpd_->contract444(&T2, &W, &T2new, 0, 1, 1, 0);
    timer_off("WmbejT2 444");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T2);
    /* T2(IA,me) * W(me,jb) --> T2(IA,jb) */
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "Wmbej");
    timer_on("WmbejT2 444");
    global_dpd_->contract444(&T2, &W, &T2new, 0, 1, 1, 1);
    timer_off("WmbejT2 444");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T2);
    /* W(ME,IA) * T2(jb,ME) --> T2(IA,jb) */
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiaJB");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMBEJ");
    timer_on("WmbejT2 444");
    global_dpd_->contract444(&W, &T2, &T2new, 1, 0, 1, 1);
    timer_off("WmbejT2 444");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T2);
    /* W(me,IA) * T2(jb,me) --> T2(IA,jb) */
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WmBeJ");
    timer_on("WmbejT2 444");
    global_dpd_->contract444(&W, &T2, &T2new, 1, 0, 1, 1);
    timer_off("WmbejT2 444");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T2);
    /* T2(IA,jb) --> T2(Ij,Ab) (part 1) */
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, prqs, 0, 5, "T2 (Ij,Ab) 1");
    global_dpd_->buf4_close(&T2new);
    /* T2(Ib,mE) * W(mE,jA) --> T2(Ib,jA) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "T2 (Ib,jA)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIbjA");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WmBEj");
    timer_on("WmbejT2 444");
    global_dpd_->contract444(&T2, &W, &T2new, 0, 1, 1, 0);
    timer_off("WmbejT2 444");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T2);
    /* W(Me,Ib) * T2(jA,Me) --> T2(Ib,jA) */
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tjAIb");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbeJ");
    timer_on("WmbejT2 444");
    global_dpd_->contract444(&W, &T2, &T2new, 1, 0, 1, 1);
    timer_off("WmbejT2 444");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T2);
    /* T2(Ib,jA) --> T2(Ij,Ab) (part 2) */
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, prsq, 0, 5, "T2 (Ij,Ab) 2");
    global_dpd_->buf4_close(&T2new);
    /* T2(Ij,Ab) (part 1) + T2(Ij,Ab) (part 2) --> New T2(Ij,Ab) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "WAbIj residual");
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "T2 (Ij,Ab) 1");
    global_dpd_->buf4_axpy(&T2, &T2new, 1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "T2 (Ij,Ab) 2");
    global_dpd_->buf4_axpy(&T2, &T2new, 1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&T2new);
  } /*** ROHF ***/
  else if(params.ref == 2) { /*** UHF ***/
    /*** AA ***/
    /* T2(IA,ME) * W(ME,JB) --> T2(IA,JB) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 20, 20, 20, 20, 0, "T2 (IA,JB)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 20, 20, 20, 20, 0, "tIAJB");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 20, 20, 20, 20, 0, "WMBEJ");
    timer_on("WmbejT2 444");
    global_dpd_->contract444(&T2, &W, &T2new, 0, 1, 1, 0);
    timer_off("WmbejT2 444");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T2);
    /* T2(IA,me) * W(me,JB) --> T2(IA,JB) */
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 30, 20, 30, 20, 0, "WmBeJ");
    timer_on("WmbejT2 444");
    global_dpd_->contract444(&T2, &W, &T2new, 0, 1, 1, 1);
    timer_off("WmbejT2 444");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T2);
    /* T2(IA,JB) --> T2(IJ,AB) */
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, prqs, 0, 5, "X(0,5) 1");
    global_dpd_->buf4_close(&T2new);
    /* P(IJ) P(AB) T2(IA,JB) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 1");
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, qprs, 0, 5, "X(0,5) 2");
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, pqsr, 0, 5, "X(0,5) 3");
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, qpsr, 0, 5, "X(0,5) 4");
    /* T2(IA,JB) - T2(JA,IB) - T2(IB,JA) + T2(JB,IA) --> T2(IA,JB) */
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 2");
    global_dpd_->buf4_axpy(&T2, &T2new, -1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 3");
    global_dpd_->buf4_axpy(&T2, &T2new, -1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 4");
    global_dpd_->buf4_axpy(&T2, &T2new, +1);
    global_dpd_->buf4_close(&T2);
    /* T2(IJ,AB) --> New T2(IJ,AB) */
    global_dpd_->buf4_init(&T2, PSIF_CC_HBAR, 0, 0, 5, 2, 7, 0, "WABIJ residual");
    global_dpd_->buf4_axpy(&T2new, &T2, 1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&T2new);

    /*** BB ***/
    /* T2(ia,me) * W(me,jb) --> T2(ia,jb) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 30, 30, 30, 30, 0, "T2 (ia,jb)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 30, 30, 30, 30, 0, "tiajb");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 30, 30, 30, 30, 0, "Wmbej");
    timer_on("WmbejT2 444");
    global_dpd_->contract444(&T2, &W, &T2new, 0, 1, 1, 0);
    timer_off("WmbejT2 444");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T2);
    /* T2(ia,ME) * W(ME,jb) --> T2(ia,jb) */
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 30, 20, 30, 20, 0, "tiaJB");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 20, 30, 20, 30, 0, "WMbEj");
    timer_on("WmbejT2 444");
    global_dpd_->contract444(&T2, &W, &T2new, 0, 1, 1, 1);
    timer_off("WmbejT2 444");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T2);
    /* T2(ia,jb) --> T2(ij,ab) */
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, prqs, 10, 15, "X(10,15) 1");
    global_dpd_->buf4_close(&T2new);
    /* P(ij) P(ab) T2(ia,jb) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 10, 15, 10, 15, 0, "X(10,15) 1");
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, qprs, 10, 15, "X(10,15) 2");
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, pqsr, 10, 15, "X(10,15) 3");
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, qpsr, 10, 15, "X(10,15) 4");
    /* T2(ij,ab) - T2(ji,ab) - T2(ij,ba) + T2(ji,ba) --> T2(ij,ab) */
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 10, 15, 10, 15, 0, "X(10,15) 2");
    global_dpd_->buf4_axpy(&T2, &T2new, -1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 10, 15, 10, 15, 0, "X(10,15) 3");
    global_dpd_->buf4_axpy(&T2, &T2new, -1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 10, 15, 10, 15, 0, "X(10,15) 4");
    global_dpd_->buf4_axpy(&T2, &T2new, +1);
    global_dpd_->buf4_close(&T2);
    /* T2(ij,ab) --> New T2(ij,ab) */
    global_dpd_->buf4_init(&T2, PSIF_CC_HBAR, 0, 10, 15, 12, 17, 0, "Wabij residual");
    global_dpd_->buf4_axpy(&T2new, &T2, 1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&T2new);

    /*** AB ***/
    /* T2(IA,ME) * W(ME,jb) --> T2(IA,jb) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 20, 30, 20, 30, 0, "T2 (IA,jb)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 20, 20, 20, 20, 0, "tIAJB");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 20, 30, 20, 30, 0, "WMbEj");
    timer_on("WmbejT2 444");
    global_dpd_->contract444(&T2, &W, &T2new, 0, 1, 1, 0);
    timer_off("WmbejT2 444");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T2);
    /* T2(IA,me) * W(me,jb) --> T2(IA,jb) */
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 30, 30, 30, 30, 0, "Wmbej");
    timer_on("WmbejT2 444");
    global_dpd_->contract444(&T2, &W, &T2new, 0, 1, 1, 1);
    timer_off("WmbejT2 444");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T2);
    /* W(ME,IA) * T2(jb,ME) --> T2(IA,jb) */
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 30, 20, 30, 20, 0, "tiaJB");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 20, 20, 20, 20, 0, "WMBEJ");
    timer_on("WmbejT2 444");
    global_dpd_->contract444(&W, &T2, &T2new, 1, 0, 1, 1);
    timer_off("WmbejT2 444");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T2);
    /* W(me,IA) * T2(jb,me) --> T2(IA,jb) */
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 30, 30, 30, 30, 0, "tiajb");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 30, 20, 30, 20, 0, "WmBeJ");
    timer_on("WmbejT2 444");
    global_dpd_->contract444(&W, &T2, &T2new, 1, 0, 1, 1);
    timer_off("WmbejT2 444");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T2);
    /* T2(IA,jb) --> T2(Ij,Ab) (part 1) */
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, prqs, 22, 28, "T2 (Ij,Ab) 1");
    global_dpd_->buf4_close(&T2new);
    /* T2(Ib,mE) * W(mE,jA) --> T2(Ib,jA) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 24, 27, 24, 27, 0, "T2 (Ib,jA)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 24, 27, 24, 27, 0, "tIbjA");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 27, 27, 27, 27, 0, "WmBEj");
    timer_on("WmbejT2 444");
    global_dpd_->contract444(&T2, &W, &T2new, 0, 1, 1, 0);
    timer_off("WmbejT2 444");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T2);
    /* W(Me,Ib) * T2(jA,Me) --> T2(Ib,jA) */
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 27, 24, 27, 24, 0, "tiBJa");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 24, 24, 24, 24, 0, "WMbeJ");
    timer_on("WmbejT2 444");
    global_dpd_->contract444(&W, &T2, &T2new, 1, 0, 1, 1);
    timer_off("WmbejT2 444");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T2);
    /* T2(Ib,jA) --> T2(Ij,Ab) (part 2) */
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, prsq, 22, 28, "T2 (Ij,Ab) 2");
    global_dpd_->buf4_close(&T2new);
    /* T2(Ij,Ab) (part 1) + T2(Ij,Ab) (part 2) --> New T2(Ij,Ab) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_HBAR, 0, 22, 28, 22, 28, 0, "WAbIj residual");
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 22, 28, 22, 28, 0, "T2 (Ij,Ab) 1");
    global_dpd_->buf4_axpy(&T2, &T2new, 1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 22, 28, 22, 28, 0, "T2 (Ij,Ab) 2");
    global_dpd_->buf4_axpy(&T2, &T2new, 1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&T2new);
  } /*** UHF ***/
}

void CT2(void)
{
  dpdfile2 tIA, tia;
  dpdbuf4 Y, C, D, T2new, T2;

  if(params.ref == 0) { /** RHF **/
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    /*** AB ***/
    /* C(mA|jE) * T(I,E) --> Y(mA,jI) */
    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (mA,jI)");
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    global_dpd_->contract424(&C, &tIA, &Y, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&C);
    /* T(m,b) * Y(mA,jI) --> T2(bA,jI) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 5, 0, 5, 0, 0, "X(5,0)");
    global_dpd_->contract244(&tIA, &Y, &T2new, 0, 0, 0, 1, 0);
    global_dpd_->buf4_close(&Y);
    /* T(bA,jI) --> Tnew(Ij,Ab) */
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, srqp, 0, 5, "X(0,5) 1");
    global_dpd_->buf4_close(&T2new);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 1");
    global_dpd_->buf4_init(&T2new, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "WAbIj residual");
    global_dpd_->buf4_axpy(&T2, &T2new, -1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&T2new);
    /* C(Mb|Ie) * T(j,e) --> Y(Mb,Ij) */
    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (Mb,Ij)");
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    global_dpd_->contract424(&C, &tIA, &Y, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&C);
    /* T(M,A) * Y(Mb,Ij) --> T2(Ab,Ij) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 5, 0, 5, 0, 0, "X(5,0)");
    global_dpd_->contract244(&tIA, &Y, &T2new, 0, 0, 0, 1, 0);
    global_dpd_->buf4_close(&Y);
    /* T(Ab,Ij) --> Tnew(Ij,Ab) */
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, rspq, 0, 5, "X(0,5) 1");
    global_dpd_->buf4_close(&T2new);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 1");
    global_dpd_->buf4_init(&T2new, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "WAbIj residual");
    global_dpd_->buf4_axpy(&T2, &T2new, -1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&T2new);
    /* D(Mb,jE) * T(I,E) --> Y(Mb,jI) */
    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Y(Mb,jI)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
    global_dpd_->contract424(&D, &tIA, &Y, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&D);
    /* T(M,A) * Y(Mb,jI) --> T2(Ab,jI) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 5, 0, 5, 0, 0, "X(5,0)");
    global_dpd_->contract244(&tIA, &Y, &T2new, 0, 0, 0, 1, 0);
    global_dpd_->buf4_close(&Y);
    /* T2(Ab,jI) --> Tnew(Ij,Ab) */
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, srpq, 0, 5, "X(0,5) 1");
    global_dpd_->buf4_close(&T2new);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 1");
    global_dpd_->buf4_init(&T2new, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "WAbIj residual");
    global_dpd_->buf4_axpy(&T2, &T2new, -1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&T2new);
    /* D(mA,Ie) * T(j,e) --> Y(mA,Ij) */
    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Y(mA,Ij)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
    global_dpd_->contract424(&D, &tIA, &Y, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&D);
    /* T(m,b) * Y(mA,Ij) --> T2(bA,Ij) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 5, 0, 5, 0, 0, "X(5,0)");
    global_dpd_->contract244(&tIA, &Y, &T2new, 0, 0, 0, 1, 0);
    global_dpd_->buf4_close(&Y);
    /* T2(bA,Ij) --> Tnew(Ij,Ab) */
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, rsqp, 0, 5, "X(0,5) 1");
    global_dpd_->buf4_close(&T2new);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 1");
    global_dpd_->buf4_init(&T2new, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "WAbIj residual");
    global_dpd_->buf4_axpy(&T2, &T2new, -1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&T2new);
    global_dpd_->file2_close(&tIA);
  }
  else if(params.ref == 1) { /** ROHF **/
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");
    /*** AA ***/
    /* C(MB||JE) * T(I,E) --> Y(MB,JI) */
    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (MB,JI)");
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
    global_dpd_->contract424(&C, &tIA, &Y, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&C);
    /* T(M,A) * Y(MB,JI) --> T(AB,JI) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 5, 0, 5, 0, 0, "X(5,0)");
    global_dpd_->contract244(&tIA, &Y, &T2new, 0, 0, 0, 1, 0);
    global_dpd_->buf4_close(&Y);
    /* T(AB,JI) --> T(IJ,AB) */
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, srpq, 0, 5, "X(0,5) 1");
    global_dpd_->buf4_close(&T2new);
    /* P(IJ) P(AB) T2(IJ,AB) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 1");
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, qprs, 0, 5, "X(0,5) 2");
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, pqsr, 0, 5, "X(0,5) 3");
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, qpsr, 0, 5, "X(0,5) 4");
    /* T2(IJ,AB) - T2(JI,AB) - T2(IJ,BA) - T2(JI,BA) */
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 2");
    global_dpd_->buf4_axpy(&T2, &T2new, -1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 3");
    global_dpd_->buf4_axpy(&T2, &T2new, -1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 4");
    global_dpd_->buf4_axpy(&T2, &T2new, 1);
    global_dpd_->buf4_close(&T2);
    /* T2(IJ,AB) --> T2new (IJ,AB) */
    global_dpd_->buf4_init(&T2, PSIF_CC_HBAR, 0, 0, 5, 2, 7, 0, "WABIJ residual");
    global_dpd_->buf4_axpy(&T2new, &T2, 1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&T2new);
    /*** BB ***/
    /* C(mb||je) * T(i,e) --> Y(mb,ji) */
    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (MB,JI)");
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
    global_dpd_->contract424(&C, &tia, &Y, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&C);
    /* T(m,a) * Y(mb,ji) --> T(ab,ji) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 5, 0, 5, 0, 0, "X(5,0)");
    global_dpd_->contract244(&tia, &Y, &T2new, 0, 0, 0, 1, 0);
    global_dpd_->buf4_close(&Y);
    /* T(ab,ji) --> T(ij,ab) */
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, srpq, 0, 5, "X(0,5) 1");
    global_dpd_->buf4_close(&T2new);
    /* P(ij) P(ab) T2(ij,ab) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 1");
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, qprs, 0, 5, "X(0,5) 2");
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, pqsr, 0, 5, "X(0,5) 3");
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, qpsr, 0, 5, "X(0,5) 4");
    /* T2(ij,ab) - T2(ji,ab) - T2(ij,ba) - T2(ji,ba) */
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 2");
    global_dpd_->buf4_axpy(&T2, &T2new, -1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 3");
    global_dpd_->buf4_axpy(&T2, &T2new, -1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 4");
    global_dpd_->buf4_axpy(&T2, &T2new, 1);
    global_dpd_->buf4_close(&T2);
    /* T2(ij,ab) --> T2new (ij,ab) */
    global_dpd_->buf4_init(&T2, PSIF_CC_HBAR, 0, 0, 5, 2, 7, 0, "Wabij residual");
    global_dpd_->buf4_axpy(&T2new, &T2, 1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&T2new);
    /*** AB ***/
    /* C(mA|jE) * T(I,E) --> Y(mA,jI) */
    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (mA,jI)");
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    global_dpd_->contract424(&C, &tIA, &Y, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&C);
    /* T(m,b) * Y(mA,jI) --> T2(bA,jI) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 5, 0, 5, 0, 0, "X(5,0)");
    global_dpd_->contract244(&tia, &Y, &T2new, 0, 0, 0, 1, 0);
    global_dpd_->buf4_close(&Y);
    /* T(bA,jI) --> Tnew(Ij,Ab) */
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, srqp, 0, 5, "X(0,5) 1");
    global_dpd_->buf4_close(&T2new);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 1");
    global_dpd_->buf4_init(&T2new, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "WAbIj residual");
    global_dpd_->buf4_axpy(&T2, &T2new, -1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&T2new);
    /* C(Mb|Ie) * T(j,e) --> Y(Mb,Ij) */
    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (Mb,Ij)");
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    global_dpd_->contract424(&C, &tia, &Y, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&C);
    /* T(M,A) * Y(Mb,Ij) --> T2(Ab,Ij) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 5, 0, 5, 0, 0, "X(5,0)");
    global_dpd_->contract244(&tIA, &Y, &T2new, 0, 0, 0, 1, 0);
    global_dpd_->buf4_close(&Y);
    /* T(Ab,Ij) --> Tnew(Ij,Ab) */
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, rspq, 0, 5, "X(0,5) 1");
    global_dpd_->buf4_close(&T2new);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 1");
    global_dpd_->buf4_init(&T2new, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "WAbIj residual");
    global_dpd_->buf4_axpy(&T2, &T2new, -1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&T2new);
    /* D(Mb,jE) * T(I,E) --> Y(Mb,jI) */
    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Y(Mb,jI)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
    global_dpd_->contract424(&D, &tIA, &Y, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&D);
    /* T(M,A) * Y(Mb,jI) --> T2(Ab,jI) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 5, 0, 5, 0, 0, "X(5,0)");
    global_dpd_->contract244(&tIA, &Y, &T2new, 0, 0, 0, 1, 0);
    global_dpd_->buf4_close(&Y);
    /* T2(Ab,jI) --> Tnew(Ij,Ab) */
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, srpq, 0, 5, "X(0,5) 1");
    global_dpd_->buf4_close(&T2new);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 1");
    global_dpd_->buf4_init(&T2new, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "WAbIj residual");
    global_dpd_->buf4_axpy(&T2, &T2new, -1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&T2new);
    /* D(mA,Ie) * T(j,e) --> Y(mA,Ij) */
    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Y(mA,Ij)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
    global_dpd_->contract424(&D, &tia, &Y, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&D);
    /* T(m,b) * Y(mA,Ij) --> T2(bA,Ij) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 5, 0, 5, 0, 0, "X(5,0)");
    global_dpd_->contract244(&tia, &Y, &T2new, 0, 0, 0, 1, 0);
    global_dpd_->buf4_close(&Y);
    /* T2(bA,Ij) --> Tnew(Ij,Ab) */
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, rsqp, 0, 5, "X(0,5) 1");
    global_dpd_->buf4_close(&T2new);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 1");
    global_dpd_->buf4_init(&T2new, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "WAbIj residual");
    global_dpd_->buf4_axpy(&T2, &T2new, -1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&T2new);
    global_dpd_->file2_close(&tIA); global_dpd_->file2_close(&tia);
  }
  else if(params.ref == 2) { /** UHF **/
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");
    /*** AA ***/
    /* C(MB||JE) * T(I,E) --> Y(MB,JI) */
    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 20, 0, 20, 0, 0, "Y (MB,JI)");
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 20, 20, 20, 20, 0, "C <IA||JB>");
    global_dpd_->contract424(&C, &tIA, &Y, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&C);
    /* T(M,A) * Y(MB,JI) --> T(AB,JI) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 5, 0, 5, 0, 0, "X(5,0)");
    global_dpd_->contract244(&tIA, &Y, &T2new, 0, 0, 0, 1, 0);
    global_dpd_->buf4_close(&Y);
    /* T(AB,JI) --> T(IJ,AB) */
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, srpq, 0, 5, "X(0,5) 1");
    global_dpd_->buf4_close(&T2new);
    /* P(IJ) P(AB) T2(IJ,AB) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 1");
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, qprs, 0, 5, "X(0,5) 2");
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, pqsr, 0, 5, "X(0,5) 3");
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, qpsr, 0, 5, "X(0,5) 4");
    /* T2(IJ,AB) - T2(JI,AB) - T2(IJ,BA) - T2(JI,BA) */
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 2");
    global_dpd_->buf4_axpy(&T2, &T2new, -1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 3");
    global_dpd_->buf4_axpy(&T2, &T2new, -1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 4");
    global_dpd_->buf4_axpy(&T2, &T2new, 1);
    global_dpd_->buf4_close(&T2);
    /* T2(IJ,AB) --> T2new (IJ,AB) */
    global_dpd_->buf4_init(&T2, PSIF_CC_HBAR, 0, 0, 5, 2, 7, 0, "WABIJ residual");
    global_dpd_->buf4_axpy(&T2new, &T2, 1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&T2new);
    /*** BB ***/
    /* C(mb||je) * T(i,e) --> Y(mb,ji) */
    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 30, 10, 30, 10, 0, "Y (mb,ji)");
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 30, 30, 30, 30, 0, "C <ia||jb>");
    global_dpd_->contract424(&C, &tia, &Y, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&C);
    /* T(m,a) * Y(mb,ji) --> T(ab,ji) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 15, 10, 15, 10, 0, "X(15,10)");
    global_dpd_->contract244(&tia, &Y, &T2new, 0, 0, 0, 1, 0);
    global_dpd_->buf4_close(&Y);
    /* T(ab,ji) --> T(ij,ab) */
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, srpq, 10, 15, "X(10,15) 1");
    global_dpd_->buf4_close(&T2new);
    /* P(ij) P(ab) T2(ij,ab) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 10, 15, 10, 15, 0, "X(10,15) 1");
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, qprs, 10, 15, "X(10,15) 2");
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, pqsr, 10, 15, "X(10,15) 3");
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, qpsr, 10, 15, "X(10,15) 4");
    /* T2(ij,ab) - T2(ji,ab) - T2(ij,ba) - T2(ji,ba) */
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 10, 15, 10, 15, 0, "X(10,15) 2");
    global_dpd_->buf4_axpy(&T2, &T2new, -1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 10, 15, 10, 15, 0, "X(10,15) 3");
    global_dpd_->buf4_axpy(&T2, &T2new, -1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 10, 15, 10, 15, 0, "X(10,15) 4");
    global_dpd_->buf4_axpy(&T2, &T2new, 1);
    global_dpd_->buf4_close(&T2);
    /* T2(ij,ab) --> T2new (ij,ab) */
    global_dpd_->buf4_init(&T2, PSIF_CC_HBAR, 0, 10, 15, 12, 17, 0, "Wabij residual");
    global_dpd_->buf4_axpy(&T2new, &T2, 1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&T2new);
    /*** AB ***/
    /* C(mA|jE) * T(I,E) --> Y(mA,jI) */
    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 27, 23, 27, 23, 0, "Y (mA,jI)");
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 27, 27, 27, 27, 0, "C <iA|jB>");
    global_dpd_->contract424(&C, &tIA, &Y, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&C);
    /* T(m,b) * Y(mA,jI) --> T2(bA,jI) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 29, 23, 29, 23, 0, "X(29,23)");
    global_dpd_->contract244(&tia, &Y, &T2new, 0, 0, 0, 1, 0);
    global_dpd_->buf4_close(&Y);
    /* T(bA,jI) --> Tnew(Ij,Ab) */
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, srqp, 22, 28, "X(22,28) 1");
    global_dpd_->buf4_close(&T2new);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 22, 28, 22, 28, 0, "X(22,28) 1");
    global_dpd_->buf4_init(&T2new, PSIF_CC_HBAR, 0, 22, 28, 22, 28, 0, "WAbIj residual");
    global_dpd_->buf4_axpy(&T2, &T2new, -1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&T2new);
    /* C(Mb|Ie) * T(j,e) --> Y(Mb,Ij) */
    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 24, 22, 24, 22, 0, "Y (Mb,Ij)");
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 24, 24, 24, 24, 0, "C <Ia|Jb>");
    global_dpd_->contract424(&C, &tia, &Y, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&C);
    /* T(M,A) * Y(Mb,Ij) --> T2(Ab,Ij) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 28, 22, 28, 22, 0, "X(28,22)");
    global_dpd_->contract244(&tIA, &Y, &T2new, 0, 0, 0, 1, 0);
    global_dpd_->buf4_close(&Y);
    /* T(Ab,Ij) --> Tnew(Ij,Ab) */
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, rspq, 22, 28, "X(22,28) 1");
    global_dpd_->buf4_close(&T2new);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 22, 28, 22, 28, 0, "X(22,28) 1");
    global_dpd_->buf4_init(&T2new, PSIF_CC_HBAR, 0, 22, 28, 22, 28, 0, "WAbIj residual");
    global_dpd_->buf4_axpy(&T2, &T2new, -1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&T2new);
    /* D(Mb,jE) * T(I,E) --> Y(Mb,jI) */
    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 24, 23, 24, 23, 0, "Y(Mb,jI)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 24, 27, 24, 27, 0, "D <Ij|Ab> (Ib,jA)");
    global_dpd_->contract424(&D, &tIA, &Y, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&D);
    /* T(M,A) * Y(Mb,jI) --> T2(Ab,jI) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 28, 23, 28, 23, 0, "X(28,23)");
    global_dpd_->contract244(&tIA, &Y, &T2new, 0, 0, 0, 1, 0);
    global_dpd_->buf4_close(&Y);
    /* T2(Ab,jI) --> Tnew(Ij,Ab) */
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, srpq, 22, 28, "X(22,28) 1");
    global_dpd_->buf4_close(&T2new);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 22, 28, 22, 28, 0, "X(22,28) 1");
    global_dpd_->buf4_init(&T2new, PSIF_CC_HBAR, 0, 22, 28, 22, 28, 0, "WAbIj residual");
    global_dpd_->buf4_axpy(&T2, &T2new, -1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&T2new);
    /* D(mA,Ie) * T(j,e) --> Y(mA,Ij) */
    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 27, 22, 27, 22, 0, "Y(mA,Ij)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 27, 24, 27, 24, 0, "D <iJ|aB> (iB,Ja)");
    global_dpd_->contract424(&D, &tia, &Y, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&D);
    /* T(m,b) * Y(mA,Ij) --> T2(bA,Ij) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 29, 22, 29, 22, 0, "X(29,22)");
    global_dpd_->contract244(&tia, &Y, &T2new, 0, 0, 0, 1, 0);
    global_dpd_->buf4_close(&Y);
    /* T2(bA,Ij) --> Tnew(Ij,Ab) */
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, rsqp, 22, 28, "X(22,28) 1");
    global_dpd_->buf4_close(&T2new);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 22, 28, 22, 28, 0, "X(22,28) 1");
    global_dpd_->buf4_init(&T2new, PSIF_CC_HBAR, 0, 22, 28, 22, 28, 0, "WAbIj residual");
    global_dpd_->buf4_axpy(&T2, &T2new, -1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&T2new);
    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);
  } /*** UHF ***/
}

void Wmnij_for_Wabij(void)
{
  dpdbuf4 A_anti, A;
  dpdbuf4 WMNIJ, Wmnij, WMnIj, W;
  dpdfile2 tIA, tia;
  dpdbuf4 Eijka, Eijka_anti, Eaijk, Eaijk_anti;
  dpdbuf4 D_anti, D, tauIJAB, tauijab, tauIjAb;

  timer_on("Wmnij");
  if(params.ref == 0) { /** RHF **/
    global_dpd_->buf4_init(&A, PSIF_CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
    global_dpd_->buf4_copy(&A, PSIF_CC_HBAR, "WMnIj");
    global_dpd_->buf4_close(&A);
  }
  else if(params.ref == 1) { /** ROHF **/
    global_dpd_->buf4_init(&A_anti, PSIF_CC_AINTS, 0, 2, 2, 0, 0, 1, "A <ij|kl>");
    global_dpd_->buf4_copy(&A_anti, PSIF_CC_HBAR, "WMNIJ");
    global_dpd_->buf4_copy(&A_anti, PSIF_CC_HBAR, "Wmnij");
    global_dpd_->buf4_close(&A_anti);
    global_dpd_->buf4_init(&A, PSIF_CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
    global_dpd_->buf4_copy(&A, PSIF_CC_HBAR, "WMnIj");
    global_dpd_->buf4_close(&A);
  }
  else if(params.ref == 2) { /*** UHF ***/
    global_dpd_->buf4_init(&A, PSIF_CC_AINTS, 0, 2, 2, 0, 0, 1, "A <IJ|KL>");
    global_dpd_->buf4_copy(&A, PSIF_CC_HBAR, "WMNIJ");
    global_dpd_->buf4_close(&A);
    global_dpd_->buf4_init(&A, PSIF_CC_AINTS, 0, 12, 12, 10, 10, 1, "A <ij|kl>");
    global_dpd_->buf4_copy(&A, PSIF_CC_HBAR, "Wmnij");
    global_dpd_->buf4_close(&A);
    global_dpd_->buf4_init(&A, PSIF_CC_AINTS, 0, 22, 22, 22, 22, 0, "A <Ij|Kl>");
    global_dpd_->buf4_copy(&A, PSIF_CC_HBAR, "WMnIj");
    global_dpd_->buf4_close(&A);
  }
  if(params.ref == 0) { /** RHF **/
    global_dpd_->buf4_init(&WMnIj, PSIF_CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj");
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->buf4_init(&Eaijk, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    global_dpd_->contract244(&tIA, &Eaijk, &WMnIj, 1, 0, 1, 1, 1);
    global_dpd_->buf4_close(&Eaijk);
    global_dpd_->buf4_init(&Eijka, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    global_dpd_->contract424(&Eijka, &tIA, &WMnIj, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&Eijka);
    global_dpd_->file2_close(&tIA);
    global_dpd_->buf4_close(&WMnIj);
  }
  else if(params.ref == 1) { /** ROHF **/
    global_dpd_->buf4_init(&WMNIJ, PSIF_CC_HBAR, 0, 2, 0, 2, 2, 0, "WMNIJ");
    global_dpd_->buf4_init(&Wmnij, PSIF_CC_HBAR, 0, 2, 0, 2, 2, 0, "Wmnij");
    global_dpd_->buf4_init(&WMnIj, PSIF_CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj");
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");
    global_dpd_->buf4_init(&Eijka_anti, PSIF_CC_EINTS, 0, 2, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
    global_dpd_->buf4_init(&Eijka, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    global_dpd_->buf4_init(&Eaijk_anti, PSIF_CC_EINTS, 0, 11, 2, 11, 0, 1, "E <ai|jk>");
    global_dpd_->buf4_init(&Eaijk, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 2, 0, 2, 0, 0, "W (MN,IJ)");
    global_dpd_->contract424(&Eijka_anti, &tIA, &W, 3, 1, 0, 1, 0);
    global_dpd_->contract244(&tIA, &Eaijk_anti, &W, 1, 0, 1, 1, 1);
    global_dpd_->buf4_axpy(&W, &WMNIJ, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 2, 0, 2, 0, 0, "W (MN,IJ)");
    global_dpd_->contract424(&Eijka_anti, &tia, &W, 3, 1, 0, 1, 0);
    global_dpd_->contract244(&tia, &Eaijk_anti, &W, 1, 0, 1, 1, 1);
    global_dpd_->buf4_axpy(&W, &Wmnij, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->contract424(&Eijka, &tia, &WMnIj, 3, 1, 0, 1, 1);
    global_dpd_->contract244(&tIA, &Eaijk, &WMnIj, 1, 0, 1, 1, 1);
    global_dpd_->buf4_close(&Eijka_anti);
    global_dpd_->buf4_close(&Eijka);
    global_dpd_->buf4_close(&Eaijk_anti);
    global_dpd_->buf4_close(&Eaijk);
    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);
    global_dpd_->buf4_close(&WMNIJ);
    global_dpd_->buf4_close(&Wmnij);
    global_dpd_->buf4_close(&WMnIj);
  }
  else if(params.ref == 2) { /*** UHF ***/
    global_dpd_->buf4_init(&WMNIJ, PSIF_CC_HBAR, 0, 2, 0, 2, 2, 0, "WMNIJ");
    global_dpd_->buf4_init(&Wmnij, PSIF_CC_HBAR, 0, 12, 10, 12, 12, 0, "Wmnij");
    global_dpd_->buf4_init(&WMnIj, PSIF_CC_HBAR, 0, 22, 22, 22, 22, 0, "WMnIj");
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");
    global_dpd_->buf4_init(&Eijka, PSIF_CC_EINTS, 0, 2, 20, 2, 20, 0, "E <IJ||KA> (I>J,KA)");
    global_dpd_->buf4_init(&Eaijk, PSIF_CC_EINTS, 0, 21, 2, 21, 0, 1, "E <AI|JK>");
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 2, 0, 2, 0, 0, "W (MN,IJ)");
    global_dpd_->contract424(&Eijka, &tIA, &W, 3, 1, 0, 1, 0);
    global_dpd_->contract244(&tIA, &Eaijk, &W, 1, 0, 1, 1, 1);
    global_dpd_->buf4_axpy(&W, &WMNIJ, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Eijka);
    global_dpd_->buf4_close(&Eaijk);
    global_dpd_->buf4_init(&Eijka, PSIF_CC_EINTS, 0, 12, 30, 12, 30, 0, "E <ij||ka> (i>j,ka)");
    global_dpd_->buf4_init(&Eaijk, PSIF_CC_EINTS, 0, 31, 12, 31, 10, 1, "E <ai|jk>");
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 12, 10, 12, 10, 0, "W (mn,ij)");
    global_dpd_->contract424(&Eijka, &tia, &W, 3, 1, 0, 1, 0);
    global_dpd_->contract244(&tia, &Eaijk, &W, 1, 0, 1, 1, 1);
    global_dpd_->buf4_axpy(&W, &Wmnij, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Eijka);
    global_dpd_->buf4_close(&Eaijk);
    global_dpd_->buf4_init(&Eijka, PSIF_CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
    global_dpd_->buf4_init(&Eaijk, PSIF_CC_EINTS, 0, 26, 22, 26, 22, 0, "E <Ai|Jk>");
    global_dpd_->contract424(&Eijka, &tia, &WMnIj, 3, 1, 0, 1, 1);
    global_dpd_->contract244(&tIA, &Eaijk, &WMnIj, 1, 0, 1, 1, 1);
    global_dpd_->buf4_close(&Eijka);
    global_dpd_->buf4_close(&Eaijk);
    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);
    global_dpd_->buf4_close(&WMNIJ);
    global_dpd_->buf4_close(&Wmnij);
    global_dpd_->buf4_close(&WMnIj);
  }
  if(params.ref == 0) { /** RHF **/
    global_dpd_->buf4_init(&WMnIj, PSIF_CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->buf4_init(&tauIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
    global_dpd_->contract444(&D, &tauIjAb, &WMnIj, 0, 0, 1, 1);
    global_dpd_->buf4_close(&tauIjAb);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&WMnIj);
  }
  else if(params.ref == 1) { /** ROHF **/
    global_dpd_->buf4_init(&WMNIJ, PSIF_CC_HBAR, 0, 2, 2, 2, 2, 0, "WMNIJ");
    global_dpd_->buf4_init(&Wmnij, PSIF_CC_HBAR, 0, 2, 2, 2, 2, 0, "Wmnij");
    global_dpd_->buf4_init(&WMnIj, PSIF_CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj");
    global_dpd_->buf4_init(&D_anti, PSIF_CC_DINTS, 0, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->buf4_init(&tauIJAB, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
    global_dpd_->buf4_init(&tauijab, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauijab");
    global_dpd_->buf4_init(&tauIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
    global_dpd_->contract444(&D_anti, &tauIJAB, &WMNIJ, 0, 0, 1, 1);
    global_dpd_->contract444(&D_anti, &tauijab, &Wmnij, 0, 0, 1, 1);
    global_dpd_->contract444(&D, &tauIjAb, &WMnIj, 0, 0, 1, 1);
    global_dpd_->buf4_close(&tauIJAB);
    global_dpd_->buf4_close(&tauijab);
    global_dpd_->buf4_close(&tauIjAb);
    global_dpd_->buf4_close(&D_anti);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&WMNIJ);
    global_dpd_->buf4_close(&Wmnij);
    global_dpd_->buf4_close(&WMnIj);
  }
  else if(params.ref == 2) { /*** UHF ***/
    global_dpd_->buf4_init(&WMNIJ, PSIF_CC_HBAR, 0, 2, 2, 2, 2, 0, "WMNIJ");
    global_dpd_->buf4_init(&Wmnij, PSIF_CC_HBAR, 0, 12, 12, 12, 12, 0, "Wmnij");
    global_dpd_->buf4_init(&WMnIj, PSIF_CC_HBAR, 0, 22, 22, 22, 22, 0, "WMnIj");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 2, 7, 2, 7, 0, "D <IJ||AB> (I>J,A>B)");
    global_dpd_->buf4_init(&tauIJAB, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
    global_dpd_->contract444(&D, &tauIJAB, &WMNIJ, 0, 0, 1, 1);
    global_dpd_->buf4_close(&tauIJAB);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 12, 17, 12, 17, 0, "D <ij||ab> (i>j,a>b)");
    global_dpd_->buf4_init(&tauijab, PSIF_CC_TAMPS, 0, 12, 17, 12, 17, 0, "tauijab");
    global_dpd_->contract444(&D, &tauijab, &Wmnij, 0, 0, 1, 1);
    global_dpd_->buf4_close(&tauijab);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    global_dpd_->buf4_init(&tauIjAb, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tauIjAb");
    global_dpd_->contract444(&D, &tauIjAb, &WMnIj, 0, 0, 1, 1);
    global_dpd_->buf4_close(&tauIjAb);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&WMNIJ);
    global_dpd_->buf4_close(&Wmnij);
    global_dpd_->buf4_close(&WMnIj);
  }
  timer_off("Wmnij");
}

void Wmbej_for_Wabij(void)
{
  dpdbuf4 WMBEJ, Wmbej, WMbEj, WmBeJ, WmBEj, WMbeJ, W;
  dpdbuf4 C, D, E, F, X, tIAjb, tiaJB, t2, Y, Z;
  dpdfile2 tIA, tia;
  int Gmb, mb, Gj, Ge, Gf, nrows, ncols, nlinks;

  timer_on("C->Wmbej");
  /* W(mb,je) <-- <mb||ej> */
  if(params.ref == 0) { /** RHF **/
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    global_dpd_->buf4_scmcopy(&C, PSIF_CC_TMP0, "WMbeJ", -1);
    global_dpd_->buf4_close(&C);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    global_dpd_->buf4_copy(&D, PSIF_CC_TMP0, "WMbEj");
    global_dpd_->buf4_close(&D);
  }
  else if(params.ref == 1) { /*** ROHF ***/
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 11, 10, 11, 0, "C <ia||jb> (ia,bj)");
    global_dpd_->buf4_scmcopy(&C, PSIF_CC_TMP0, "WMBEJ", -1);
    global_dpd_->buf4_scmcopy(&C, PSIF_CC_TMP0, "Wmbej", -1);
    global_dpd_->buf4_close(&C);
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    global_dpd_->buf4_scmcopy(&C, PSIF_CC_TMP0, "WmBEj", -1);
    global_dpd_->buf4_scmcopy(&C, PSIF_CC_TMP0, "WMbeJ", -1);
    global_dpd_->buf4_close(&C);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    global_dpd_->buf4_copy(&D, PSIF_CC_TMP0, "WMbEj");
    global_dpd_->buf4_copy(&D, PSIF_CC_TMP0, "WmBeJ");
    global_dpd_->buf4_close(&D);
  }
  else if(params.ref == 2) { /*** UHF ***/
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 20, 21, 20, 21, 0, "C <IA||JB> (IA,BJ)");
    global_dpd_->buf4_scmcopy(&C, PSIF_CC_TMP0, "WMBEJ", -1);
    global_dpd_->buf4_close(&C);
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 30, 31, 30, 31, 0, "C <ia||jb> (ia,bj)");
    global_dpd_->buf4_scmcopy(&C, PSIF_CC_TMP0, "Wmbej", -1);
    global_dpd_->buf4_close(&C);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 24, 26, 24, 26, 0, "D <Ij|Ab> (Ib,Aj)");
    global_dpd_->buf4_scmcopy(&D, PSIF_CC_TMP0, "WMbEj", 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 27, 25, 27, 25, 0, "D <iJ|aB> (iB,aJ)");
    global_dpd_->buf4_scmcopy(&D, PSIF_CC_TMP0, "WmBeJ", 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 27, 27, 27, 27, 0, "C <iA|jB>");
    global_dpd_->buf4_scmcopy(&C, PSIF_CC_TMP0, "WmBEj", -1);
    global_dpd_->buf4_close(&C);
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 24, 24, 24, 24, 0, "C <Ia|Jb>");
    global_dpd_->buf4_scmcopy(&C, PSIF_CC_TMP0, "WMbeJ", -1);
    global_dpd_->buf4_close(&C);
  }
  timer_off("C->Wmbej");
  timer_on("F->Wmbej");
  if(params.ref == 0) { /** RHF **/
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    global_dpd_->buf4_init(&WMbEj, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "WMbEj");
    global_dpd_->contract424(&F, &tIA, &WMbEj, 3, 1, 0, 1, 1); /* should be OOC-capable in libdpd */
    global_dpd_->buf4_close(&WMbEj);
    global_dpd_->buf4_close(&F);
    /*
    dpd_buf4_init(&F, CC_FINTS, 0, 11, 5, 11, 5, 0, "F <ai|bc>");
    dpd_buf4_init(&Z, CC_TMP0, 0, 11, 11, 11, 11, 0, "Z(bM,eJ)");
    dpd_contract424(&F, &tIA, &Z, 3, 1, 0, -1, 0);
    dpd_buf4_sort(&Z, CC_TMP0, qpsr, 10, 10, "Z(Mb,Je)");
    dpd_buf4_close(&Z);
    dpd_buf4_init(&Z, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(Mb,Je)");
    dpd_buf4_init(&WMbeJ, CC_TMP0, 0, 10, 10, 10, 10, 0, "WMbeJ");
    dpd_buf4_axpy(&Z, &WMbeJ, 1.0);
    dpd_buf4_close(&WMbeJ);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&F);
    */
    /* W(Mb,Je) <-- t(J,F) <Mb|Fe> */
    /* OOC code added to replace above on 3/23/05, TDC */
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "WMbeJ");
    global_dpd_->file2_mat_init(&tIA);
    global_dpd_->file2_mat_rd(&tIA);
    for(Gmb=0; Gmb < moinfo.nirreps; Gmb++) {
      global_dpd_->buf4_mat_irrep_row_init(&W, Gmb);
      global_dpd_->buf4_mat_irrep_row_init(&F, Gmb);
      for(mb=0; mb < F.params->rowtot[Gmb]; mb++) {
	global_dpd_->buf4_mat_irrep_row_rd(&W, Gmb, mb);
	global_dpd_->buf4_mat_irrep_row_rd(&F, Gmb, mb);
	for(Gj=0; Gj < moinfo.nirreps; Gj++) {
	  Gf = Gj;  /* T1 is totally symmetric */
	  Ge = Gmb ^ Gf; /* <mb|fe> is totally symmetric */
	  nrows = moinfo.occpi[Gj];
	  ncols = moinfo.virtpi[Ge];
	  nlinks = moinfo.virtpi[Gf];
	  if(nrows && ncols && nlinks)
	    C_DGEMM('n','n',nrows,ncols,nlinks,-1.0,tIA.matrix[Gj][0],nlinks,
		    &F.matrix[Gmb][0][F.col_offset[Gmb][Gf]],ncols,1.0,
		    &W.matrix[Gmb][0][W.col_offset[Gmb][Gj]],ncols);
	}
	global_dpd_->buf4_mat_irrep_row_wrt(&W, Gmb, mb);
      }
      global_dpd_->buf4_mat_irrep_row_close(&F, Gmb);
      global_dpd_->buf4_mat_irrep_row_close(&W, Gmb);
    }
    global_dpd_->file2_mat_close(&tIA);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&F);
    global_dpd_->file2_close(&tIA);
  }
  else if(params.ref == 1) { /** ROHF **/
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 1, "F <ia|bc>");
    global_dpd_->buf4_init(&WMBEJ, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "WMBEJ");
    global_dpd_->contract424(&F, &tIA, &WMBEJ, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&WMBEJ);
    global_dpd_->buf4_init(&Wmbej, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "Wmbej");
    global_dpd_->contract424(&F, &tia, &Wmbej, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&Wmbej);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    global_dpd_->buf4_init(&WMbEj, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "WMbEj");
    global_dpd_->contract424(&F, &tia, &WMbEj, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&WMbEj);
    global_dpd_->buf4_init(&WmBeJ, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "WmBeJ");
    global_dpd_->contract424(&F, &tIA, &WmBeJ, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&WmBeJ);
    global_dpd_->buf4_init(&WMbeJ, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "WMbeJ");
    global_dpd_->contract244(&tIA, &F, &WMbeJ, 1, 2, 1, -1, 1);
    global_dpd_->buf4_close(&WMbeJ);
    global_dpd_->buf4_init(&WmBEj, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "WmBEj");
    global_dpd_->contract244(&tia, &F, &WmBEj, 1, 2, 1, -1, 1);
    global_dpd_->buf4_close(&WmBEj);
    global_dpd_->buf4_close(&F);
    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);
  }
  else if(params.ref == 2) { /** UHF **/
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 20, 21, 20, 21, 0, "WMBEJ");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 20, 5, 20, 5, 1, "F <IA|BC>");
    global_dpd_->contract424(&F, &tIA, &W, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 30, 31, 30, 31, 0, "Wmbej");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 30, 15, 30, 15, 1, "F <ia|bc>");
    global_dpd_->contract424(&F, &tia, &W, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 24, 26, 24, 26, 0, "WMbEj");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
    global_dpd_->contract424(&F, &tia, &W, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 27, 25, 27, 25, 0, "WmBeJ");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 27, 29, 27, 29, 0, "F <iA|bC>");
    global_dpd_->contract424(&F, &tIA, &W, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 24, 24, 24, 24, 0, "WMbeJ");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
    global_dpd_->contract244(&tIA, &F, &W, 1, 2, 1, -1, 1);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 27, 27, 27, 27, 0, "WmBEj");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 27, 29, 27, 29, 0, "F <iA|bC>");
    global_dpd_->contract244(&tia, &F, &W, 1, 2, 1, -1, 1);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&W);
    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);
  }
  timer_off("F->Wmbej");
  timer_on("E->Wmbej");
  if(params.ref == 0) { /** RHF **/
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    global_dpd_->buf4_init(&WMbEj, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "WMbEj");
    global_dpd_->contract424(&E, &tIA, &WMbEj, 3, 0, 1, -1, 1);
    global_dpd_->buf4_close(&WMbEj);
    global_dpd_->buf4_close(&E);

    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    global_dpd_->buf4_init(&WMbeJ, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "WMbeJ");
    global_dpd_->contract424(&E, &tIA, &WMbeJ, 1, 0, 1, 1, 1);
    global_dpd_->buf4_close(&WMbeJ);
    global_dpd_->buf4_close(&E);
    global_dpd_->file2_close(&tIA);
  }
  else if(params.ref == 1) { /** ROHF **/
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 0, 11, 2, 11, 0, "E <ij||ka> (i>j,ak)");
    global_dpd_->buf4_init(&WMBEJ, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "WMBEJ");
    global_dpd_->contract424(&E, &tIA, &WMBEJ, 1, 0, 1, 1, 1);
    global_dpd_->buf4_close(&WMBEJ);
    global_dpd_->buf4_init(&Wmbej, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "Wmbej");
    global_dpd_->contract424(&E, &tia, &Wmbej, 1, 0, 1, 1, 1);
    global_dpd_->buf4_close(&Wmbej);
    global_dpd_->buf4_close(&E);
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    global_dpd_->buf4_init(&WMbEj, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "WMbEj");
    global_dpd_->contract424(&E, &tia, &WMbEj, 3, 0, 1, -1, 1);
    global_dpd_->buf4_close(&WMbEj);
    global_dpd_->buf4_init(&WmBeJ, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "WmBeJ");
    global_dpd_->contract424(&E, &tIA, &WmBeJ, 3, 0, 1, -1, 1);
    global_dpd_->buf4_close(&WmBeJ);
    global_dpd_->buf4_close(&E);
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    global_dpd_->buf4_init(&WMbeJ, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "WMbeJ");
    global_dpd_->contract424(&E, &tia, &WMbeJ, 1, 0, 1, 1, 1);
    global_dpd_->buf4_close(&WMbeJ);
    global_dpd_->buf4_init(&WmBEj, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "WmBEj");
    global_dpd_->contract424(&E, &tIA, &WmBEj, 1, 0, 1, 1, 1);
    global_dpd_->buf4_close(&WmBEj);
    global_dpd_->buf4_close(&E);
    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);
  }
  else if(params.ref == 2) { /** UHF **/
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 20, 21, 20, 21, 0, "WMBEJ");
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 0, 21, 2, 21, 0, "E <IJ||KA> (I>J,AK)");
    global_dpd_->contract424(&E, &tIA, &W, 1, 0, 1, 1, 1);
    global_dpd_->buf4_close(&E);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 30, 31, 30, 31, 0, "Wmbej");
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 10, 31, 12, 31, 0, "E <ij||ka> (i>j,ak)");
    global_dpd_->contract424(&E, &tia, &W, 1, 0, 1, 1, 1);
    global_dpd_->buf4_close(&E);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 24, 26, 24, 26, 0, "WMbEj");
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 22, 26, 22, 26, 0, "E <Ij|Ak>");
    global_dpd_->contract424(&E, &tia, &W, 1, 0, 1, -1, 1);
    global_dpd_->buf4_close(&E);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 27, 25, 27, 25, 0, "WmBeJ");
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 23, 25, 23, 25, 0, "E <iJ|aK>");
    global_dpd_->contract424(&E, &tIA, &W, 1, 0, 1, -1, 1);
    global_dpd_->buf4_close(&E);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 24, 24, 24, 24, 0, "WMbeJ");
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
    global_dpd_->contract424(&E, &tia, &W, 1, 0, 1, 1, 1);
    global_dpd_->buf4_close(&E);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 27, 27, 27, 27, 0, "WmBEj");
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 23, 27, 23, 27, 0, "E <iJ|kA>");
    global_dpd_->contract424(&E, &tIA, &W, 1, 0, 1, 1, 1);
    global_dpd_->buf4_close(&E);
    global_dpd_->buf4_close(&W);
    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);
  }
  timer_off("E->Wmbej");
  /* Convert to (ME,JB) for remaining terms */
  timer_on("sort Wmbej");
  if(params.ref == 0) { /** RHF **/
    global_dpd_->buf4_init(&WMbEj, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "WMbEj");
    global_dpd_->buf4_sort(&WMbEj, PSIF_CC_HBAR, prsq, 10, 10, "WMbEj");
    global_dpd_->buf4_close(&WMbEj);
    global_dpd_->buf4_init(&WMbeJ, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "WMbeJ");
    global_dpd_->buf4_sort(&WMbeJ, PSIF_CC_HBAR, psrq, 10, 10, "WMbeJ");
    global_dpd_->buf4_close(&WMbeJ);
  }
  else if(params.ref == 1) { /** ROHF **/
    global_dpd_->buf4_init(&WMBEJ, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "WMBEJ");
    global_dpd_->buf4_sort(&WMBEJ, PSIF_CC_HBAR, prsq, 10, 10, "WMBEJ");
    global_dpd_->buf4_close(&WMBEJ);
    global_dpd_->buf4_init(&Wmbej, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "Wmbej");
    global_dpd_->buf4_sort(&Wmbej, PSIF_CC_HBAR, prsq, 10, 10, "Wmbej");
    global_dpd_->buf4_close(&Wmbej);
    global_dpd_->buf4_init(&WMbEj, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "WMbEj");
    global_dpd_->buf4_sort(&WMbEj, PSIF_CC_HBAR, prsq, 10, 10, "WMbEj");
    global_dpd_->buf4_close(&WMbEj);
    global_dpd_->buf4_init(&WmBeJ, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "WmBeJ");
    global_dpd_->buf4_sort(&WmBeJ, PSIF_CC_HBAR, prsq, 10, 10, "WmBeJ");
    global_dpd_->buf4_close(&WmBeJ);
    global_dpd_->buf4_init(&WMbeJ, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "WMbeJ");
    global_dpd_->buf4_sort(&WMbeJ, PSIF_CC_HBAR, psrq, 10, 10, "WMbeJ");
    global_dpd_->buf4_close(&WMbeJ);
    global_dpd_->buf4_init(&WmBEj, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "WmBEj");
    global_dpd_->buf4_sort(&WmBEj, PSIF_CC_HBAR, psrq, 10, 10, "WmBEj");
    global_dpd_->buf4_close(&WmBEj);
  }
  else if(params.ref == 2) { /** UHF **/
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 20, 21, 20, 21, 0, "WMBEJ");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, prsq, 20, 20, "WMBEJ");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 30, 31, 30, 31, 0, "Wmbej");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, prsq, 30, 30, "Wmbej");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 24, 26, 24, 26, 0, "WMbEj");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, prsq, 20, 30, "WMbEj");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 27, 25, 27, 25, 0, "WmBeJ");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, prsq, 30, 20, "WmBeJ");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 24, 24, 24, 24, 0, "WMbeJ");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, psrq, 24, 24, "WMbeJ");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 27, 27, 27, 27, 0, "WmBEj");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, psrq, 27, 27, "WmBEj");
    global_dpd_->buf4_close(&W);
  }
  timer_off("sort Wmbej");
  timer_on("X->Wmbej");
  if(params.ref == 0) { /** RHF **/
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    /*** ABAB ***/
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbEj");
    global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D 2<ij|ab> - <ij|ba> (ia,jb)");
    global_dpd_->contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tjAIb");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    global_dpd_->contract444(&D, &t2, &W, 0, 0, -0.5, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ia,bj)");
    global_dpd_->contract244(&tIA, &D, &Y, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbEj");
    global_dpd_->contract424(&Y, &tIA, &W, 3, 0, 0, -1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y);
    /*** ABBA ***/
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbeJ");
    global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIbjA");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
    global_dpd_->contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    global_dpd_->contract244(&tIA, &D, &Y, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbeJ");
    global_dpd_->contract424(&Y, &tIA, &W, 3, 0, 0, 1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y);
    global_dpd_->file2_close(&tIA);
  }
  else if(params.ref == 1) { /** ROHF **/
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");
    /*** AAAA ***/
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMBEJ");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
    global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
    global_dpd_->contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    global_dpd_->contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij||ab> (ia,bj)");
    global_dpd_->contract244(&tIA, &D, &Y, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMBEJ");
    global_dpd_->contract424(&Y, &tIA, &W, 3, 0, 0, -1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y);
    /*** BBBB ***/
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "Wmbej");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
    global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
    global_dpd_->contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiaJB");
    global_dpd_->contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij||ab> (ia,bj)");
    global_dpd_->contract244(&tia, &D, &Y, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "Wmbej");
    global_dpd_->contract424(&Y, &tia, &W, 3, 0, 0, -1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y);
    /*** ABAB ***/
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbEj");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
    global_dpd_->contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
    global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiaJB");
    global_dpd_->contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ia,bj)");
    global_dpd_->contract244(&tia, &D, &Y, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbEj");
    global_dpd_->contract424(&Y, &tia, &W, 3, 0, 0, -1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y);
    /*** BABA ***/
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WmBeJ");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
    global_dpd_->contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
    global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    global_dpd_->contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ia,bj)");
    global_dpd_->contract244(&tIA, &D, &Y, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WmBeJ");
    global_dpd_->contract424(&Y, &tIA, &W, 3, 0, 0, -1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y);
    /*** ABBA ***/
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbeJ");
    global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIbjA");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
    global_dpd_->contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    global_dpd_->contract244(&tIA, &D, &Y, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbeJ");
    global_dpd_->contract424(&Y, &tia, &W, 3, 0, 0, 1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y);
    /*** BAAB ***/
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WmBEj");
    global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tjAIb");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
    global_dpd_->contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    global_dpd_->contract244(&tia, &D, &Y, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WmBEj");
    global_dpd_->contract424(&Y, &tIA, &W, 3, 0, 0, 1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y);
    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);
  }
  else if(params.ref == 2) { /** UHF **/
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");
    /*** AAAA ***/
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 20, 20, 20, 20, 0, "WMBEJ");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 20, 20, 20, 20, 0, "D <IJ||AB> (IA,JB)");
    global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 20, 20, 20, 20, 0, "tIAJB");
    global_dpd_->contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 20, 30, 20, 30, 0, "D <Ij|Ab> (IA,jb)");
    global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
    global_dpd_->contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 20, 0, 20, 0, 0, "Y (ME,JN)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 20, 21, 20, 21, 0, "D <IJ||AB> (IA,BJ)");
    global_dpd_->contract244(&tIA, &D, &Y, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 20, 20, 20, 20, 0, "WMBEJ");
    global_dpd_->contract424(&Y, &tIA, &W, 3, 0, 0, -1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y);
    /*** BBBB ***/
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 30, 30, 30, 30, 0, "Wmbej");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 30, 30, 30, 30, 0, "D <ij||ab> (ia,jb)");
    global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 30, 30, 30, 30, 0, "tiajb");
    global_dpd_->contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 30, 20, 30, 20, 0, "D <Ij|Ab> (ia,JB)");
    global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 30, 20, 30, 20, 0, "tiaJB");
    global_dpd_->contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 30, 10, 30, 10, 0, "Y (me,jn)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 30, 31, 30, 31, 0, "D <ij||ab> (ia,bj)");
    global_dpd_->contract244(&tia, &D, &Y, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 30, 30, 30, 30, 0, "Wmbej");
    global_dpd_->contract424(&Y, &tia, &W, 3, 0, 0, -1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y);
    /*** ABAB ***/
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 20, 30, 20, 30, 0, "WMbEj");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 20, 30, 20, 30, 0, "D <Ij|Ab> (IA,jb)");
    global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 30, 30, 30, 30, 0, "tiajb");
    global_dpd_->contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 20, 20, 20, 20, 0, "D <IJ||AB> (IA,JB)");
    global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 30, 20, 30, 20, 0, "tiaJB");
    global_dpd_->contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 20, 10, 20, 10, 0, "Y (ME,jn)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 20, 31, 20, 31, 0, "D <Ij|Ab> (IA,bj)");
    global_dpd_->contract244(&tia, &D, &Y, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 20, 30, 20, 30, 0, "WMbEj");
    global_dpd_->contract424(&Y, &tia, &W, 3, 0, 0, -1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y);
    /*** BABA ***/
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 30, 20, 30, 20, 0, "WmBeJ");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 30, 20, 30, 20, 0, "D <Ij|Ab> (ia,JB)");
    global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 20, 20, 20, 20, 0, "tIAJB");
    global_dpd_->contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 30, 30, 30, 30, 0, "D <ij||ab> (ia,jb)");
    global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
    global_dpd_->contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 30, 0, 30, 0, 0, "Y (me,JN)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 30, 21, 30, 21, 0, "D <Ij|Ab> (ia,BJ)");
    global_dpd_->contract244(&tIA, &D, &Y, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 30, 20, 30, 20, 0, "WmBeJ");
    global_dpd_->contract424(&Y, &tIA, &W, 3, 0, 0, -1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y);
    /*** ABBA ***/
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 24, 24, 24, 24, 0, "WMbeJ");
    global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 24, 27, 24, 27, 0, "tIbjA");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 24, 27, 24, 27, 0, "D <Ij|Ab> (Ib,jA)");
    global_dpd_->contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 24, 22, 24, 22, 0, "Y (Me,Jn)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 24, 26, 24, 26, 0, "D <Ij|Ab> (Ib,Aj)");
    global_dpd_->contract244(&tIA, &D, &Y, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 24, 24, 24, 24, 0, "WMbeJ");
    global_dpd_->contract424(&Y, &tia, &W, 3, 0, 0, 1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y);
    /*** BAAB ***/
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 27, 27, 27, 27, 0, "WmBEj");
    global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 27, 24, 27, 24, 0, "tiBJa");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 27, 24, 27, 24, 0, "D <iJ|aB> (iB,Ja)");
    global_dpd_->contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 27, 23, 27, 23, 0, "Y (mE,jN)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 27, 25, 27, 25, 0, "D <iJ|aB> (iB,aJ)");
    global_dpd_->contract244(&tia, &D, &Y, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 27, 27, 27, 27, 0, "WmBEj");
    global_dpd_->contract424(&Y, &tIA, &W, 3, 0, 0, 1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y);
    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);
  }
  timer_off("X->Wmbej");
}

void purge_Wabij(void) {
  dpdfile4 W;
  int *occpi, *virtpi;
  int h, a, b, e, f, i, j, m, n;
  int    A, B, E, F, I, J, M, N;
  int mn, ei, ma, ef, me, jb, mb, ij, ab;
  int asym, bsym, esym, fsym, isym, jsym, msym, nsym;
  int *occ_off, *vir_off;
  int *occ_sym, *vir_sym;
  int *openpi, nirreps;

  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi; virtpi = moinfo.virtpi;
  occ_off = moinfo.occ_off; vir_off = moinfo.vir_off;
  occ_sym = moinfo.occ_sym; vir_sym = moinfo.vir_sym;
  openpi = moinfo.openpi;

  /* Purge Wabij (ij,ab) matrix elements */
  global_dpd_->file4_init(&W, PSIF_CC_HBAR, 0, 2, 7, "WABIJ residual");
  for(h=0; h < nirreps; h++) {
    global_dpd_->file4_mat_irrep_init(&W, h);
    global_dpd_->file4_mat_irrep_rd(&W, h);
    for(mn=0; mn < W.params->rowtot[h]; mn++) {
      m = W.params->roworb[h][mn][0];
      n = W.params->roworb[h][mn][1];
      msym = W.params->psym[m];
      nsym = W.params->qsym[n];
      M = m - occ_off[msym];
      N = n - occ_off[nsym];
      for(ab=0; ab < W.params->coltot[h]; ab++) {
        a = W.params->colorb[h][ab][0];
        b = W.params->colorb[h][ab][1];
        asym = W.params->rsym[a];
        bsym = W.params->ssym[b];
        A = a - vir_off[asym];
        B = b - vir_off[bsym];
        if ( (M >= (occpi[msym] - openpi[msym])) ||
             (N >= (occpi[nsym] - openpi[nsym]))  ||
             (A >= (virtpi[asym] - openpi[asym])) ||
             (B >= (virtpi[bsym] - openpi[bsym])) )
                W.matrix[h][mn][ab] = 0.0;
      }
    }
    global_dpd_->file4_mat_irrep_wrt(&W, h);
    global_dpd_->file4_mat_irrep_close(&W, h);
  }
  global_dpd_->file4_close(&W);

  /* Purge Wabij (ij,ab) matrix elements */
  global_dpd_->file4_init(&W, PSIF_CC_HBAR, 0, 2, 7, "Wabij residual");
  for(h=0; h < nirreps; h++) {
    global_dpd_->file4_mat_irrep_init(&W, h);
    global_dpd_->file4_mat_irrep_rd(&W, h);
    for(mn=0; mn < W.params->rowtot[h]; mn++) {
      m = W.params->roworb[h][mn][0];
      n = W.params->roworb[h][mn][1];
      msym = W.params->psym[m];
      nsym = W.params->qsym[n];
      M = m - occ_off[msym];
      N = n - occ_off[nsym];
      for(ab=0; ab < W.params->coltot[h]; ab++) {
        a = W.params->colorb[h][ab][0];
        b = W.params->colorb[h][ab][1];
        asym = W.params->rsym[a];
        bsym = W.params->ssym[b];
        A = a - vir_off[asym];
        B = b - vir_off[bsym];
        if ( (M >= (occpi[msym] - openpi[msym])) ||
             (N >= (occpi[nsym] - openpi[nsym])) ||
             (A >= (virtpi[asym] - openpi[asym])) ||
             (B >= (virtpi[bsym] - openpi[bsym])) )
                W.matrix[h][mn][ab] = 0.0;
      }
    }
    global_dpd_->file4_mat_irrep_wrt(&W, h);
    global_dpd_->file4_mat_irrep_close(&W, h);
  }
  global_dpd_->file4_close(&W);

  /* Purge Wabij (ij,ab) matrix elements */
  global_dpd_->file4_init(&W, PSIF_CC_HBAR, 0, 0, 5, "WAbIj residual");
  for(h=0; h < nirreps; h++) {
    global_dpd_->file4_mat_irrep_init(&W, h);
    global_dpd_->file4_mat_irrep_rd(&W, h);
    for(mn=0; mn < W.params->rowtot[h]; mn++) {
      m = W.params->roworb[h][mn][0];
      n = W.params->roworb[h][mn][1];
      msym = W.params->psym[m];
      nsym = W.params->qsym[n];
      M = m - occ_off[msym];
      N = n - occ_off[nsym];
      for(ab=0; ab < W.params->coltot[h]; ab++) {
        a = W.params->colorb[h][ab][0];
        b = W.params->colorb[h][ab][1];
        asym = W.params->rsym[a];
        bsym = W.params->ssym[b];
        A = a - vir_off[asym];
        B = b - vir_off[bsym];
        if ( (M >= (occpi[msym] - openpi[msym])) ||
             (N >= (occpi[nsym] - openpi[nsym])) ||
             (A >= (virtpi[asym] - openpi[asym])) ||
             (B >= (virtpi[bsym] - openpi[bsym])) )
                W.matrix[h][mn][ab] = 0.0;
      }
    }
    global_dpd_->file4_mat_irrep_wrt(&W, h);
    global_dpd_->file4_mat_irrep_close(&W, h);
  }
  global_dpd_->file4_close(&W);
  return;
}

}} // namespace psi::cchbar
