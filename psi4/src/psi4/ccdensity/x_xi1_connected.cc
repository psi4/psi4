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
    \ingroup CCDENSITY
    \brief Enter brief description of file here
*/
#include <cstdio>
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

/* In the case (xi_connected == 0) only diagrams that involve intermediate
   states that are at least triply excited are included,e.g.,
   xi1 += <0|L Hbar|T,Q><T,Q|R|S,D>

   However, if (xi_connected), then Hbar must be connected to R and to the
   density, but triply excited intermediate states are not required.
   Additional terms to xi_1 look like
     xi_1 += <0|L Hbar R|g> with Hbar connected to both R and g
   The new terms all have doubly excited intermediate states,
     xi1 += <0|L Hbar|D><D|R|S>  (Hbar connected to R and S)

   We remove all terms from the xi_1 and xi_2 equations that do not have
   Hbar connected to R and to S,D.  The exceptions are that we keep
   <ij||ab> in xi_2 and Fia in xi_1.  See below.

   The new terms may be beautifully evaluted by:
     xi_1 +=  (E <0|L|D><D|R|S>) => E*<Limae|Rme> minus all the diagrams
        of <0|L Hbar R|S> where Hbar is not connected to S and R.

   We do _not_ substract <Lme|Rme> Fia, because this term adds to the
   normal xi_1 term <Lmnef|Rmnef> to make (1)Fia.  This constant term
   (along with <ij||ab> in xi_2) causes cclambda to solve the ground-state
   lambda equations implicitly at the same time as zeta.

   We set R0=0, in the sense that the ground state density code now
   acts on only zeta, spat out by cclambda, not some linear
   combination (like R0 * L + Zeta).

   We remove all terms in the excited state density code that do not have
   R connected to the density (many of these terms involve L2R1_OV).

   This trick was provided compliments of Dr. John Stanton.

   RAK 4/04
 */

/* double aug_xi_check(dpdfile2 *HIA, dpdfile2 *Hia); */

void x_xi1_connected(void)
{
  dpdfile2 L1, XIA, Xia, HIA, Hia, I1, R1, F1, IME, Ime;
  int L_irr, R_irr, G_irr;
  dpdbuf4 D, R2, L2, H2, I2;
  double dot, tval;

  L_irr = params.L_irr;
  R_irr = params.R_irr;
  G_irr = params.G_irr;

  if (params.ref == 0) { /* RHF */
    global_dpd_->file2_init(&HIA, PSIF_EOM_TMP0, G_irr, 0, 1, "HIA");

    global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 0, 5, 0, 5, 0, "2LIjAb - LIjbA");
    global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
    global_dpd_->dot24(&R1, &L2, &HIA, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&R1);
    global_dpd_->buf4_close(&L2);

    tval =  params.cceom_energy;
    /* outfile->Printf("\nenergy: %15.10lf\n",tval); */
    global_dpd_->file2_scm(&HIA, tval);

    /* -= (Fme Rme) Lia */
    if (R_irr == 0) {
      global_dpd_->file2_init(&F1, PSIF_CC_OEI, 0, 0, 1, "FME");
      global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
      dot = 2.0 * global_dpd_->file2_dot(&F1, &R1);
      global_dpd_->file2_close(&R1);
      global_dpd_->file2_close(&F1);

      global_dpd_->file2_init(&L1, PSIF_CC_GL, L_irr, 0, 1, "LIA");
      global_dpd_->file2_axpy(&L1, &HIA, -dot, 0);
      global_dpd_->file2_close(&L1);
    }

    /* -= - (Rme Lmnea) Fin */
    global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    global_dpd_->file2_init(&F1, PSIF_CC_OEI, 0, 0, 0, "FMI");
    global_dpd_->contract222(&F1, &I1, &HIA, 0, 1, 1.0, 1.0);
    global_dpd_->file2_close(&F1);
    global_dpd_->file2_close(&I1);

    /* -= Rme Lmief Ffa */
    global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    global_dpd_->file2_init(&F1, PSIF_CC_OEI, 0, 1, 1, "FAE");
    global_dpd_->contract222(&I1, &F1, &HIA, 0, 1, -1.0, 1.0);
    global_dpd_->file2_close(&F1);
    global_dpd_->file2_close(&I1);

    /* -= Rme Lmnef Wifan */
    global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "2 W(ME,jb) + W(Me,Jb)");
    global_dpd_->buf4_sort(&H2, PSIF_EOM_TMP_XI, prqs, 0, 5, "2 W(ME,jb) + W(Me,Jb) (Mj,Eb)");
    global_dpd_->buf4_close(&H2);

    global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    global_dpd_->buf4_init(&H2, PSIF_EOM_TMP_XI, 0, 0, 5, 0, 5, 0, "2 W(ME,jb) + W(Me,Jb) (Mj,Eb)");
    global_dpd_->dot24(&I1, &H2, &HIA, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&H2);
    global_dpd_->file2_close(&I1);

    /* -= Limae ( Rmf Fef - Rne Fnm + Rnf Wnefm ) */
    global_dpd_->file2_init(&IME, PSIF_EOM_TMP_XI, R_irr, 0, 1, "IME");

    global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
    global_dpd_->file2_init(&F1, PSIF_CC_OEI, 0, 1, 1, "FAE");
    global_dpd_->contract222(&R1, &F1, &IME, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&F1);
    global_dpd_->file2_close(&R1);

    global_dpd_->file2_init(&F1, PSIF_CC_OEI, 0, 0, 0, "FMI");
    global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
    global_dpd_->contract222(&F1, &R1, &IME, 1, 1, -1.0, 1.0);
    global_dpd_->file2_close(&R1);
    global_dpd_->file2_close(&F1);

    global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
    global_dpd_->buf4_init(&H2, PSIF_EOM_TMP_XI, 0, 0, 5, 0, 5, 0, "2 W(ME,jb) + W(Me,Jb) (Mj,Eb)");
    global_dpd_->dot13(&R1, &H2, &IME, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&H2);
    global_dpd_->file2_close(&R1);

    /* HIA -= LIAME IME + LIAme Ime */
    global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 0, 5, 0, 5, 0, "2LIjAb - LIjbA");
    global_dpd_->dot24(&IME, &L2, &HIA, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&L2);

    global_dpd_->file2_close(&IME);

    /* add to Xi1 */
    /* aug_xi_check(&HIA, &Hia); */

    global_dpd_->file2_init(&XIA, PSIF_EOM_XI, G_irr, 0, 1, "XIA");
    global_dpd_->file2_axpy(&HIA, &XIA, 1.0, 0);
    global_dpd_->file2_close(&XIA);
    global_dpd_->file2_close(&HIA);
  }


  else if (params.ref == 1) { /* ROHF */
    global_dpd_->file2_init(&HIA, PSIF_EOM_TMP0, G_irr, 0, 1, "HIA");
    global_dpd_->file2_init(&Hia, PSIF_EOM_TMP0, G_irr, 0, 1, "Hia");

    global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 0, 5, 2, 7, 0, "LIJAB");
    global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
    global_dpd_->dot24(&R1, &L2, &HIA, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&R1);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "Ria");
    global_dpd_->dot24(&R1, &L2, &HIA, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&R1);
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 0, 5, 2, 7, 0, "Lijab");
    global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "Ria");
    global_dpd_->dot24(&R1, &L2, &Hia, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&R1);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 0, 5, 0, 5, 0, "LiJaB");
    global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
    global_dpd_->dot24(&R1, &L2, &Hia, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&R1);
    global_dpd_->buf4_close(&L2);

    tval =  params.cceom_energy;
    outfile->Printf("\nenergy: %15.10lf\n",tval);
    global_dpd_->file2_scm(&HIA, tval);
    global_dpd_->file2_scm(&Hia, tval);

    /* -= (Fme Rme) Lia */
    if (R_irr == 0) {
    global_dpd_->file2_init(&F1, PSIF_CC_OEI, 0, 0, 1, "FME");
    global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
    dot = global_dpd_->file2_dot(&F1, &R1);
    global_dpd_->file2_close(&R1);
    global_dpd_->file2_close(&F1);
    global_dpd_->file2_init(&F1, PSIF_CC_OEI, 0, 0, 1, "Fme");
    global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "Ria");
    dot += global_dpd_->file2_dot(&F1, &R1);
    global_dpd_->file2_close(&R1);
    global_dpd_->file2_close(&F1);

    global_dpd_->file2_init(&L1, PSIF_CC_GL, L_irr, 0, 1, "LIA");
    global_dpd_->file2_axpy(&L1, &HIA, -dot, 0);
    global_dpd_->file2_close(&L1);
    global_dpd_->file2_init(&L1, PSIF_CC_GL, L_irr, 0, 1, "Lia");
    global_dpd_->file2_axpy(&L1, &Hia, -dot, 0);
    global_dpd_->file2_close(&L1);
    }

    /* -= - (Rme Lmnea) Fin */
    global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    global_dpd_->file2_init(&F1, PSIF_CC_OEI, 0, 0, 0, "FMI");
    global_dpd_->contract222(&F1, &I1, &HIA, 0, 1, 1.0, 1.0);
    global_dpd_->file2_close(&F1);
    global_dpd_->file2_close(&I1);

    global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 0, 1, "L2R1_ov");
    global_dpd_->file2_init(&F1, PSIF_CC_OEI, 0, 0, 0, "Fmi");
    global_dpd_->contract222(&F1, &I1, &Hia, 0, 1, 1.0, 1.0);
    global_dpd_->file2_close(&F1);
    global_dpd_->file2_close(&I1);

    /* -= Rme Lmief Ffa */
    global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    global_dpd_->file2_init(&F1, PSIF_CC_OEI, 0, 1, 1, "FAE");
    global_dpd_->contract222(&I1, &F1, &HIA, 0, 1, -1.0, 1.0);
    global_dpd_->file2_close(&F1);
    global_dpd_->file2_close(&I1);

    global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 0, 1, "L2R1_ov");
    global_dpd_->file2_init(&F1, PSIF_CC_OEI, 0, 1, 1, "Fae");
    global_dpd_->contract222(&I1, &F1, &Hia, 0, 1, -1.0, 1.0);
    global_dpd_->file2_close(&F1);
    global_dpd_->file2_close(&I1);

    /* -= Rme Lmnef Wifan */
    global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMBEJ");
    global_dpd_->buf4_sort(&H2, PSIF_EOM_TMP_XI, prqs, 0, 5, "WMBEJ (MJ,EB)");
    global_dpd_->buf4_close(&H2);
    global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "Wmbej");
    global_dpd_->buf4_sort(&H2, PSIF_EOM_TMP_XI, prqs, 0, 5, "Wmbej (mj,eb)");
    global_dpd_->buf4_close(&H2);
    global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbEj");
    global_dpd_->buf4_sort(&H2, PSIF_EOM_TMP_XI, prqs, 0, 5, "WMbEj (Mj,Eb)");
    global_dpd_->buf4_close(&H2);
    global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WmBeJ");
    global_dpd_->buf4_sort(&H2, PSIF_EOM_TMP_XI, prqs, 0, 5, "WmBeJ (mJ,eB)");
    global_dpd_->buf4_close(&H2);

    global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    global_dpd_->buf4_init(&H2, PSIF_EOM_TMP_XI, 0, 0, 5, 0, 5, 0, "WMBEJ (MJ,EB)");
    global_dpd_->dot24(&I1, &H2, &HIA, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&H2);
    global_dpd_->buf4_init(&H2, PSIF_EOM_TMP_XI, 0, 0, 5, 0, 5, 0, "WmBeJ (mJ,eB)");
    global_dpd_->dot24(&I1, &H2, &Hia, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&H2);
    global_dpd_->file2_close(&I1);

    global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 0, 1, "L2R1_ov");
    global_dpd_->buf4_init(&H2, PSIF_EOM_TMP_XI, 0, 0, 5, 0, 5, 0, "Wmbej (mj,eb)");
    global_dpd_->dot24(&I1, &H2, &Hia, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&H2);
    global_dpd_->buf4_init(&H2, PSIF_EOM_TMP_XI, 0, 0, 5, 0, 5, 0, "WMbEj (Mj,Eb)");
    global_dpd_->dot24(&I1, &H2, &HIA, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&H2);
    global_dpd_->file2_close(&I1);

    /* -= Limae ( Rmf Fef - Rne Fnm + Rnf Wnefm ) */
    global_dpd_->file2_init(&IME, PSIF_EOM_TMP_XI, R_irr, 0, 1, "IME");
    global_dpd_->file2_init(&Ime, PSIF_EOM_TMP_XI, R_irr, 0, 1, "Ime");

    global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
    global_dpd_->file2_init(&F1, PSIF_CC_OEI, 0, 1, 1, "FAE");
    global_dpd_->contract222(&R1, &F1, &IME, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&F1);
    global_dpd_->file2_close(&R1);
    global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "Ria");
    global_dpd_->file2_init(&F1, PSIF_CC_OEI, 0, 1, 1, "Fae");
    global_dpd_->contract222(&R1, &F1, &Ime, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&F1);
    global_dpd_->file2_close(&R1);

    global_dpd_->file2_init(&F1, PSIF_CC_OEI, 0, 0, 0, "FMI");
    global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
    global_dpd_->contract222(&F1, &R1, &IME, 1, 1, -1.0, 1.0);
    global_dpd_->file2_close(&R1);
    global_dpd_->file2_close(&F1);
    global_dpd_->file2_init(&F1, PSIF_CC_OEI, 0, 0, 0, "Fmi");
    global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "Ria");
    global_dpd_->contract222(&F1, &R1, &Ime, 1, 1, -1.0, 1.0);
    global_dpd_->file2_close(&R1);
    global_dpd_->file2_close(&F1);

    global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
    global_dpd_->buf4_init(&H2, PSIF_EOM_TMP_XI, 0, 0, 5, 0, 5, 0, "WMBEJ (MJ,EB)");
    global_dpd_->dot13(&R1, &H2, &IME, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&H2);
    global_dpd_->buf4_init(&H2, PSIF_EOM_TMP_XI, 0, 0, 5, 0, 5, 0, "WMbEj (Mj,Eb)");
    global_dpd_->dot13(&R1, &H2, &Ime, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&H2);
    global_dpd_->file2_close(&R1);

    global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "Ria");
    global_dpd_->buf4_init(&H2, PSIF_EOM_TMP_XI, 0, 0, 5, 0, 5, 0, "Wmbej (mj,eb)");
    global_dpd_->dot13(&R1, &H2, &Ime, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&H2);
    global_dpd_->buf4_init(&H2, PSIF_EOM_TMP_XI, 0, 0, 5, 0, 5, 0, "WmBeJ (mJ,eB)");
    global_dpd_->dot13(&R1, &H2, &IME, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&H2);
    global_dpd_->file2_close(&R1);

    /* HIA -= LIAME IME + LIAme Ime */
    global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 0, 5, 2, 7, 0, "LIJAB");
    global_dpd_->dot24(&IME, &L2, &HIA, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    global_dpd_->dot24(&Ime, &L2, &HIA, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 0, 5, 2, 7, 0, "Lijab");
    global_dpd_->dot24(&Ime, &L2, &Hia, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 0, 5, 0, 5, 0, "LiJaB");
    global_dpd_->dot24(&IME, &L2, &Hia, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&L2);

    global_dpd_->file2_close(&IME);
    global_dpd_->file2_close(&Ime);

    /* add to Xi1 */
    /* aug_xi_check(&HIA, &Hia); */

    global_dpd_->file2_init(&XIA, PSIF_EOM_XI, G_irr, 0, 1, "XIA");
    global_dpd_->file2_init(&Xia, PSIF_EOM_XI, G_irr, 0, 1, "Xia");
    global_dpd_->file2_axpy(&HIA, &XIA, 1.0, 0);
    global_dpd_->file2_axpy(&Hia, &Xia, 1.0, 0);
    global_dpd_->file2_close(&XIA);
    global_dpd_->file2_close(&Xia);

    global_dpd_->file2_close(&HIA);
    global_dpd_->file2_close(&Hia);
  }


  else { /* UHF */
    global_dpd_->file2_init(&HIA, PSIF_EOM_TMP0, G_irr, 0, 1, "HIA");
    global_dpd_->file2_init(&Hia, PSIF_EOM_TMP0, G_irr, 2, 3, "Hia");

    global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 0, 5, 2, 7, 0, "LIJAB");
    global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
    global_dpd_->dot24(&R1, &L2, &HIA, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&R1);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 22, 28, 22, 28, 0, "LIjAb");
    global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 2, 3, "Ria");
    global_dpd_->dot24(&R1, &L2, &HIA, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&R1);
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 10, 15, 12, 17, 0, "Lijab");
    global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 2, 3, "Ria");
    global_dpd_->dot24(&R1, &L2, &Hia, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&R1);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 23, 29, 23, 29, 0, "LiJaB");
    global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
    global_dpd_->dot24(&R1, &L2, &Hia, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&R1);
    global_dpd_->buf4_close(&L2);

    tval =  params.cceom_energy;
    outfile->Printf("\nenergy: %15.10lf\n",tval);
    global_dpd_->file2_scm(&HIA, tval);
    global_dpd_->file2_scm(&Hia, tval);

    /* -= (Fme Rme) Lia */
    if (R_irr == 0) {
    global_dpd_->file2_init(&F1, PSIF_CC_OEI, 0, 0, 1, "FME");
    global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
    dot = global_dpd_->file2_dot(&F1, &R1);
    global_dpd_->file2_close(&R1);
    global_dpd_->file2_close(&F1);
    global_dpd_->file2_init(&F1, PSIF_CC_OEI, 0, 2, 3, "Fme");
    global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 2, 3, "Ria");
    dot += global_dpd_->file2_dot(&F1, &R1);
    global_dpd_->file2_close(&R1);
    global_dpd_->file2_close(&F1);

    global_dpd_->file2_init(&L1, PSIF_CC_GL, L_irr, 0, 1, "LIA");
    global_dpd_->file2_axpy(&L1, &HIA, -dot, 0);
    global_dpd_->file2_close(&L1);
    global_dpd_->file2_init(&L1, PSIF_CC_GL, L_irr, 2, 3, "Lia");
    global_dpd_->file2_axpy(&L1, &Hia, -dot, 0);
    global_dpd_->file2_close(&L1);
    }

    /* -= - (Rme Lmnea) Fin */
    global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    global_dpd_->file2_init(&F1, PSIF_CC_OEI, 0, 0, 0, "FMI");
    global_dpd_->contract222(&F1, &I1, &HIA, 0, 1, 1.0, 1.0);
    global_dpd_->file2_close(&F1);
    global_dpd_->file2_close(&I1);

    global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 2, 3, "L2R1_ov");
    global_dpd_->file2_init(&F1, PSIF_CC_OEI, 0, 2, 2, "Fmi");
    global_dpd_->contract222(&F1, &I1, &Hia, 0, 1, 1.0, 1.0);
    global_dpd_->file2_close(&F1);
    global_dpd_->file2_close(&I1);

    /* -= Rme Lmief Ffa */
    global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    global_dpd_->file2_init(&F1, PSIF_CC_OEI, 0, 1, 1, "FAE");
    global_dpd_->contract222(&I1, &F1, &HIA, 0, 1, -1.0, 1.0);
    global_dpd_->file2_close(&F1);
    global_dpd_->file2_close(&I1);

    global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 2, 3, "L2R1_ov");
    global_dpd_->file2_init(&F1, PSIF_CC_OEI, 0, 3, 3, "Fae");
    global_dpd_->contract222(&I1, &F1, &Hia, 0, 1, -1.0, 1.0);
    global_dpd_->file2_close(&F1);
    global_dpd_->file2_close(&I1);

    /* -= Rme Lmnef Wifan */
    global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 20, 20, 20, 20, 0, "WMBEJ");
    global_dpd_->buf4_sort(&H2, PSIF_EOM_TMP_XI, prqs, 0, 5, "WMBEJ (MJ,EB)");
    global_dpd_->buf4_close(&H2);
    global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 30, 30, 30, 30, 0, "Wmbej");
    global_dpd_->buf4_sort(&H2, PSIF_EOM_TMP_XI, prqs, 10, 15, "Wmbej (mj,eb)");
    global_dpd_->buf4_close(&H2);
    global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 20, 30, 20, 30, 0, "WMbEj");
    global_dpd_->buf4_sort(&H2, PSIF_EOM_TMP_XI, prqs, 22, 28, "WMbEj (Mj,Eb)");
    global_dpd_->buf4_close(&H2);
    global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 30, 20, 30, 20, 0, "WmBeJ");
    global_dpd_->buf4_sort(&H2, PSIF_EOM_TMP_XI, prqs, 23, 29, "WmBeJ (mJ,eB)");
    global_dpd_->buf4_close(&H2);

    global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    global_dpd_->buf4_init(&H2, PSIF_EOM_TMP_XI, 0, 0, 5, 0, 5, 0, "WMBEJ (MJ,EB)");
    global_dpd_->dot24(&I1, &H2, &HIA, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&H2);
    global_dpd_->buf4_init(&H2, PSIF_EOM_TMP_XI, 0, 23, 29, 23, 29, 0, "WmBeJ (mJ,eB)");
    global_dpd_->dot24(&I1, &H2, &Hia, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&H2);
    global_dpd_->file2_close(&I1);
    global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 2, 3, "L2R1_ov");
    global_dpd_->buf4_init(&H2, PSIF_EOM_TMP_XI, 0, 10, 15, 10, 15, 0, "Wmbej (mj,eb)");
    global_dpd_->dot24(&I1, &H2, &Hia, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&H2);
    global_dpd_->buf4_init(&H2, PSIF_EOM_TMP_XI, 0, 22, 28, 22, 28, 0, "WMbEj (Mj,Eb)");
    global_dpd_->dot24(&I1, &H2, &HIA, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&H2);
    global_dpd_->file2_close(&I1);

    /* -= Limae ( Rmf Fef - Rne Fnm + Rnf Wnefm ) */
    global_dpd_->file2_init(&IME, PSIF_EOM_TMP_XI, R_irr, 0, 1, "IME");
    global_dpd_->file2_init(&Ime, PSIF_EOM_TMP_XI, R_irr, 2, 3, "Ime");

    global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
    global_dpd_->file2_init(&F1, PSIF_CC_OEI, 0, 1, 1, "FAE");
    global_dpd_->contract222(&R1, &F1, &IME, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&F1);
    global_dpd_->file2_close(&R1);
    global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 2, 3, "Ria");
    global_dpd_->file2_init(&F1, PSIF_CC_OEI, 0, 3, 3, "Fae");
    global_dpd_->contract222(&R1, &F1, &Ime, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&F1);
    global_dpd_->file2_close(&R1);

    global_dpd_->file2_init(&F1, PSIF_CC_OEI, 0, 0, 0, "FMI");
    global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
    global_dpd_->contract222(&F1, &R1, &IME, 1, 1, -1.0, 1.0);
    global_dpd_->file2_close(&R1);
    global_dpd_->file2_close(&F1);
    global_dpd_->file2_init(&F1, PSIF_CC_OEI, 0, 2, 2, "Fmi");
    global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 2, 3, "Ria");
    global_dpd_->contract222(&F1, &R1, &Ime, 1, 1, -1.0, 1.0);
    global_dpd_->file2_close(&R1);
    global_dpd_->file2_close(&F1);

    global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
    global_dpd_->buf4_init(&H2, PSIF_EOM_TMP_XI, 0, 0, 5, 0, 5, 0, "WMBEJ (MJ,EB)");
    global_dpd_->dot13(&R1, &H2, &IME, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&H2);
    global_dpd_->buf4_init(&H2, PSIF_EOM_TMP_XI, 0, 22, 28, 22, 28, 0, "WMbEj (Mj,Eb)");
    global_dpd_->dot13(&R1, &H2, &Ime, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&H2);
    global_dpd_->file2_close(&R1);
    global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 2, 3, "Ria");
    global_dpd_->buf4_init(&H2, PSIF_EOM_TMP_XI, 0, 10, 15, 10, 15, 0, "Wmbej (mj,eb)");
    global_dpd_->dot13(&R1, &H2, &Ime, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&H2);
    global_dpd_->buf4_init(&H2, PSIF_EOM_TMP_XI, 0, 23, 29, 23, 29, 0, "WmBeJ (mJ,eB)");
    global_dpd_->dot13(&R1, &H2, &IME, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&H2);
    global_dpd_->file2_close(&R1);

    /* HIA -= LIAME IME + LIAme Ime */
    global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 0, 5, 2, 7, 0, "LIJAB");
    global_dpd_->dot24(&IME, &L2, &HIA, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 22, 28, 22, 28, 0, "LIjAb");
    global_dpd_->dot24(&Ime, &L2, &HIA, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 10, 15, 12, 17, 0, "Lijab");
    global_dpd_->dot24(&Ime, &L2, &Hia, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 23, 29, 23, 29, 0, "LiJaB");
    global_dpd_->dot24(&IME, &L2, &Hia, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&L2);

    global_dpd_->file2_close(&IME);
    global_dpd_->file2_close(&Ime);

    /* aug_xi_check(&HIA, &Hia); */

    /* add to Xi1 */
    global_dpd_->file2_init(&XIA, PSIF_EOM_XI, G_irr, 0, 1, "XIA");
    global_dpd_->file2_init(&Xia, PSIF_EOM_XI, G_irr, 2, 3, "Xia");
    global_dpd_->file2_axpy(&HIA, &XIA, 1.0, 0);
    global_dpd_->file2_axpy(&Hia, &Xia, 1.0, 0);
    global_dpd_->file2_close(&XIA);
    global_dpd_->file2_close(&Xia);

    global_dpd_->file2_close(&HIA);
    global_dpd_->file2_close(&Hia);
  }
}

/*
double aug_xi_check(dpdfile2 *HIA, dpdfile2 *Hia)
{
  double tvalA, tvalB;
  tvalA = tvalB = 0.0;

  if (params.ref == 0) {
    tvalA = dpd_file2_dot_self(HIA);
    outfile->Printf( "<HIA|HIA> = %15.10lf\n", tvalA);
  }
  else {
    tvalB = dpd_file2_dot_self(Hia);
    outfile->Printf( "<HIA|HIA> = %15.10lf\n", tvalA);
  }
  outfile->Printf( "<H1|H1> = %15.10lf\n", tvalA + tvalB);
  return;
}
*/

}} // namespace psi::ccdensity
