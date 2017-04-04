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
#include <strings.h>
#include <string.h>
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

/* onepdm(): Computes the one-particle density matrix for CC wave functions.
**
** The spin-orbital expressions for the onepdm components are:
**
** D_ij = -1/2 t_im^ef L^jm_ef - t_i^e L^j_e
**
** D_ab = 1/2 L^mn_ae t_mn^be + L^m_a t_m^b
**
** D_ia = t_i^a + (t_im^ae - t_i^e t_m^a) L^m_e
**        - 1/2 L^mn_ef (t_in^ef t_m^a + t_i^e t_mn^af)
**
** D_ai = L^i_a
**
** [cf. Gauss and Stanton, JCP 103, 3561-3577 (1995).]
**
** TDC, July 2002
*/

void onepdm(struct RHO_Params rho_params)
{
  dpdfile2 D, D1, T1, L1, Z;
  dpdbuf4 T2, L2;
  double trace=0.0, dot_AI, dot_IA, dot_ai, dot_ia;
  double factor=0.0;
  bool L2_T2_T_F=true;

  // L2 * T2 * T * F is absent for CC2 Lagrangian
  if(params.wfn == "CC2" && params.dertype == 1) L2_T2_T_F = false;  

  if(params.ref == 0 || params.ref == 1) { /** RHF/ROHF **/
    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, rho_params.DIJ_lbl);
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
    global_dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 0, 7, 2, 7, 0, "LIJAB");
    global_dpd_->contract442(&T2, &L2, &D, 0, 0, -1.0, 0.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 0, 5, 0, 5, 0, "LIjAb");
    global_dpd_->contract442(&T2, &L2, &D, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "LIA");
    global_dpd_->contract222(&T1, &L1, &D, 0, 0, -1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->file2_close(&T1);
    trace += global_dpd_->file2_trace(&D);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, rho_params.Dij_lbl);
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tijab");
    global_dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 0, 7, 2, 7, 0, "Lijab");
    global_dpd_->contract442(&T2, &L2, &D, 0, 0, -1.0, 0.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
    global_dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 0, 5, 0, 5, 0, "LiJaB");
    global_dpd_->contract442(&T2, &L2, &D, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
    global_dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "Lia");
    global_dpd_->contract222(&T1, &L1, &D, 0, 0, -1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->file2_close(&T1);
    trace += global_dpd_->file2_trace(&D);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, rho_params.DAB_lbl);
    global_dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 2, 5, 2, 7, 0, "LIJAB");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    global_dpd_->contract442(&L2, &T2, &D, 3, 3, 1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 0, 5, 0, 5, 0, "LiJaB");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
    global_dpd_->contract442(&L2, &T2, &D, 3, 3, 1.0, 1.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "LIA");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract222(&L1, &T1, &D, 1, 1, 1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->file2_close(&T1);
    trace += global_dpd_->file2_trace(&D);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, rho_params.Dab_lbl);
    global_dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 2, 5, 2, 7, 0, "Lijab");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tijab");
    global_dpd_->contract442(&L2, &T2, &D, 3, 3, 1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 0, 5, 0, 5, 0, "LIjAb");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->contract442(&L2, &T2, &D, 3, 3, 1.0, 1.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "Lia");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
    global_dpd_->contract222(&L1, &T1, &D, 1, 1, 1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->file2_close(&T1);
    trace += global_dpd_->file2_trace(&D);
    global_dpd_->file2_close(&D);
    /* add the T3 contributions to CCSD(T) onepdm calculated
     in cctriples*/
    if(params.wfn == "CCSD_T" && params.ref == 0) {
      global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, rho_params.DIJ_lbl);
      global_dpd_->file2_init(&D1, PSIF_CC_OEI, 0, 0, 0, "DIJ(T)");
      global_dpd_->file2_axpy(&D1, &D, 1.0, 0);
      global_dpd_->file2_close(&D);
      global_dpd_->file2_close(&D1);

      global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, rho_params.Dij_lbl);
      global_dpd_->file2_init(&D1, PSIF_CC_OEI, 0, 0, 0, "DIJ(T)");
      global_dpd_->file2_axpy(&D1, &D, 1.0, 0);
      global_dpd_->file2_close(&D);
      global_dpd_->file2_close(&D1);

      global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, rho_params.DAB_lbl);
      global_dpd_->file2_init(&D1, PSIF_CC_OEI, 0, 1, 1, "DAB(T)");
      global_dpd_->file2_axpy(&D1, &D, 1.0, 0);
      global_dpd_->file2_close(&D);
      global_dpd_->file2_close(&D1);

      global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, rho_params.Dab_lbl);
      global_dpd_->file2_init(&D1, PSIF_CC_OEI, 0, 1, 1, "DAB(T)");
      global_dpd_->file2_axpy(&D1, &D, 1.0, 0);
      global_dpd_->file2_close(&D);
      global_dpd_->file2_close(&D1);
    }

    /*outfile->Printf( "\n\tTrace of onepdm = %20.15f\n", trace);*/

    /* This term is * L0 = 0 for excited states */
    if (rho_params.L_ground) {
      global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
      global_dpd_->file2_copy(&T1, PSIF_CC_OEI, rho_params.DIA_lbl);
      global_dpd_->file2_close(&T1);
    }
    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DIA_lbl);

    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
    global_dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "LIA");
    global_dpd_->dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "Lia");
    global_dpd_->dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "LIA");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&Z, PSIF_CC_TMP0, 0, 0, 0, "Z(I,M)");
    global_dpd_->contract222(&T1, &L1, &Z, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_close(&Z);

    if(L2_T2_T_F){
    /*  D(I,A) << L2(MN,EF) T2(IN,EF) T(M,A) + L2(Mn,Ef) T2(In,Ef) T(M,A) */
    global_dpd_->file2_init(&Z, PSIF_CC_TMP0, 0, 0, 0, "Z(I,M)");
    global_dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 0, 7, 2, 7, 0, "LIJAB");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
    global_dpd_->contract442(&T2, &L2, &Z, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 0, 5, 0, 5, 0, "LIjAb");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->contract442(&T2, &L2, &Z, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    global_dpd_->file2_close(&Z);
    global_dpd_->file2_close(&T1);
    /* T2(MN,AF) L2(MN,EF) T(I,E) + T2(Mn,Af) L2(Mn,Ef) T(I,E) */
    global_dpd_->file2_init(&Z, PSIF_CC_TMP0, 0, 1, 1, "Z(A,E)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    global_dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 2, 5, 2, 7, 0, "LIJAB");
    global_dpd_->contract442(&T2, &L2, &Z, 2, 2, 1.0, 0.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 0, 5, 0, 5, 0, "LIjAb");
    global_dpd_->contract442(&T2, &L2, &Z, 2, 2, 1.0, 1.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract222(&T1, &Z, &D, 0, 0, -1.0, 1.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_close(&Z);
    global_dpd_->file2_close(&D);
    }

    /* This term is * L0 = 0 for excited states */
    if (rho_params.L_ground) {
      global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
      global_dpd_->file2_copy(&T1, PSIF_CC_OEI, rho_params.Dia_lbl);
      global_dpd_->file2_close(&T1);
    }
    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.Dia_lbl);

    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 2, 7, 0, "tijab");
    global_dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "Lia");
    global_dpd_->dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
    global_dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "LIA");
    global_dpd_->dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "Lia");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
    global_dpd_->file2_init(&Z, PSIF_CC_TMP0, 0, 0, 0, "Z(i,m)");
    global_dpd_->contract222(&T1, &L1, &Z, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_close(&Z);

    if (L2_T2_T_F){
    global_dpd_->file2_init(&Z, PSIF_CC_TMP0, 0, 0, 0, "Z(i,m)");
    global_dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 0, 7, 2, 7, 0, "Lijab");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tijab");
    global_dpd_->contract442(&T2, &L2, &Z, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 0, 5, 0, 5, 0, "LiJaB");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
    global_dpd_->contract442(&T2, &L2, &Z, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
    global_dpd_->contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    global_dpd_->file2_close(&Z);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_init(&Z, PSIF_CC_TMP0, 0, 1, 1, "Z(a,e)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tijab");
    global_dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 2, 5, 2, 7, 0, "Lijab");
    global_dpd_->contract442(&T2, &L2, &Z, 2, 2, 1.0, 0.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
    global_dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 0, 5, 0, 5, 0, "LiJaB");
    global_dpd_->contract442(&T2, &L2, &Z, 2, 2, 1.0, 1.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
    global_dpd_->contract222(&T1, &Z, &D, 0, 0, -1.0, 1.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_close(&Z);
    global_dpd_->file2_close(&D);
    }

    /* Note that these blocks are still stored occ/vir */
    global_dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "LIA");
    global_dpd_->file2_copy(&L1, PSIF_CC_OEI, rho_params.DAI_lbl);
    global_dpd_->file2_close(&L1);

    global_dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "Lia");
    global_dpd_->file2_copy(&L1, PSIF_CC_OEI, rho_params.Dai_lbl);
    global_dpd_->file2_close(&L1);

    /* Check overlaps */
    /*
    dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
    dot_IA = dpd_file2_dot_self(&D);
    dpd_file2_close(&D);
    dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.Dia_lbl);
    dot_ia = dpd_file2_dot_self(&D);
    dpd_file2_close(&D);
    dpd_file2_init(&D, CC_OEI, 0, 1, 0, rho_params.DAI_lbl);
    dot_AI = dpd_file2_dot_self(&D);
    dpd_file2_close(&D);
    dpd_file2_init(&D, CC_OEI, 0, 1, 0, rho_params.Dai_lbl);
    dot_ai = dpd_file2_dot_self(&D);
    dpd_file2_close(&D);
    outfile->Printf("\tOverlaps of onepdm after ground-state parts added.\n");
    outfile->Printf("\t<DIA|DIA> = %15.10lf     <Dia|Dia> = %15.10lf\n", dot_IA, dot_ia);
    outfile->Printf("\t<DAI|DAI> = %15.10lf     <Dai|Dai> = %15.10lf\n", dot_AI, dot_ai);
    outfile->Printf("\t<Dpq|Dqp> = %15.10lf\n", dot_IA+dot_ia+dot_AI+dot_ai);
    */
  }
  else if(params.ref == 2) { /** UHF **/

    if(params.wfn == "CCSD_T" && params.dertype == 1) {
      /* For CCSD(T) gradients, some density contributions are
	 calculated in cctriples */
      factor = 1.0;
      global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, "DIJ");
      global_dpd_->file2_copy(&D, PSIF_CC_OEI, rho_params.DIJ_lbl);
      global_dpd_->file2_close(&D);
      global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 2, 2, "Dij");
      global_dpd_->file2_copy(&D, PSIF_CC_OEI, rho_params.Dij_lbl);
      global_dpd_->file2_close(&D);
      global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, "DAB");
      global_dpd_->file2_copy(&D, PSIF_CC_OEI, rho_params.DAB_lbl);
      global_dpd_->file2_close(&D);
      global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 3, 3, "Dab");
      global_dpd_->file2_copy(&D, PSIF_CC_OEI, rho_params.Dab_lbl);
      global_dpd_->file2_close(&D);
    }
    else factor = 0.0;

    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, rho_params.DIJ_lbl);
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
    global_dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 0, 7, 2, 7, 0, "LIJAB");
    global_dpd_->contract442(&T2, &L2, &D, 0, 0, -1.0, factor);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    global_dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 22, 28, 22, 28, 0, "LIjAb");
    global_dpd_->contract442(&T2, &L2, &D, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "LIA");
    global_dpd_->contract222(&T1, &L1, &D, 0, 0, -1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->file2_close(&T1);
    trace += global_dpd_->file2_trace(&D);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 2, 2, rho_params.Dij_lbl);
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 17, 12, 17, 0, "tijab");
    global_dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 10, 17, 12, 17, 0, "Lijab");
    global_dpd_->contract442(&T2, &L2, &D, 0, 0, -1.0, factor);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    global_dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 23, 29, 23, 29, 0, "LiJaB");
    global_dpd_->contract442(&T2, &L2, &D, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
    global_dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 2, 3, "Lia");
    global_dpd_->contract222(&T1, &L1, &D, 0, 0, -1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->file2_close(&T1);
    trace += global_dpd_->file2_trace(&D);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, rho_params.DAB_lbl);
    global_dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 2, 5, 2, 7, 0, "LIJAB");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    global_dpd_->contract442(&L2, &T2, &D, 3, 3, 1.0, factor);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 23, 29, 23, 29, 0, "LiJaB");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    global_dpd_->contract442(&L2, &T2, &D, 3, 3, 1.0, 1.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "LIA");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract222(&L1, &T1, &D, 1, 1, 1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->file2_close(&T1);
    trace += global_dpd_->file2_trace(&D);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 3, 3, rho_params.Dab_lbl);
    global_dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 12, 15, 12, 17, 0, "Lijab");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 12, 15, 12, 17, 0, "tijab");
    global_dpd_->contract442(&L2, &T2, &D, 3, 3, 1.0, factor);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 22, 28, 22, 28, 0, "LIjAb");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    global_dpd_->contract442(&L2, &T2, &D, 3, 3, 1.0, 1.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 2, 3, "Lia");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
    global_dpd_->contract222(&L1, &T1, &D, 1, 1, 1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->file2_close(&T1);
    trace += global_dpd_->file2_trace(&D);
    global_dpd_->file2_close(&D);

    /*outfile->Printf( "\n\tTrace of onepdm = %20.15f\n", trace);*/

    /* This term is * L0 = 0 for excited states */
    if (rho_params.L_ground) {
      global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
      global_dpd_->file2_copy(&T1, PSIF_CC_OEI, rho_params.DIA_lbl);
      global_dpd_->file2_close(&T1);
    }
    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DIA_lbl);

    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
    global_dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "LIA");
    global_dpd_->dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    global_dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 2, 3, "Lia");
    global_dpd_->dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "LIA");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&Z, PSIF_CC_TMP0, 0, 0, 0, "Z(I,M)");
    global_dpd_->contract222(&T1, &L1, &Z, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_close(&Z);

    if (L2_T2_T_F){
    /*  D(I,A) << L2(MN,EF) T2(IN,EF) T(M,A) + L2(Mn,Ef) T2(In,Ef) T(M,A) */
    global_dpd_->file2_init(&Z, PSIF_CC_TMP0, 0, 0, 0, "Z(I,M)");
    global_dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 0, 7, 2, 7, 0, "LIJAB");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
    global_dpd_->contract442(&T2, &L2, &Z, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 22, 28, 22, 28, 0, "LIjAb");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    global_dpd_->contract442(&T2, &L2, &Z, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    global_dpd_->file2_close(&Z);
    global_dpd_->file2_close(&T1);
    /* T2(MN,AF) L2(MN,EF) T(I,E) + T2(Mn,Af) L2(Mn,Ef) T(I,E) */
    global_dpd_->file2_init(&Z, PSIF_CC_TMP0, 0, 1, 1, "Z(A,E)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    global_dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 2, 5, 2, 7, 0, "LIJAB");
    global_dpd_->contract442(&T2, &L2, &Z, 2, 2, 1.0, 0.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    global_dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 22, 28, 22, 28, 0, "LIjAb");
    global_dpd_->contract442(&T2, &L2, &Z, 2, 2, 1.0, 1.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract222(&T1, &Z, &D, 0, 0, -1.0, 1.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_close(&Z);
    global_dpd_->file2_close(&D);
    }

    /* This term is * L0 = 0 for excited states */
    if (rho_params.L_ground) {
      global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
      global_dpd_->file2_copy(&T1, PSIF_CC_OEI, rho_params.Dia_lbl);
      global_dpd_->file2_close(&T1);
    }
    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 2, 3, rho_params.Dia_lbl);

    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 15, 12, 17, 0, "tijab");
    global_dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 2, 3, "Lia");
    global_dpd_->dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    global_dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "LIA");
    global_dpd_->dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 2, 3, "Lia");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
    global_dpd_->file2_init(&Z, PSIF_CC_TMP0, 0, 2, 2, "Z(i,m)");
    global_dpd_->contract222(&T1, &L1, &Z, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_close(&Z);

    if (L2_T2_T_F){
    global_dpd_->file2_init(&Z, PSIF_CC_TMP0, 0, 2, 2, "Z(i,m)");
    global_dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 10, 17, 12, 17, 0, "Lijab");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 17, 12, 17, 0, "tijab");
    global_dpd_->contract442(&T2, &L2, &Z, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 23, 29, 23, 29, 0, "LiJaB");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    global_dpd_->contract442(&T2, &L2, &Z, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
    global_dpd_->contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    global_dpd_->file2_close(&Z);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_init(&Z, PSIF_CC_TMP0, 0, 3, 3, "Z(a,e)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 12, 15, 12, 17, 0, "tijab");
    global_dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 12, 15, 12, 17, 0, "Lijab");
    global_dpd_->contract442(&T2, &L2, &Z, 2, 2, 1.0, 0.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    global_dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 23, 29, 23, 29, 0, "LiJaB");
    global_dpd_->contract442(&T2, &L2, &Z, 2, 2, 1.0, 1.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
    global_dpd_->contract222(&T1, &Z, &D, 0, 0, -1.0, 1.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_close(&Z);
    global_dpd_->file2_close(&D);
    }

    /* Note that these blocks are still stored occ/vir */
    global_dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "LIA");
    global_dpd_->file2_copy(&L1, PSIF_CC_OEI, rho_params.DAI_lbl);
    global_dpd_->file2_close(&L1);

    global_dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 2, 3, "Lia");
    global_dpd_->file2_copy(&L1, PSIF_CC_OEI, rho_params.Dai_lbl);
    global_dpd_->file2_close(&L1);

    /* Check overlaps */
    /*
    dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
    dot_IA = dpd_file2_dot_self(&D);
    dpd_file2_close(&D);
    dpd_file2_init(&D, CC_OEI, 0, 2, 3, rho_params.Dia_lbl);
    dot_ia = dpd_file2_dot_self(&D);
    dpd_file2_close(&D);
    dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
    dot_AI = dpd_file2_dot_self(&D);
    dpd_file2_close(&D);
    dpd_file2_init(&D, CC_OEI, 0, 2, 3, rho_params.Dai_lbl);
    dot_ai = dpd_file2_dot_self(&D);
    dpd_file2_close(&D);
    outfile->Printf("\tOverlaps of onepdm after ground-state parts added.\n");
    outfile->Printf("\t<DIA|DIA> = %15.10lf     <Dia|Dia> = %15.10lf\n", dot_IA, dot_ia);
    outfile->Printf("\t<DAI|DAI> = %15.10lf     <Dai|Dai> = %15.10lf\n", dot_AI, dot_ai);
    outfile->Printf("\t<Dpq|Dqp> = %15.10lf\n", dot_IA+dot_ia+dot_AI+dot_ai);
    */
  }
}

}} // namespace psi::ccdensity
