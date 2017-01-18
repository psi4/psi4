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
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

    /* Iai(): Build the virtual-occupied block of the orbital Lagrangian
    ** using the expression given in lag.c.  Note that we include an
    ** additional term here referred to as the reference contribution,
    ** 2*fai.  This comes from the general spin-orbital SCF gradient
    ** expression, but for unperturbed canonical Hartree-Fock orbitals
    ** (i.e., RHF and UHF only) this contribution is zero.  However, since
    ** the code to include the terms is trivial, we go ahead and do the
    ** work for all reference types.
    **
    ** Note that the code currently produces only the unique I_IJ terms,
    ** but the actual contractions still need to be spin-adapted for
    ** greater efficiency.
    **
    ** TDC, 2/2008
    */

    void Iai(struct RHO_Params rho_params)
    {
      dpdfile2 F, D, I;
      dpdbuf4 G, Eints, Dints, Cints, Fints, Bints;

      if(params.ref == 0) { /** RHF **/
	/* I'AI <-- sum_J fAJ (DIJ + DJI) + sum_B fAB (DIB + DBI) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 0, "I'AI");

	global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 1, "fIA");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, rho_params.DIJ_lbl);
	global_dpd_->contract222(&F, &D, &I, 1, 0, 1.0, 0.0);
	global_dpd_->contract222(&F, &D, &I, 1, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);

	/* Add reference contribution: I'AI <-- 2 fAI */
	global_dpd_->file2_axpy(&F, &I, 2.0, 1);
	global_dpd_->file2_close(&F);

	global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 1, 1, "fAB");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
	global_dpd_->contract222(&F, &D, &I, 0, 0, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
	global_dpd_->contract222(&F, &D, &I, 0, 0, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->file2_close(&F);

	global_dpd_->file2_close(&I);
      }
      else if(params.ref == 1) { /** ROHF **/

	/* I'AI <-- sum_J fAJ (DIJ + DJI) + sum_B fAB (DIB + DBI) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 0, "I'AI");

	global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 1, "fIA");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, rho_params.DIJ_lbl);
	global_dpd_->contract222(&F, &D, &I, 1, 0, 1.0, 0.0);
	global_dpd_->contract222(&F, &D, &I, 1, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);

	/* Add reference contribution: I'AI <-- 2 fAI */
	global_dpd_->file2_axpy(&F, &I, 2.0, 1);
	global_dpd_->file2_close(&F);

	global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 1, 1, "fAB");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
	global_dpd_->contract222(&F, &D, &I, 0, 0, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
	global_dpd_->contract222(&F, &D, &I, 0, 0, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->file2_close(&F);

	global_dpd_->file2_close(&I);

	/* I'ai <-- sum_j faj (Dij + Dji) + sum_b fab (Dib + Dbi) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 0, "I'ai");

	global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 1, "fia");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, rho_params.Dij_lbl);
	global_dpd_->contract222(&F, &D, &I, 1, 0, 1.0, 0.0);
	global_dpd_->contract222(&F, &D, &I, 1, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);

	/* Add reference contribution: I'ai <-- 2 fai */
	global_dpd_->file2_axpy(&F, &I, 2.0, 1);
	global_dpd_->file2_close(&F);

	global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 1, 1, "fab");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.Dia_lbl);
	global_dpd_->contract222(&F, &D, &I, 0, 0, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.Dai_lbl);
	global_dpd_->contract222(&F, &D, &I, 0, 0, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->file2_close(&F);

	global_dpd_->file2_close(&I);
      }
      else if(params.ref == 2) { /** UHF **/

	/* I'AI <-- sum_J fAJ (DIJ + DJI) + sum_B fAB (DIB + DBI) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 0, "I'AI");

	global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 1, "fIA");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, rho_params.DIJ_lbl);
	global_dpd_->contract222(&F, &D, &I, 1, 0, 1.0, 0.0);
	global_dpd_->contract222(&F, &D, &I, 1, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);

	/* Add reference contribution: I'AI <-- 2 fAI */
	global_dpd_->file2_axpy(&F, &I, 2.0, 1);
	global_dpd_->file2_close(&F);

	global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 1, 1, "fAB");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
	global_dpd_->contract222(&F, &D, &I, 0, 0, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
	global_dpd_->contract222(&F, &D, &I, 0, 0, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->file2_close(&F);

	global_dpd_->file2_close(&I);

	/* I'ai <-- sum_j faj (Dij + Dji) + sum_b fab (Dib + Dbi) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 3, 2, "I'ai");

	global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 2, 3, "fia");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 2, 2, rho_params.Dij_lbl);
	global_dpd_->contract222(&F, &D, &I, 1, 0, 1.0, 0.0);
	global_dpd_->contract222(&F, &D, &I, 1, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);

	/* Add reference contribution: I'ai <-- 2 fai */
	global_dpd_->file2_axpy(&F, &I, 2.0, 1);
	global_dpd_->file2_close(&F);

	global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 3, 3, "fab");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 2, 3, rho_params.Dia_lbl);
	global_dpd_->contract222(&F, &D, &I, 0, 0, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 2, 3, rho_params.Dai_lbl);
	global_dpd_->contract222(&F, &D, &I, 0, 0, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->file2_close(&F);

	global_dpd_->file2_close(&I);

      }



      if(params.ref == 0) { /** RHF **/
	/* I'AI <-- sum_JK <AJ||IK> (D_JK + D_KJ) + sum_jk <Aj|Ik> (D_jk + D_kj) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 0, "I'AI");

	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, rho_params.DIJ_lbl);
	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 1, "E <ai|jk>");
	global_dpd_->dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
	global_dpd_->dot24(&D, &Eints, &I, 1, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&Eints);
	global_dpd_->file2_close(&D);

	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, rho_params.Dij_lbl);
	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
	global_dpd_->dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
	global_dpd_->dot24(&D, &Eints, &I, 1, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&Eints);
	global_dpd_->file2_close(&D);

	global_dpd_->file2_close(&I);
      }
      else if(params.ref == 1) { /** ROHF **/

	/* I'AI <-- sum_JK <AJ||IK> (D_JK + D_KJ) + sum_jk <Aj|Ik> (D_jk + D_kj) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 0, "I'AI");

	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, rho_params.DIJ_lbl);
	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 1, "E <ai|jk>");
	global_dpd_->dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
	global_dpd_->dot24(&D, &Eints, &I, 1, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&Eints);
	global_dpd_->file2_close(&D);

	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, rho_params.Dij_lbl);
	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
	global_dpd_->dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
	global_dpd_->dot24(&D, &Eints, &I, 1, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&Eints);
	global_dpd_->file2_close(&D);

	global_dpd_->file2_close(&I);

	/* I'ai <-- sum_jk <aj||ik> (D_jk + D_kj) + sum_jk <aJ|iK> (D_JK + D_KJ) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 0, "I'ai");

	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, rho_params.Dij_lbl);
	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 1, "E <ai|jk>");
	global_dpd_->dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
	global_dpd_->dot24(&D, &Eints, &I, 1, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&Eints);
	global_dpd_->file2_close(&D);

	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, rho_params.DIJ_lbl);
	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
	global_dpd_->dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
	global_dpd_->dot24(&D, &Eints, &I, 1, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&Eints);
	global_dpd_->file2_close(&D);

	global_dpd_->file2_close(&I);
      }
      else if(params.ref == 2) { /** UHF **/

	/* I'AI <-- sum_JK <AJ||IK> (D_JK + D_KJ) + sum_jk <Aj|Ik> (D_jk + D_kj) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 0, "I'AI");

	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, rho_params.DIJ_lbl);
	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 21, 0, 21, 0, 1, "E <AI|JK>");
	global_dpd_->dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
	global_dpd_->dot24(&D, &Eints, &I, 1, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&Eints);
	global_dpd_->file2_close(&D);

	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 2, 2, rho_params.Dij_lbl);
	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 26, 22, 26, 22, 0, "E <Ai|Jk>");
	global_dpd_->dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
	global_dpd_->dot24(&D, &Eints, &I, 1, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&Eints);
	global_dpd_->file2_close(&D);

	global_dpd_->file2_close(&I);

	/* I'ai <-- sum_jk <aj||ik> (D_jk + D_kj) + sum_jk <aJ|iK> (D_JK + D_KJ) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 3, 2, "I'ai");

	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 2, 2, rho_params.Dij_lbl);
	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 31, 10, 31, 10, 1, "E <ai|jk>");
	global_dpd_->dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
	global_dpd_->dot24(&D, &Eints, &I, 1, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&Eints);
	global_dpd_->file2_close(&D);

	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, rho_params.DIJ_lbl);
	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 24, 22, 24, 22, 0, "E <Ia|Jk>");
	global_dpd_->dot13(&D, &Eints, &I, 0, 0, 1.0, 1.0);
	global_dpd_->dot13(&D, &Eints, &I, 1, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&Eints);
	global_dpd_->file2_close(&D);

	global_dpd_->file2_close(&I);

      }



      if(params.ref == 0) { /** RHF **/
	/* I'AI <-- - sum_JB <JA||IB> (D_JB + D_BJ) + sum_jb <Ij|Ab> (D_jb + D_bj) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 0, "I'AI");

	global_dpd_->buf4_init(&Cints, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
	global_dpd_->dot14(&D, &Cints, &I, 0, 0, -1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
	global_dpd_->dot14(&D, &Cints, &I, 0, 0, -1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->buf4_close(&Cints);

	global_dpd_->buf4_init(&Dints, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.Dia_lbl);
	global_dpd_->dot24(&D, &Dints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.Dai_lbl);
	global_dpd_->dot24(&D, &Dints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->buf4_close(&Dints);

	global_dpd_->file2_close(&I);
      }
      else if(params.ref == 1) { /** ROHF **/

	/* I'AI <-- - sum_JB <JA||IB> (D_JB + D_BJ) + sum_jb <Ij|Ab> (D_jb + D_bj) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 0, "I'AI");

	global_dpd_->buf4_init(&Cints, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
	global_dpd_->dot14(&D, &Cints, &I, 0, 0, -1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
	global_dpd_->dot14(&D, &Cints, &I, 0, 0, -1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->buf4_close(&Cints);

	global_dpd_->buf4_init(&Dints, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.Dia_lbl);
	global_dpd_->dot24(&D, &Dints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.Dai_lbl);
	global_dpd_->dot24(&D, &Dints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->buf4_close(&Dints);

	global_dpd_->file2_close(&I);

	/* I'ai <-- - sum_jb <ja||ib> (D_jb + D_bj) + sum_JB <iJ|aB> (D_JB + D_BJ) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 0, "I'ai");

	global_dpd_->buf4_init(&Cints, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.Dia_lbl);
	global_dpd_->dot14(&D, &Cints, &I, 0, 0, -1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.Dai_lbl);
	global_dpd_->dot14(&D, &Cints, &I, 0, 0, -1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->buf4_close(&Cints);

	global_dpd_->buf4_init(&Dints, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
	global_dpd_->dot24(&D, &Dints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
	global_dpd_->dot24(&D, &Dints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->buf4_close(&Dints);

	global_dpd_->file2_close(&I);
      }
      else if(params.ref == 2) { /** UHF **/

	/* I'AI <-- - sum_JB <JA||IB> (D_JB + D_BJ) + sum_jb <Ij|Ab> (D_jb + D_bj) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 0, "I'AI");

	global_dpd_->buf4_init(&Cints, PSIF_CC_CINTS, 0, 20, 20, 20, 20, 0, "C <IA||JB>");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
	global_dpd_->dot14(&D, &Cints, &I, 0, 0, -1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
	global_dpd_->dot14(&D, &Cints, &I, 0, 0, -1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->buf4_close(&Cints);

	global_dpd_->buf4_init(&Dints, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 2, 3, rho_params.Dia_lbl);
	global_dpd_->dot24(&D, &Dints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 2, 3, rho_params.Dai_lbl);
	global_dpd_->dot24(&D, &Dints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->buf4_close(&Dints);

	global_dpd_->file2_close(&I);

	/* I'ai <-- - sum_jb <ja||ib> (D_jb + D_bj) + sum_JB <iJ|aB> (D_JB + D_BJ) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 3, 2, "I'ai");

	global_dpd_->buf4_init(&Cints, PSIF_CC_CINTS, 0, 30, 30, 30, 30, 0, "C <ia||jb>");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 2, 3, rho_params.Dia_lbl);
	global_dpd_->dot14(&D, &Cints, &I, 0, 0, -1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 2, 3, rho_params.Dai_lbl);
	global_dpd_->dot14(&D, &Cints, &I, 0, 0, -1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->buf4_close(&Cints);

	global_dpd_->buf4_init(&Dints, PSIF_CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
	global_dpd_->dot24(&D, &Dints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
	global_dpd_->dot24(&D, &Dints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->buf4_close(&Dints);

	global_dpd_->file2_close(&I);

      }



      if(params.ref == 0) { /** RHF **/
	/* I'AI <-- sum_BJ <IJ||AB> (D_BJ + D_JB) + sum_bj <Ij|Ab> (D_bj + D_jb) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 0, "I'AI");

	global_dpd_->buf4_init(&Dints, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij||ab>");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
	global_dpd_->dot24(&D, &Dints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
	global_dpd_->dot24(&D, &Dints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->buf4_close(&Dints);

	global_dpd_->buf4_init(&Dints, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.Dai_lbl);
	global_dpd_->dot24(&D, &Dints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.Dia_lbl);
	global_dpd_->dot24(&D, &Dints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->buf4_close(&Dints);

	global_dpd_->file2_close(&I);
      }
      else if(params.ref == 1) { /** ROHF **/

	/* I'AI <-- sum_BJ <IJ||AB> (D_BJ + D_JB) + sum_bj <Ij|Ab> (D_bj + D_jb) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 0, "I'AI");

	global_dpd_->buf4_init(&Dints, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij||ab>");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
	global_dpd_->dot24(&D, &Dints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
	global_dpd_->dot24(&D, &Dints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->buf4_close(&Dints);

	global_dpd_->buf4_init(&Dints, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.Dai_lbl);
	global_dpd_->dot24(&D, &Dints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.Dia_lbl);
	global_dpd_->dot24(&D, &Dints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->buf4_close(&Dints);

	global_dpd_->file2_close(&I);

	/* I'ai <-- sum_bj <ij||ab> (D_bj + D_jb) + sum_BJ <iJ|aB> (D_BJ + D_JB) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 0, "I'ai");

	global_dpd_->buf4_init(&Dints, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij||ab>");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.Dai_lbl);
	global_dpd_->dot24(&D, &Dints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.Dia_lbl);
	global_dpd_->dot24(&D, &Dints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->buf4_close(&Dints);

	global_dpd_->buf4_init(&Dints, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
	global_dpd_->dot24(&D, &Dints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
	global_dpd_->dot24(&D, &Dints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->buf4_close(&Dints);

	global_dpd_->file2_close(&I);
      }
      else if(params.ref == 2) { /** UHF **/

	/* I'AI <-- sum_BJ <IJ||AB> (D_BJ + D_JB) + sum_bj <Ij|Ab> (D_bj + D_jb) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 0, "I'AI");

	global_dpd_->buf4_init(&Dints, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <IJ||AB>");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
	global_dpd_->dot24(&D, &Dints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
	global_dpd_->dot24(&D, &Dints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->buf4_close(&Dints);

	global_dpd_->buf4_init(&Dints, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 2, 3, rho_params.Dai_lbl);
	global_dpd_->dot24(&D, &Dints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 2, 3, rho_params.Dia_lbl);
	global_dpd_->dot24(&D, &Dints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->buf4_close(&Dints);

	global_dpd_->file2_close(&I);

	/* I'ai <-- sum_bj <ij||ab> (D_bj + D_jb) + sum_BJ <iJ|aB> (D_BJ + D_JB) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 3, 2, "I'ai");

	global_dpd_->buf4_init(&Dints, PSIF_CC_DINTS, 0, 10, 15, 10, 15, 0, "D <ij||ab>");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 2, 3, rho_params.Dai_lbl);
	global_dpd_->dot24(&D, &Dints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 2, 3, rho_params.Dia_lbl);
	global_dpd_->dot24(&D, &Dints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->buf4_close(&Dints);

	global_dpd_->buf4_init(&Dints, PSIF_CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
	global_dpd_->dot24(&D, &Dints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
	global_dpd_->dot24(&D, &Dints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->buf4_close(&Dints);

	global_dpd_->file2_close(&I);

      }



      if(params.ref == 0) { /** RHF **/
	/* I'AI <-- sum_BC <IC||AB> (D_BC + D_CB) + sum_bc <Ib|Ac>(D_bc + D_cb) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 0, "I'AI");

	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, rho_params.DAB_lbl);
	global_dpd_->buf4_init(&Fints, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 1, "F <ia|bc>");
	global_dpd_->dot24(&D, &Fints, &I, 1, 1, 1.0, 1.0);
	global_dpd_->dot24(&D, &Fints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->buf4_close(&Fints);
	global_dpd_->file2_close(&D);

	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, rho_params.Dab_lbl);
	global_dpd_->buf4_init(&Fints, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
	global_dpd_->dot24(&D, &Fints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->dot24(&D, &Fints, &I, 1, 1, 1.0, 1.0);
	global_dpd_->buf4_close(&Fints);
	global_dpd_->file2_close(&D);

	global_dpd_->file2_close(&I);
      }
      else if(params.ref == 1) { /** ROHF **/

	/* I'AI <-- sum_BC <IC||AB> (D_BC + D_CB) + sum_bc <Ib|Ac>(D_bc + D_cb) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 0, "I'AI");

	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, rho_params.DAB_lbl);
	global_dpd_->buf4_init(&Fints, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 1, "F <ia|bc>");
	global_dpd_->dot24(&D, &Fints, &I, 1, 1, 1.0, 1.0);
	global_dpd_->dot24(&D, &Fints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->buf4_close(&Fints);
	global_dpd_->file2_close(&D);

	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, rho_params.Dab_lbl);
	global_dpd_->buf4_init(&Fints, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
	global_dpd_->dot24(&D, &Fints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->dot24(&D, &Fints, &I, 1, 1, 1.0, 1.0);
	global_dpd_->buf4_close(&Fints);
	global_dpd_->file2_close(&D);

	global_dpd_->file2_close(&I);

	/* I'ai <-- sum_bc <ic||ab> (D_bc + D_cb) + sum_BC <iB|aC>(D_BC + D_CB) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 0, "I'ai");

	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, rho_params.Dab_lbl);
	global_dpd_->buf4_init(&Fints, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 1, "F <ia|bc>");
	global_dpd_->dot24(&D, &Fints, &I, 1, 1, 1.0, 1.0);
	global_dpd_->dot24(&D, &Fints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->buf4_close(&Fints);
	global_dpd_->file2_close(&D);

	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, rho_params.DAB_lbl);
	global_dpd_->buf4_init(&Fints, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
	global_dpd_->dot24(&D, &Fints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->dot24(&D, &Fints, &I, 1, 1, 1.0, 1.0);
	global_dpd_->buf4_close(&Fints);
	global_dpd_->file2_close(&D);

	global_dpd_->file2_close(&I);
      }
      else if(params.ref == 2) { /** UHF **/

	/* I'AI <-- sum_BC <IC||AB> (D_BC + D_CB) + sum_bc <Ib|Ac>(D_bc + D_cb) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 0, "I'AI");

	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, rho_params.DAB_lbl);
	global_dpd_->buf4_init(&Fints, PSIF_CC_FINTS, 0, 20, 5, 20, 5, 1, "F <IA|BC>");
	global_dpd_->dot24(&D, &Fints, &I, 1, 1, 1.0, 1.0);
	global_dpd_->dot24(&D, &Fints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->buf4_close(&Fints);
	global_dpd_->file2_close(&D);

	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 3, 3, rho_params.Dab_lbl);
	global_dpd_->buf4_init(&Fints, PSIF_CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
	global_dpd_->dot24(&D, &Fints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->dot24(&D, &Fints, &I, 1, 1, 1.0, 1.0);
	global_dpd_->buf4_close(&Fints);
	global_dpd_->file2_close(&D);

	global_dpd_->file2_close(&I);

	/* I'ai <-- sum_bc <ic||ab> (D_bc + D_cb) + sum_BC <iB|aC>(D_BC + D_CB) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 3, 2, "I'ai");

	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 3, 3, rho_params.Dab_lbl);
	global_dpd_->buf4_init(&Fints, PSIF_CC_FINTS, 0, 30, 15, 30, 15, 1, "F <ia|bc>");
	global_dpd_->dot24(&D, &Fints, &I, 1, 1, 1.0, 1.0);
	global_dpd_->dot24(&D, &Fints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->buf4_close(&Fints);
	global_dpd_->file2_close(&D);

	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, rho_params.DAB_lbl);
	global_dpd_->buf4_init(&Fints, PSIF_CC_FINTS, 0, 27, 29, 27, 29, 0, "F <iA|bC>");
	global_dpd_->dot24(&D, &Fints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->dot24(&D, &Fints, &I, 1, 1, 1.0, 1.0);
	global_dpd_->buf4_close(&Fints);
	global_dpd_->file2_close(&D);

	global_dpd_->file2_close(&I);

      }



      if(params.ref == 0) { /** RHF **/
	/* I'AI <-- sum_JKL <AJ||KL> G(IJ,KL) + 2 sum_jKl <Aj|Kl> G(Ij,Kl) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 0, "I'AI");

	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 0, 0, 0, 0, "GIjKl");
	global_dpd_->buf4_scmcopy(&G, PSIF_CC_GAMMA, "2 Gijkl - Gijlk", 2);
	global_dpd_->buf4_sort_axpy(&G, PSIF_CC_GAMMA, pqsr, 0, 0, "2 Gijkl - Gijlk", -1);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 0, 0, 0, 0, "2 Gijkl - Gijlk");
	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
	global_dpd_->contract442(&Eints, &G, &I, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&Eints);
	global_dpd_->buf4_close(&G);

	global_dpd_->file2_close(&I);
      }
      else if(params.ref == 1) { /** ROHF **/

	/* I'AI <-- sum_JKL <AJ||KL> G(IJ,KL) + 2 sum_jKl <Aj|Kl> G(Ij,Kl) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 0, "I'AI");

	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 11, 2, 11, 0, 1, "E <ai|jk>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 2, 2, 2, 0, "GIJKL");
	global_dpd_->contract442(&Eints, &G, &I, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Eints);

	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 0, 0, 0, 0, "GIjKl");
	global_dpd_->contract442(&Eints, &G, &I, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Eints);

	global_dpd_->file2_close(&I);

	/* I'ai <-- sum_jkl <aj||kl> G(ij,kl) + 2 sum_JkL <Lk|Ja> G(Lk,Ji) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 0, "I'ai");

	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 11, 2, 11, 0, 1, "E <ai|jk>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 2, 2, 2, 0, "Gijkl");
	global_dpd_->contract442(&Eints, &G, &I, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Eints);

	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 0, 0, 0, 0, "GIjKl");
	global_dpd_->contract442(&Eints, &G, &I, 3, 3, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Eints);

	global_dpd_->file2_close(&I);
      }
      else if(params.ref == 2) { /** UHF **/

	/* I'AI <-- sum_JKL <AJ||KL> G(IJ,KL) + 2 sum_jKl <Aj|Kl> G(Ij,Kl) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 0, "I'AI");

	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 21, 2, 21, 0, 1, "E <AI|JK>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 2, 2, 2, 0, "GIJKL");
	global_dpd_->contract442(&Eints, &G, &I, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Eints);

	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 26, 22, 26, 22, 0, "E <Ai|Jk>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 22, 22, 22, 22, 0, "GIjKl");
	global_dpd_->contract442(&Eints, &G, &I, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Eints);

	global_dpd_->file2_close(&I);

	/* I'ai <-- sum_jkl <aj||kl> G(ij,kl) + 2 sum_JkL <Lk|Ja> G(Lk,Ji) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 3, 2, "I'ai");

	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 31, 12, 31, 10, 1, "E <ai|jk>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 12, 12, 12, 0, "Gijkl");
	global_dpd_->contract442(&Eints, &G, &I, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Eints);

	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 22, 22, 22, 22, 0, "GIjKl");
	global_dpd_->contract442(&Eints, &G, &I, 3, 3, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Eints);

	global_dpd_->file2_close(&I);

      }



      if(params.ref == 0) { /** RHF **/
	/* I'AI <-- 2 sum_JKB <JA||KB> G(JI,KB) + 2 sum_jkB <jA|kB> G(jI,kB) +
	   2 sum_jKb <Kj|Ab> (Aj,Kb) G(Ij,Kb) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 0, "I'AI");

	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 10, 0, 10, 0, "GIjKa");
	global_dpd_->buf4_scmcopy(&G, PSIF_CC_GAMMA, "2 Gijka - Gjika", 2);
	global_dpd_->buf4_sort_axpy(&G, PSIF_CC_GAMMA, qprs, 0, 10, "2 Gijka - Gjika", -1);
	global_dpd_->buf4_close(&G);

	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 10, 0, 10, 0, "2 Gijka - Gjika");
	global_dpd_->buf4_init(&Cints, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
	global_dpd_->contract442(&Cints, &G, &I, 1, 1, 2.0, 1.0);
	global_dpd_->buf4_close(&Cints);
	global_dpd_->buf4_init(&Dints, PSIF_CC_DINTS, 0, 11, 10, 11, 10, 0, "D <ij|ab> (aj,ib)");
	global_dpd_->contract442(&Dints, &G, &I, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&G);

	global_dpd_->file2_close(&I);
      }
      else if(params.ref == 1) { /** ROHF **/

	/* I'AI <-- 2 sum_JKB <JA||KB> G(JI,KB) + 2 sum_jkB <jA|kB> G(jI,kB) +
	   2 sum_jKb <Kj|Ab> (Aj,Kb) G(Ij,Kb) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 0, "I'AI");

	global_dpd_->buf4_init(&Cints, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 10, 2, 10, 0, "GIJKA");
	global_dpd_->contract442(&Cints, &G, &I, 1, 1, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Cints);

	global_dpd_->buf4_init(&Cints, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 10, 0, 10, 0, "GiJkA");
	global_dpd_->contract442(&Cints, &G, &I, 1, 1, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Cints);

	global_dpd_->buf4_init(&Dints, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
	global_dpd_->buf4_sort(&Dints, PSIF_CC_TMP0, rqps, 11, 10, "D <ij|ab> (aj,ib)");
	global_dpd_->buf4_close(&Dints);
	global_dpd_->buf4_init(&Dints, PSIF_CC_TMP0, 0, 11, 10, 11, 10, 0, "D <ij|ab> (aj,ib)");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 10, 0, 10, 0, "GIjKa");
	global_dpd_->contract442(&Dints, &G, &I, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Dints);

	global_dpd_->file2_close(&I);

	/* I'ai <-- 2 sum_jkb <ja||kb> G(ji,kb) + 2 sum_JKb <jA|Kb> G(Ji,Kb) +
	   2 sum_JkB <kJ|aB> (aJ,kB) G(iJ,kB) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 0, "I'ai");

	global_dpd_->buf4_init(&Cints, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 10, 2, 10, 0, "Gijka");
	global_dpd_->contract442(&Cints, &G, &I, 1, 1, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Cints);

	global_dpd_->buf4_init(&Cints, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 10, 0, 10, 0, "GIjKa");
	global_dpd_->contract442(&Cints, &G, &I, 1, 1, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Cints);

	/* This sorted D-group is formed in the last code block */
	global_dpd_->buf4_init(&Dints, PSIF_CC_TMP0, 0, 11, 10, 11, 10, 0, "D <ij|ab> (aj,ib)");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 10, 0, 10, 0, "GiJkA");
	global_dpd_->contract442(&Dints, &G, &I, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Dints);

	global_dpd_->file2_close(&I);
      }
      else if(params.ref == 2) { /** UHF **/

	/* I'AI <-- 2 sum_JKB <JA||KB> G(JI,KB) + 2 sum_jkB <jA|kB> G(jI,kB) +
	   2 sum_jKb <Kj|Ab> (Aj,Kb) G(Ij,Kb) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 0, "I'AI");

	global_dpd_->buf4_init(&Cints, PSIF_CC_CINTS, 0, 20, 20, 20, 20, 0, "C <IA||JB>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 20, 2, 20, 0, "GIJKA");
	global_dpd_->contract442(&Cints, &G, &I, 1, 1, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Cints);

	global_dpd_->buf4_init(&Cints, PSIF_CC_CINTS, 0, 27, 27, 27, 27, 0, "C <iA|jB>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 23, 27, 23, 27, 0, "GiJkA");
	global_dpd_->contract442(&Cints, &G, &I, 1, 1, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Cints);

	global_dpd_->buf4_init(&Dints, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
	global_dpd_->buf4_sort(&Dints, PSIF_CC_TMP0, rqps, 26, 24, "D <Ij|Ab> (Aj,Ib)");
	global_dpd_->buf4_close(&Dints);
	global_dpd_->buf4_init(&Dints, PSIF_CC_TMP0, 0, 26, 24, 26, 24, 0, "D <Ij|Ab> (Aj,Ib)");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 22, 24, 22, 24, 0, "GIjKa");
	global_dpd_->contract442(&Dints, &G, &I, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Dints);

	global_dpd_->file2_close(&I);

	/* I'ai <-- 2 sum_jkb <ja||kb> G(ji,kb) + 2 sum_JKb <Ja|Kb> G(Ji,Kb) +
	   2 sum_JkB <kJ|aB> (aJ,kB) G(iJ,kB) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 3, 2, "I'ai");

	global_dpd_->buf4_init(&Cints, PSIF_CC_CINTS, 0, 30, 30, 30, 30, 0, "C <ia||jb>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 30, 12, 30, 0, "Gijka");
	global_dpd_->contract442(&Cints, &G, &I, 1, 1, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Cints);

	global_dpd_->buf4_init(&Cints, PSIF_CC_CINTS, 0, 24, 24, 24, 24, 0, "C <Ia|Jb>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 22, 24, 22, 24, 0, "GIjKa");
	global_dpd_->contract442(&Cints, &G, &I, 1, 1, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Cints);

	global_dpd_->buf4_init(&Dints, PSIF_CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");
	global_dpd_->buf4_sort(&Dints, PSIF_CC_TMP0, rqps, 25, 27, "D <iJ|aB> (aJ,iB)");
	global_dpd_->buf4_close(&Dints);
	global_dpd_->buf4_init(&Dints, PSIF_CC_TMP0, 0, 25, 27, 25, 27, 0, "D <iJ|aB> (aJ,iB)");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 23, 27, 23, 27, 0, "GiJkA");
	global_dpd_->contract442(&Dints, &G, &I, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Dints);

	global_dpd_->file2_close(&I);

      }



      if(params.ref == 0) { /** RHF **/
	/* I'AI <-- sum_BJK <JK||AB> G(JK,IB) + 2 sum_Jkb <Jk|Ab> G(Jk,Ib) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 0, "I'AI");

	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 10, 0, 10, 0, "GIjKa");
	global_dpd_->buf4_scmcopy(&G, PSIF_CC_GAMMA, "2 Gijka - Gjika", 2);
	global_dpd_->buf4_sort_axpy(&G, PSIF_CC_GAMMA, qprs, 0, 10, "2 Gijka - Gjika", -1);
	global_dpd_->buf4_close(&G);

	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 10, 0, 10, 0, "2 Gijka - Gjika");

	global_dpd_->buf4_init(&Dints, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
	global_dpd_->contract442(&Dints, &G, &I, 2, 2, 2.0, 1.0);
	global_dpd_->buf4_close(&Dints);
	global_dpd_->buf4_close(&G);

	global_dpd_->file2_close(&I);
      }
      else if(params.ref == 1) { /** ROHF **/

	/* I'AI <-- sum_BJK <JK||AB> G(JK,IB) + 2 sum_Jkb <Jk|Ab> G(Jk,Ib) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 0, "I'AI");

	global_dpd_->buf4_init(&Dints, PSIF_CC_DINTS, 0, 2, 5, 2, 5, 0, "D <ij||ab> (i>j,ab)");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 2, 10, 2, 10, 0, "GIJKA");
	global_dpd_->contract442(&Dints, &G, &I, 2, 2, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Dints);

	global_dpd_->buf4_init(&Dints, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 10, 0, 10, 0, "GIjKa");
	global_dpd_->contract442(&Dints, &G, &I, 2, 2, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Dints);

	global_dpd_->file2_close(&I);

	/* I'ai <-- sum_bjk <jk||ab> G(jk,ib) + 2 sum_jKB <jK|aB> G(jK,iB) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 0, "I'ai");

	global_dpd_->buf4_init(&Dints, PSIF_CC_DINTS, 0, 2, 5, 2, 5, 0, "D <ij||ab> (i>j,ab)");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 2, 10, 2, 10, 0, "Gijka");
	global_dpd_->contract442(&Dints, &G, &I, 2, 2, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Dints);

	global_dpd_->buf4_init(&Dints, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 10, 0, 10, 0, "GiJkA");
	global_dpd_->contract442(&Dints, &G, &I, 2, 2, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Dints);

	global_dpd_->file2_close(&I);
      }
      else if(params.ref == 2) { /** UHF **/

	/* I'AI <-- sum_BJK <JK||AB> G(JK,IB) + 2 sum_Jkb <Jk|Ab> G(Jk,Ib) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 0, "I'AI");

	global_dpd_->buf4_init(&Dints, PSIF_CC_DINTS, 0, 2, 5, 2, 5, 0, "D <IJ||AB> (I>J,AB)");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 2, 20, 2, 20, 0, "GIJKA");
	global_dpd_->contract442(&Dints, &G, &I, 2, 2, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Dints);

	global_dpd_->buf4_init(&Dints, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 22, 24, 22, 24, 0, "GIjKa");
	global_dpd_->contract442(&Dints, &G, &I, 2, 2, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Dints);

	global_dpd_->file2_close(&I);

	/* I'ai <-- sum_bjk <jk||ab> G(jk,ib) + 2 sum_jKB <jK|aB> G(jK,iB) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 3, 2, "I'ai");

	global_dpd_->buf4_init(&Dints, PSIF_CC_DINTS, 0, 12, 15, 12, 15, 0, "D <ij||ab> (i>j,ab)");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 12, 30, 12, 30, 0, "Gijka");
	global_dpd_->contract442(&Dints, &G, &I, 2, 2, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Dints);

	global_dpd_->buf4_init(&Dints, PSIF_CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 23, 27, 23, 27, 0, "GiJkA");
	global_dpd_->contract442(&Dints, &G, &I, 2, 2, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Dints);

	global_dpd_->file2_close(&I);

      }



      if(params.ref == 0) { /** RHF **/
	/* I'AI <-- sum_JBC <JA||BC> G(JI,BC) + 2 sum_jbC <jA|bC> G(jI,bC) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 0, "I'AI");

	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
	global_dpd_->buf4_scmcopy(&G, PSIF_CC_GAMMA, "2 Gijab - Gijba", 2);
	global_dpd_->buf4_sort_axpy(&G, PSIF_CC_GAMMA, pqsr, 0, 5, "2 Gijab - Gijba", -1);
	global_dpd_->buf4_close(&G);

	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 5, 0, 5, 0, "2 Gijab - Gijba");
	global_dpd_->buf4_init(&Fints, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
	global_dpd_->contract442(&Fints, &G, &I, 1, 1, 2.0, 1.0);
	global_dpd_->buf4_close(&Fints);
	global_dpd_->buf4_close(&G);

	global_dpd_->file2_close(&I);
      }
      else if(params.ref == 1) { /** ROHF **/

	/* I'AI <-- sum_JBC <JA||BC> G(JI,BC) + 2 sum_jbC <jA|bC> G(jI,bC) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 0, "I'AI");

	global_dpd_->buf4_init(&Fints, PSIF_CC_FINTS, 0, 10, 7, 10, 5, 1, "F <ia|bc>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 7, 2, 7, 0, "GIJAB");
	global_dpd_->contract442(&Fints, &G, &I, 1, 1, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Fints);

	global_dpd_->buf4_init(&Fints, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
	global_dpd_->buf4_sort(&G, PSIF_CC_TMP0, qpsr, 0, 5, "GiJaB");
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_init(&G, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "GiJaB");
	global_dpd_->contract442(&Fints, &G, &I, 1, 1, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Fints);

	global_dpd_->file2_close(&I);

	/* I'ai <-- sum_jbc <ja||bc> G(ji,bc) + 2 sum_JBc <Ja|Bc> G(Ji,Bc) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 0, "I'ai");

	global_dpd_->buf4_init(&Fints, PSIF_CC_FINTS, 0, 10, 7, 10, 5, 1, "F <ia|bc>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 7, 2, 7, 0, "Gijab");
	global_dpd_->contract442(&Fints, &G, &I, 1, 1, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Fints);

	global_dpd_->buf4_init(&Fints, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
	global_dpd_->contract442(&Fints, &G, &I, 1, 1, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Fints);

	global_dpd_->file2_close(&I);
      }
      else if(params.ref == 2) { /** UHF **/

	/* I'AI <-- sum_JBC <JA||BC> G(JI,BC) + 2 sum_jbC <jA|bC> G(jI,bC) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 0, "I'AI");

	global_dpd_->buf4_init(&Fints, PSIF_CC_FINTS, 0, 20, 7, 20, 5, 1, "F <IA|BC>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 7, 2, 7, 0, "GIJAB");
	global_dpd_->contract442(&Fints, &G, &I, 1, 1, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Fints);

	global_dpd_->buf4_init(&Fints, PSIF_CC_FINTS, 0, 27, 29, 27, 29, 0, "F <iA|bC>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 22, 28, 22, 28, 0, "GIjAb");
	global_dpd_->buf4_sort(&G, PSIF_CC_TMP0, qpsr, 23, 29, "GiJaB");
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_init(&G, PSIF_CC_TMP0, 0, 23, 29, 23, 29, 0, "GiJaB");
	global_dpd_->contract442(&Fints, &G, &I, 1, 1, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Fints);

	global_dpd_->file2_close(&I);

	/* I'ai <-- sum_jbc <ja||bc> G(ji,bc) + 2 sum_JBc <Ja|Bc> G(Ji,Bc) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 3, 2, "I'ai");

	global_dpd_->buf4_init(&Fints, PSIF_CC_FINTS, 0, 30, 17, 30, 15, 1, "F <ia|bc>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 17, 12, 17, 0, "Gijab");
	global_dpd_->contract442(&Fints, &G, &I, 1, 1, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Fints);

	global_dpd_->buf4_init(&Fints, PSIF_CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 22, 28, 22, 28, 0, "GIjAb");
	global_dpd_->contract442(&Fints, &G, &I, 1, 1, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Fints);

	global_dpd_->file2_close(&I);

      }



      if(params.ref == 0) { /** RHF **/
	/* I'AI <-- 2 sum_BJC <JC||AB> G(JC,IB) + 2 sum_bJc <Jc||Ab> G(Jc,Ib) +
	   2 sum_bjC <jC|bA> GjCIb(jC,bI) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 0, "I'AI");

	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIBJA");
	global_dpd_->buf4_init(&Fints, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 1, "F <ia|bc>");
	global_dpd_->contract442(&Fints, &G, &I, 2, 2, 2.0, 1.0);
	global_dpd_->buf4_close(&Fints);
	global_dpd_->buf4_close(&G);

	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbJa");
	global_dpd_->buf4_init(&Fints, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
	global_dpd_->contract442(&Fints, &G, &I, 2, 2, 2.0, 1.0);
	global_dpd_->buf4_close(&Fints);
	global_dpd_->buf4_close(&G);

	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbjA");
	global_dpd_->buf4_sort(&G, PSIF_CC_TMP0, rsqp, 10, 11, "GIbjA (jA,bI)");
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_init(&G, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "GIbjA (jA,bI)");
	global_dpd_->buf4_init(&Fints, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
	global_dpd_->contract442(&Fints, &G, &I, 3, 3, -2.0, 1.0);
	global_dpd_->buf4_close(&Fints);
	global_dpd_->buf4_close(&G);

	global_dpd_->file2_close(&I);
      }
      else if(params.ref == 1) { /** ROHF **/

	/* I'AI <-- 2 sum_BJC <JC||AB> G(JC,IB) + 2 sum_bJc <Jc||Ab> G(Jc,Ib) +
	   2 sum_bjC <jC|bA> GjCIb(jC,bI) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 0, "I'AI");

	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIBJA");
	global_dpd_->buf4_init(&Fints, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 1, "F <ia|bc>");
	global_dpd_->contract442(&Fints, &G, &I, 2, 2, 2.0, 1.0);
	global_dpd_->buf4_close(&Fints);
	global_dpd_->buf4_close(&G);

	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbJa");
	global_dpd_->buf4_init(&Fints, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
	global_dpd_->contract442(&Fints, &G, &I, 2, 2, 2.0, 1.0);
	global_dpd_->buf4_close(&Fints);
	global_dpd_->buf4_close(&G);

	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbjA");
	global_dpd_->buf4_sort(&G, PSIF_CC_TMP0, rsqp, 10, 11, "GIbjA (jA,bI)");
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_init(&G, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "GIbjA (jA,bI)");
	global_dpd_->buf4_init(&Fints, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
	global_dpd_->contract442(&Fints, &G, &I, 3, 3, -2.0, 1.0);
	global_dpd_->buf4_close(&Fints);
	global_dpd_->buf4_close(&G);

	global_dpd_->file2_close(&I);

	/* I'ai <-- 2 sum_bjc <jc||ab> G(jc,ib) + 2 sum_BjC <jC||aB> G(jC,iB) +
	   2 sum_BJc <Jc|Ba> GJciB(Jc,Bi) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 0, "I'ai");

	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "Gibja");
	global_dpd_->buf4_init(&Fints, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 1, "F <ia|bc>");
	global_dpd_->contract442(&Fints, &G, &I, 2, 2, 2.0, 1.0);
	global_dpd_->buf4_close(&Fints);
	global_dpd_->buf4_close(&G);

	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GiBjA");
	global_dpd_->buf4_init(&Fints, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
	global_dpd_->contract442(&Fints, &G, &I, 2, 2, 2.0, 1.0);
	global_dpd_->buf4_close(&Fints);
	global_dpd_->buf4_close(&G);

	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GiBJa");
	global_dpd_->buf4_sort(&G, PSIF_CC_TMP0, rsqp, 10, 11, "GiBJa (Ja,Bi)");
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_init(&G, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "GiBJa (Ja,Bi)");
	global_dpd_->buf4_init(&Fints, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
	global_dpd_->contract442(&Fints, &G, &I, 3, 3, -2.0, 1.0);
	global_dpd_->buf4_close(&Fints);
	global_dpd_->buf4_close(&G);

	global_dpd_->file2_close(&I);
      }
      else if(params.ref == 2) { /** UHF **/

	/* I'AI <-- 2 sum_BJC <JC||AB> G(JC,IB) + 2 sum_bJc <Jc||Ab> G(Jc,Ib) +
	   2 sum_bjC <jC|bA> GjCIb(jC,bI) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 0, "I'AI");

	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 20, 20, 20, 20, 0, "GIBJA");
	global_dpd_->buf4_init(&Fints, PSIF_CC_FINTS, 0, 20, 5, 20, 5, 1, "F <IA|BC>");
	global_dpd_->contract442(&Fints, &G, &I, 2, 2, 2.0, 1.0);
	global_dpd_->buf4_close(&Fints);
	global_dpd_->buf4_close(&G);

	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 24, 24, 24, 24, 0, "GIbJa");
	global_dpd_->buf4_init(&Fints, PSIF_CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
	global_dpd_->contract442(&Fints, &G, &I, 2, 2, 2.0, 1.0);
	global_dpd_->buf4_close(&Fints);
	global_dpd_->buf4_close(&G);

	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 24, 27, 24, 27, 0, "GIbjA");
	global_dpd_->buf4_sort(&G, PSIF_CC_TMP0, rsqp, 27, 25, "GIbjA (jA,bI)");
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_init(&G, PSIF_CC_TMP0, 0, 27, 25, 27, 25, 0, "GIbjA (jA,bI)");
	global_dpd_->buf4_init(&Fints, PSIF_CC_FINTS, 0, 27, 29, 27, 29, 0, "F <iA|bC>");
	global_dpd_->contract442(&Fints, &G, &I, 3, 3, -2.0, 1.0);
	global_dpd_->buf4_close(&Fints);
	global_dpd_->buf4_close(&G);

	global_dpd_->file2_close(&I);

	/* I'ai <-- 2 sum_bjc <jc||ab> G(jc,ib) + 2 sum_BjC <jC||aB> G(jC,iB) +
	   2 sum_BJc <Jc|Ba> GJciB(Jc,Bi) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 3, 2, "I'ai");

	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 30, 30, 30, 30, 0, "Gibja");
	global_dpd_->buf4_init(&Fints, PSIF_CC_FINTS, 0, 30, 15, 30, 15, 1, "F <ia|bc>");
	global_dpd_->contract442(&Fints, &G, &I, 2, 2, 2.0, 1.0);
	global_dpd_->buf4_close(&Fints);
	global_dpd_->buf4_close(&G);

	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 27, 27, 27, 27, 0, "GiBjA");
	global_dpd_->buf4_init(&Fints, PSIF_CC_FINTS, 0, 27, 29, 27, 29, 0, "F <iA|bC>");
	global_dpd_->contract442(&Fints, &G, &I, 2, 2, 2.0, 1.0);
	global_dpd_->buf4_close(&Fints);
	global_dpd_->buf4_close(&G);

	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 27, 24, 27, 24, 0, "GiBJa");
	global_dpd_->buf4_sort(&G, PSIF_CC_TMP0, rsqp, 24, 26, "GiBJa (Ja,Bi)");
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_init(&G, PSIF_CC_TMP0, 0, 24, 26, 24, 26, 0, "GiBJa (Ja,Bi)");
	global_dpd_->buf4_init(&Fints, PSIF_CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
	global_dpd_->contract442(&Fints, &G, &I, 3, 3, -2.0, 1.0);
	global_dpd_->buf4_close(&Fints);
	global_dpd_->buf4_close(&G);

	global_dpd_->file2_close(&I);

      }



      if(params.ref == 0) { /** RHF **/
	/* I'AI <-- sum_BCD <AB||CD> G(IB,CD) + 2 sum_bCd <Ab|Cd> G(Ib,Cd) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 0, "I'AI");

	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 11, 5, 11, 5, 0, "GCiAb");
	global_dpd_->buf4_scmcopy(&G, PSIF_CC_GAMMA, "2 Gciab - Gciba", 2);
	global_dpd_->buf4_sort_axpy(&G, PSIF_CC_GAMMA, pqsr, 11, 5, "2 Gciab - Gciba", -1);
	global_dpd_->buf4_close(&G);

	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 11, 5, 11, 5, 0, "2 Gciab - Gciba");
	global_dpd_->buf4_sort(&G, PSIF_CC_GAMMA, qpsr, 10, 5, "2 Giabc - Giacb");
	global_dpd_->buf4_close(&G);

	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 5, 10, 5, 0, "2 Giabc - Giacb");
	global_dpd_->buf4_init(&Bints, PSIF_CC_BINTS, 0, 5, 5, 5, 5, 0, "B <ab|cd>");
	global_dpd_->contract442(&Bints, &G, &I, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&Bints);
	global_dpd_->buf4_close(&G);

	global_dpd_->file2_close(&I);
      }
      else if(params.ref == 1) { /** ROHF **/

	/* I'AI <-- sum_BCD <AB||CD> G(IB,CD) + 2 sum_bCd <Ab|Cd> G(Ib,Cd) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 0, "I'AI");

	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 11, 7, 11, 7, 0, "GCIAB");
	global_dpd_->buf4_sort(&G, PSIF_CC_TMP0, qprs, 10, 7, "GICAB");
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_init(&G, PSIF_CC_TMP0, 0, 10, 7, 10, 7, 0, "GICAB");
	global_dpd_->buf4_init(&Bints, PSIF_CC_BINTS, 0, 5, 7, 5, 5, 1, "B <ab|cd>");
	global_dpd_->contract442(&Bints, &G, &I, 0, 0, -2.0, 1.0);
	global_dpd_->buf4_close(&Bints);
	global_dpd_->buf4_close(&G);

	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 11, 5, 11, 5, 0, "GcIaB");
	global_dpd_->buf4_sort(&G, PSIF_CC_TMP0, qpsr, 10, 5, "GIcAb");
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_init(&G, PSIF_CC_TMP0, 0, 10, 5, 10, 5, 0, "GIcAb");
	global_dpd_->buf4_init(&Bints, PSIF_CC_BINTS, 0, 5, 5, 5, 5, 0, "B <ab|cd>");
	global_dpd_->contract442(&Bints, &G, &I, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&Bints);
	global_dpd_->buf4_close(&G);

	global_dpd_->file2_close(&I);

	/* I'ai <-- sum_bcd <ab||cd> G(ib,cd) + 2 sum_BcD <aB|cD> G(iB,cD) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 0, "I'ai");

	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 11, 7, 11, 7, 0, "Gciab");
	global_dpd_->buf4_sort(&G, PSIF_CC_TMP0, qprs, 10, 7, "Gicab");
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_init(&G, PSIF_CC_TMP0, 0, 10, 7, 10, 7, 0, "Gicab");
	global_dpd_->buf4_init(&Bints, PSIF_CC_BINTS, 0, 5, 7, 5, 5, 1, "B <ab|cd>");
	global_dpd_->contract442(&Bints, &G, &I, 0, 0, -2.0, 1.0);
	global_dpd_->buf4_close(&Bints);
	global_dpd_->buf4_close(&G);

	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 11, 5, 11, 5, 0, "GCiAb");
	global_dpd_->buf4_sort(&G, PSIF_CC_TMP0, qpsr, 10, 5, "GiCaB");
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_init(&G, PSIF_CC_TMP0, 0, 10, 5, 10, 5, 0, "GiCaB");
	global_dpd_->buf4_init(&Bints, PSIF_CC_BINTS, 0, 5, 5, 5, 5, 0, "B <ab|cd>");
	global_dpd_->contract442(&Bints, &G, &I, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&Bints);
	global_dpd_->buf4_close(&G);

	global_dpd_->file2_close(&I);
      }
      else if(params.ref == 2) { /** UHF **/

	/* I'AI <-- sum_BCD <AB||CD> G(IB,CD) + 2 sum_bCd <Ab|Cd> G(Ib,Cd) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 0, "I'AI");

	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 21, 7, 21, 7, 0, "GCIAB");
	global_dpd_->buf4_sort(&G, PSIF_CC_TMP0, qprs, 20, 7, "GICAB");
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_init(&G, PSIF_CC_TMP0, 0, 20, 7, 20, 7, 0, "GICAB");
	global_dpd_->buf4_init(&Bints, PSIF_CC_BINTS, 0, 5, 7, 5, 5, 1, "B <AB|CD>");
	global_dpd_->contract442(&Bints, &G, &I, 0, 0, -2.0, 1.0);
	global_dpd_->buf4_close(&Bints);
	global_dpd_->buf4_close(&G);

	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 25, 29, 25, 29, 0, "GcIaB");
	global_dpd_->buf4_sort(&G, PSIF_CC_TMP0, qpsr, 24, 28, "GIcAb");
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_init(&G, PSIF_CC_TMP0, 0, 24, 28, 24, 28, 0, "GIcAb");
	global_dpd_->buf4_init(&Bints, PSIF_CC_BINTS, 0, 28, 28, 28, 28, 0, "B <Ab|Cd>");
	global_dpd_->contract442(&Bints, &G, &I, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&Bints);
	global_dpd_->buf4_close(&G);

	global_dpd_->file2_close(&I);

	/* I'ai <-- sum_bcd <ab||cd> G(ib,cd) + 2 sum_BcD <Cd|Ba> G(Cd,Bi) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 3, 2, "I'ai");

	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 31, 17, 31, 17, 0, "Gciab");
	global_dpd_->buf4_sort(&G, PSIF_CC_TMP0, qprs, 30, 17, "Gicab");
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_init(&G, PSIF_CC_TMP0, 0, 30, 17, 30, 17, 0, "Gicab");
	global_dpd_->buf4_init(&Bints, PSIF_CC_BINTS, 0, 15, 17, 15, 15, 1, "B <ab|cd>");
	global_dpd_->contract442(&Bints, &G, &I, 0, 0, -2.0, 1.0);
	global_dpd_->buf4_close(&Bints);
	global_dpd_->buf4_close(&G);

	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 26, 28, 26, 28, 0, "GCiAb");
	global_dpd_->buf4_sort(&G, PSIF_CC_TMP0, rspq, 28, 26, "GAbCi");
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_init(&G, PSIF_CC_TMP0, 0, 28, 26, 28, 26, 0, "GAbCi");
	global_dpd_->buf4_init(&Bints, PSIF_CC_BINTS, 0, 28, 28, 28, 28, 0, "B <Ab|Cd>");
	global_dpd_->contract442(&Bints, &G, &I, 3, 3, 2.0, 1.0);
	global_dpd_->buf4_close(&Bints);
	global_dpd_->buf4_close(&G);

	global_dpd_->file2_close(&I);

      }

    }

  }} // namespace psi::ccdensity
