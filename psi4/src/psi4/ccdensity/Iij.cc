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

    /* Iij(): Build the occupied-occupied block of the orbital Lagrangian
    ** using the expression given in lag.c.  Note that we include an
    ** addition term here referred to as the reference contribution,
    ** 2*fij.  This comes from the general spin-orbital SCF gradient
    ** expression and is present for all reference types (though for
    ** canonical unperturbed orbitals, only the diagonal elements
    ** contribute, of course).
    **
    ** Note that the code currently produces only the unique I_IJ terms,
    ** but the actual contractions still need to be spin-adapted for
    ** greater efficiency.
    **
    ** TDC, 2/2008
    */

    void Iij(struct RHO_Params rho_params)
    {
      dpdfile2 I, F, D;
      dpdbuf4 G, Aints, Fints, Dints, Cints, Eints;

      if(params.ref == 0) { /** RHF **/
	/* I'IJ <-- sum_K fIK (DJK + DKJ) + sum_A fIA (DJA + DAJ) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'IJ");

	global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 0, "fIJ");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, rho_params.DIJ_lbl);
	global_dpd_->contract222(&F, &D, &I, 0, 0, 1.0, 0.0);
	global_dpd_->contract222(&F, &D, &I, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);

	/* Add reference contribution: I'IJ <-- 2 fIJ */
	global_dpd_->file2_axpy(&F, &I, 2.0, 0);
	global_dpd_->file2_close(&F);

	global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 1, "fIA");
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

	/* I'IJ <-- sum_K fIK (DJK + DKJ) + sum_A fIA (DJA + DAJ) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'IJ");

	global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 0, "fIJ");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, rho_params.DIJ_lbl);
	global_dpd_->contract222(&F, &D, &I, 0, 0, 1.0, 0.0);
	global_dpd_->contract222(&F, &D, &I, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);

	/* Add reference contribution: I'IJ <-- 2 fIJ */
	global_dpd_->file2_axpy(&F, &I, 2.0, 0);
	global_dpd_->file2_close(&F);

	global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 1, "fIA");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
	global_dpd_->contract222(&F, &D, &I, 0, 0, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
	global_dpd_->contract222(&F, &D, &I, 0, 0, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->file2_close(&F);

	global_dpd_->file2_close(&I);

	/* I'ij <-- sum_k fik (Djk + Dkj)  + sum_a fia (Dja + Daj) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'ij");

	global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 0, "fij");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, rho_params.Dij_lbl);
	global_dpd_->contract222(&F, &D, &I, 0, 0, 1.0, 0.0);
	global_dpd_->contract222(&F, &D, &I, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);

	/* Add reference contribution: I'ij <-- 2 fij */
	global_dpd_->file2_axpy(&F, &I, 2.0, 0);
	global_dpd_->file2_close(&F);

	global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 1, "fia");
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

	/* I'IJ <-- sum_K fIK (DJK + DKJ) + sum_A fIA (DJA + DAJ) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'IJ");

	global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 0, "fIJ");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, rho_params.DIJ_lbl);
	global_dpd_->contract222(&F, &D, &I, 0, 0, 1.0, 0.0);
	global_dpd_->contract222(&F, &D, &I, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);

	/* Add reference contribution: I'IJ <-- 2 fIJ */
	global_dpd_->file2_axpy(&F, &I, 2.0, 0);
	global_dpd_->file2_close(&F);

	global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 1, "fIA");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
	global_dpd_->contract222(&F, &D, &I, 0, 0, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
	global_dpd_->contract222(&F, &D, &I, 0, 0, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->file2_close(&F);

	global_dpd_->file2_close(&I);

	/* I'ij <-- sum_k fik (Djk + Dkj)  + sum_a fia (Dja + Daj) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 2, 2, "I'ij");

	global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 2, 2, "fij");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 2, 2, rho_params.Dij_lbl);
	global_dpd_->contract222(&F, &D, &I, 0, 0, 1.0, 0.0);
	global_dpd_->contract222(&F, &D, &I, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);

	/* Add reference contribution: I'ij <-- 2 fij */
	global_dpd_->file2_axpy(&F, &I, 2.0, 0);
	global_dpd_->file2_close(&F);

	global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 2, 3, "fia");
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
	/* I'IJ <-- sum_KL <IK||JL> (D_KL + D_LK) + sum_kl <Ik|Jl> (D_kl + D_lk) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'IJ");

	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, rho_params.DIJ_lbl);
	global_dpd_->buf4_init(&Aints, PSIF_CC_AINTS, 0, 0, 0, 0, 0, 1, "A <ij|kl>");
	global_dpd_->dot24(&D, &Aints, &I, 0, 0, 1.0, 1.0);
	global_dpd_->dot24(&D, &Aints, &I, 1, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&Aints);
	global_dpd_->file2_close(&D);

	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, rho_params.Dij_lbl);
	global_dpd_->buf4_init(&Aints, PSIF_CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
	global_dpd_->dot24(&D, &Aints, &I, 0, 0, 1.0, 1.0);
	global_dpd_->dot24(&D, &Aints, &I, 1, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&Aints);
	global_dpd_->file2_close(&D);

	global_dpd_->file2_close(&I);
      }
      else if(params.ref == 1) { /** ROHF **/

	/* I'IJ <-- sum_KL <IK||JL> (D_KL + D_LK) + sum_kl <Ik|Jl> (D_kl + D_lk) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'IJ");

	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, rho_params.DIJ_lbl);
	global_dpd_->buf4_init(&Aints, PSIF_CC_AINTS, 0, 0, 0, 0, 0, 1, "A <ij|kl>");
	global_dpd_->dot24(&D, &Aints, &I, 0, 0, 1.0, 1.0);
	global_dpd_->dot24(&D, &Aints, &I, 1, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&Aints);
	global_dpd_->file2_close(&D);

	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, rho_params.Dij_lbl);
	global_dpd_->buf4_init(&Aints, PSIF_CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
	global_dpd_->dot24(&D, &Aints, &I, 0, 0, 1.0, 1.0);
	global_dpd_->dot24(&D, &Aints, &I, 1, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&Aints);
	global_dpd_->file2_close(&D);

	global_dpd_->file2_close(&I);

	/* I'ij <-- sum_kl <ik||jl> (D_kl + D_lk) + sum_KL <iK|jL> (D_KL + D_LK) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'ij");

	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, rho_params.Dij_lbl);
	global_dpd_->buf4_init(&Aints, PSIF_CC_AINTS, 0, 0, 0, 0, 0, 1, "A <ij|kl>");
	global_dpd_->dot24(&D, &Aints, &I, 0, 0, 1.0, 1.0);
	global_dpd_->dot24(&D, &Aints, &I, 1, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&Aints);
	global_dpd_->file2_close(&D);

	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, rho_params.DIJ_lbl);
	global_dpd_->buf4_init(&Aints, PSIF_CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
	global_dpd_->dot24(&D, &Aints, &I, 0, 0, 1.0, 1.0);
	global_dpd_->dot24(&D, &Aints, &I, 1, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&Aints);
	global_dpd_->file2_close(&D);

	global_dpd_->file2_close(&I);
      }
      else if(params.ref == 2) { /** UHF **/

	/* I'IJ <-- sum_KL <IK||JL> (D_KL + D_LK) + sum_kl <Ik|Jl> (D_kl + D_lk) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'IJ");

	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, rho_params.DIJ_lbl);
	global_dpd_->buf4_init(&Aints, PSIF_CC_AINTS, 0, 0, 0, 0, 0, 1, "A <IJ|KL>");
	global_dpd_->dot24(&D, &Aints, &I, 0, 0, 1.0, 1.0);
	global_dpd_->dot24(&D, &Aints, &I, 1, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&Aints);
	global_dpd_->file2_close(&D);

	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 2, 2, rho_params.Dij_lbl);
	global_dpd_->buf4_init(&Aints, PSIF_CC_AINTS, 0, 22, 22, 22, 22, 0, "A <Ij|Kl>");
	global_dpd_->dot24(&D, &Aints, &I, 0, 0, 1.0, 1.0);
	global_dpd_->dot24(&D, &Aints, &I, 1, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&Aints);
	global_dpd_->file2_close(&D);

	global_dpd_->file2_close(&I);

	/* I'ij <-- sum_kl <ik||jl> (D_kl + D_lk) + sum_KL <iK|jL> (D_KL + D_LK) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 2, 2, "I'ij");

	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 2, 2, rho_params.Dij_lbl);
	global_dpd_->buf4_init(&Aints, PSIF_CC_AINTS, 0, 10, 10, 10, 10, 1, "A <ij|kl>");
	global_dpd_->dot24(&D, &Aints, &I, 0, 0, 1.0, 1.0);
	global_dpd_->dot24(&D, &Aints, &I, 1, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&Aints);
	global_dpd_->file2_close(&D);

	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, rho_params.DIJ_lbl);
	global_dpd_->buf4_init(&Aints, PSIF_CC_AINTS, 0, 22, 22, 22, 22, 0, "A <Ij|Kl>");
	global_dpd_->dot13(&D, &Aints, &I, 0, 0, 1.0, 1.0);
	global_dpd_->dot13(&D, &Aints, &I, 1, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&Aints);
	global_dpd_->file2_close(&D);

	global_dpd_->file2_close(&I);

      }

      if(params.ref == 0) { /** RHF **/
	/* I'IJ <-- sum_KA <IK||JA> (D_KA + D_AK) + sum_ka <Ik|Ja> (D_ka + D_ak) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'IJ");

	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 0, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
	global_dpd_->dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
	global_dpd_->dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->buf4_close(&Eints);

	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.Dia_lbl);
	global_dpd_->dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.Dai_lbl);
	global_dpd_->dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->buf4_close(&Eints);

	global_dpd_->file2_close(&I);
      }
      else if(params.ref == 1) { /** ROHF **/

	/* I'IJ <-- sum_KA <IK||JA> (D_KA + D_AK) + sum_ka <Ik|Ja> (D_ka + D_ak) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'IJ");

	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 0, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
	global_dpd_->dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
	global_dpd_->dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->buf4_close(&Eints);

	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.Dia_lbl);
	global_dpd_->dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.Dai_lbl);
	global_dpd_->dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->buf4_close(&Eints);

	global_dpd_->file2_close(&I);

	/* I'ij <-- sum_ka <ik||ja> (D_ka + D_ak) + sum_KA <iK|jA> (D_KA + D_AK) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'ij");

	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 0, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.Dia_lbl);
	global_dpd_->dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.Dai_lbl);
	global_dpd_->dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->buf4_close(&Eints);

	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
	global_dpd_->dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
	global_dpd_->dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->buf4_close(&Eints);

	global_dpd_->file2_close(&I);
      }
      else if(params.ref == 2) { /** UHF **/

	/* I'IJ <-- sum_KA <IK||JA> (D_KA + D_AK) + sum_ka <Ik|Ja> (D_ka + D_ak) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'IJ");

	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 0, 20, 2, 20, 0, "E <IJ||KA> (I>J,KA)");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
	global_dpd_->dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
	global_dpd_->dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->buf4_close(&Eints);

	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 2, 3, rho_params.Dia_lbl);
	global_dpd_->dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 2, 3, rho_params.Dai_lbl);
	global_dpd_->dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->buf4_close(&Eints);

	global_dpd_->file2_close(&I);

	/* I'ij <-- sum_ka <ik||ja> (D_ka + D_ak) + sum_KA <iK|jA> (D_KA + D_AK) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 2, 2, "I'ij");

	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 10, 30, 12, 30, 0, "E <ij||ka> (i>j,ka)");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 2, 3, rho_params.Dia_lbl);
	global_dpd_->dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 2, 3, rho_params.Dai_lbl);
	global_dpd_->dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->buf4_close(&Eints);

	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 23, 27, 23, 27, 0, "E <iJ|kA>");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
	global_dpd_->dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
	global_dpd_->dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->buf4_close(&Eints);

	global_dpd_->file2_close(&I);
      }

      if(params.ref == 0) { /** RHF **/
	/* I'IJ <-- sum_AK <JK||IA> (D_AK + D_KA) + sum_ak <Jk|Ia> (D_ak + D_ka) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'IJ");

	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 0, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
	global_dpd_->dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
	global_dpd_->dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->buf4_close(&Eints);

	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.Dia_lbl);
	global_dpd_->dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.Dai_lbl);
	global_dpd_->dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->buf4_close(&Eints);

	global_dpd_->file2_close(&I);
      }
      else if(params.ref == 1) { /** ROHF **/

	/* I'IJ <-- sum_AK <JK||IA> (D_AK + D_KA) + sum_ak <Jk|Ia> (D_ak + D_ka) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'IJ");

	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 0, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
	global_dpd_->dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
	global_dpd_->dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->buf4_close(&Eints);

	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.Dia_lbl);
	global_dpd_->dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.Dai_lbl);
	global_dpd_->dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->buf4_close(&Eints);

	global_dpd_->file2_close(&I);

	/* I'ij <-- sum_ak <jk||ia> (D_ak + D_ka) + sum_AK <jK|iA> (D_AK + D_KA) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'ij");

	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 0, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.Dia_lbl);
	global_dpd_->dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.Dai_lbl);
	global_dpd_->dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->buf4_close(&Eints);

	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
	global_dpd_->dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
	global_dpd_->dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->buf4_close(&Eints);

	global_dpd_->file2_close(&I);
      }
      else if(params.ref == 2) { /** UHF **/

	/* I'IJ <-- sum_AK <JK||IA> (D_AK + D_KA) + sum_ak <Jk|Ia> (D_ak + D_ka) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'IJ");

	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 0, 20, 2, 20, 0, "E <IJ||KA> (I>J,KA)");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
	global_dpd_->dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
	global_dpd_->dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->buf4_close(&Eints);

	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 2, 3, rho_params.Dia_lbl);
	global_dpd_->dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 2, 3, rho_params.Dai_lbl);
	global_dpd_->dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->buf4_close(&Eints);

	global_dpd_->file2_close(&I);

	/* I'ij <-- sum_ak <jk||ia> (D_ak + D_ka) + sum_AK <jK|iA> (D_AK + D_KA) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 2, 2, "I'ij");

	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 10, 30, 12, 30, 0, "E <ij||ka> (i>j,ka)");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 2, 3, rho_params.Dia_lbl);
	global_dpd_->dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 2, 3, rho_params.Dai_lbl);
	global_dpd_->dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->buf4_close(&Eints);

	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 23, 27, 23, 27, 0, "E <iJ|kA>");
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
	global_dpd_->dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
	global_dpd_->dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0);
	global_dpd_->file2_close(&D);
	global_dpd_->buf4_close(&Eints);

	global_dpd_->file2_close(&I);

      }

      if(params.ref == 0) { /** RHF **/
	/* I'IJ <-- sum_AB <IA||JB> (D_AB + D_BA) + sum_ab <Ia|Jb> (D_ab + D_ba) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'IJ");

	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, rho_params.DAB_lbl);
	global_dpd_->buf4_init(&Cints, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
	global_dpd_->dot24(&D, &Cints, &I, 0, 0, 1.0, 1.0);
	global_dpd_->dot24(&D, &Cints, &I, 1, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&Cints);
	global_dpd_->file2_close(&D);

	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, rho_params.Dab_lbl);
	global_dpd_->buf4_init(&Cints, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
	global_dpd_->dot24(&D, &Cints, &I, 0, 0, 1.0, 1.0);
	global_dpd_->dot24(&D, &Cints, &I, 1, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&Cints);
	global_dpd_->file2_close(&D);

	global_dpd_->file2_close(&I);
      }
      else if(params.ref == 1) { /** ROHF **/

	/* I'IJ <-- sum_AB <IA||JB> (D_AB + D_BA) + sum_ab <Ia|Jb> (D_ab + D_ba) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'IJ");

	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, rho_params.DAB_lbl);
	global_dpd_->buf4_init(&Cints, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
	global_dpd_->dot24(&D, &Cints, &I, 0, 0, 1.0, 1.0);
	global_dpd_->dot24(&D, &Cints, &I, 1, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&Cints);
	global_dpd_->file2_close(&D);

	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, rho_params.Dab_lbl);
	global_dpd_->buf4_init(&Cints, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
	global_dpd_->dot24(&D, &Cints, &I, 0, 0, 1.0, 1.0);
	global_dpd_->dot24(&D, &Cints, &I, 1, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&Cints);
	global_dpd_->file2_close(&D);

	global_dpd_->file2_close(&I);

	/* I'ij <-- sum_ab <ia||jb> (D_ab + D_ba) + sum_AB <iA|jB> (D_AB + D_BA) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'ij");

	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, rho_params.Dab_lbl);
	global_dpd_->buf4_init(&Cints, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
	global_dpd_->dot24(&D, &Cints, &I, 0, 0, 1.0, 1.0);
	global_dpd_->dot24(&D, &Cints, &I, 1, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&Cints);
	global_dpd_->file2_close(&D);

	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, rho_params.DAB_lbl);
	global_dpd_->buf4_init(&Cints, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
	global_dpd_->dot24(&D, &Cints, &I, 0, 0, 1.0, 1.0);
	global_dpd_->dot24(&D, &Cints, &I, 1, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&Cints);
	global_dpd_->file2_close(&D);

	global_dpd_->file2_close(&I);
      }
      else if(params.ref == 2) { /** UHF **/

	/* I'IJ <-- sum_AB <IA||JB> (D_AB + D_BA) + sum_ab <Ia|Jb> (D_ab + D_ba) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'IJ");

	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, rho_params.DAB_lbl);
	global_dpd_->buf4_init(&Cints, PSIF_CC_CINTS, 0, 20, 20, 20, 20, 0, "C <IA||JB>");
	global_dpd_->dot24(&D, &Cints, &I, 0, 0, 1.0, 1.0);
	global_dpd_->dot24(&D, &Cints, &I, 1, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&Cints);
	global_dpd_->file2_close(&D);

	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 3, 3, rho_params.Dab_lbl);
	global_dpd_->buf4_init(&Cints, PSIF_CC_CINTS, 0, 24, 24, 24, 24, 0, "C <Ia|Jb>");
	global_dpd_->dot24(&D, &Cints, &I, 0, 0, 1.0, 1.0);
	global_dpd_->dot24(&D, &Cints, &I, 1, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&Cints);
	global_dpd_->file2_close(&D);

	global_dpd_->file2_close(&I);

	/* I'ij <-- sum_ab <ia||jb> (D_ab + D_ba) + sum_AB <iA|jB> (D_AB + D_BA) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 2, 2, "I'ij");

	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 3, 3, rho_params.Dab_lbl);
	global_dpd_->buf4_init(&Cints, PSIF_CC_CINTS, 0, 30, 30, 30, 30, 0, "C <ia||jb>");
	global_dpd_->dot24(&D, &Cints, &I, 0, 0, 1.0, 1.0);
	global_dpd_->dot24(&D, &Cints, &I, 1, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&Cints);
	global_dpd_->file2_close(&D);

	global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, rho_params.DAB_lbl);
	global_dpd_->buf4_init(&Cints, PSIF_CC_CINTS, 0, 27, 27, 27, 27, 0, "C <iA|jB>");
	global_dpd_->dot24(&D, &Cints, &I, 0, 0, 1.0, 1.0);
	global_dpd_->dot24(&D, &Cints, &I, 1, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&Cints);
	global_dpd_->file2_close(&D);

	global_dpd_->file2_close(&I);

      }

      if(params.ref == 0) { /** RHF **/
	/* I'IJ <-- 2 [ 2 (Gjklm - Gjkml) <ik|lm>] */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'IJ");

	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 0, 0, 0, 0, "GIjKl");
	global_dpd_->buf4_scmcopy(&G, PSIF_CC_GAMMA, "2 Gijkl - Gijlk", 2);
	global_dpd_->buf4_sort_axpy(&G, PSIF_CC_GAMMA, pqsr, 0, 0, "2 Gijkl - Gijlk", -1);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 0, 0, 0, 0, "2 Gijkl - Gijlk");
	global_dpd_->buf4_init(&Aints, PSIF_CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
	global_dpd_->contract442(&Aints, &G, &I, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&Aints);
	global_dpd_->buf4_close(&G);
	global_dpd_->file2_close(&I);
      }
      else if(params.ref == 1) { /** ROHF **/

	/* I'IJ <-- sum_KLM <IK||LM> G(JK,LM) + 2 sum_kLm <Ik|Lm> G(Jk,Lm) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'IJ");

	global_dpd_->buf4_init(&Aints, PSIF_CC_AINTS, 0, 0, 2, 0, 0, 1, "A <ij|kl>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 2, 2, 2, 0, "GIJKL");
	global_dpd_->contract442(&Aints, &G, &I, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Aints);

	global_dpd_->buf4_init(&Aints, PSIF_CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 0, 0, 0, 0, "GIjKl");
	global_dpd_->contract442(&Aints, &G, &I, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Aints);

	global_dpd_->file2_close(&I);

	/* I'ij <-- sum_klm <ik||lm> G(jk,lm) + 2 sum_KlM <Ki|Ml> G(Kj,Ml) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'ij");

	global_dpd_->buf4_init(&Aints, PSIF_CC_AINTS, 0, 0, 2, 0, 0, 1, "A <ij|kl>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 2, 2, 2, 0, "Gijkl");
	global_dpd_->contract442(&Aints, &G, &I, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Aints);

	global_dpd_->buf4_init(&Aints, PSIF_CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 0, 0, 0, 0, "GIjKl");
	global_dpd_->contract442(&Aints, &G, &I, 1, 1, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Aints);

	global_dpd_->file2_close(&I);

      }
      else if(params.ref == 2) { /** UHF **/

	/* I'IJ <-- sum_KLM <IK||LM> G(JK,LM) + 2 sum_kLm <Ik|Lm> G(Jk,Lm) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'IJ");

	global_dpd_->buf4_init(&Aints, PSIF_CC_AINTS, 0, 0, 2, 0, 0, 1, "A <IJ|KL>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 2, 2, 2, 0, "GIJKL");
	global_dpd_->contract442(&Aints, &G, &I, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Aints);

	global_dpd_->buf4_init(&Aints, PSIF_CC_AINTS, 0, 22, 22, 22, 22, 0, "A <Ij|Kl>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 22, 22, 22, 22, 0, "GIjKl");
	global_dpd_->contract442(&Aints, &G, &I, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Aints);

	global_dpd_->file2_close(&I);

	/* I'ij <-- sum_klm <ik||lm> G(jk,lm) + 2 sum_KlM <Ki|Ml> G(Kj,Ml) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 2, 2, "I'ij");

	global_dpd_->buf4_init(&Aints, PSIF_CC_AINTS, 0, 10, 12, 10, 10, 1, "A <ij|kl>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 12, 12, 12, 0, "Gijkl");
	global_dpd_->contract442(&Aints, &G, &I, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Aints);

	global_dpd_->buf4_init(&Aints, PSIF_CC_AINTS, 0, 22, 22, 22, 22, 0, "A <Ij|Kl>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 22, 22, 22, 22, 0, "GIjKl");
	global_dpd_->contract442(&Aints, &G, &I, 1, 1, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Aints);

	global_dpd_->file2_close(&I);

      }

      if(params.ref == 0) { /** RHF **/
	/* I'IJ <-- sum_ABC <IA||BC> G(JA,BC) + 2 sum_AbC <aI|bC> G(aJ,bC) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'IJ");

	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 11, 5, 11, 5, 0, "GCiAb");
	global_dpd_->buf4_scmcopy(&G, PSIF_CC_GAMMA, "2 Gciab - Gciba", 2);
	global_dpd_->buf4_sort_axpy(&G, PSIF_CC_GAMMA, pqsr, 11, 5, "2 Gciab - Gciba", -1);
	global_dpd_->buf4_close(&G);

	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 11, 5, 11, 5, 0, "2 Gciab - Gciba");
	global_dpd_->buf4_sort(&G, PSIF_CC_GAMMA, qpsr, 10, 5, "2 Giabc - Giacb");
	global_dpd_->buf4_close(&G);

	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 5, 10, 5, 0, "2 Giabc - Giacb");
	global_dpd_->buf4_init(&Fints, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
	global_dpd_->contract442(&Fints, &G, &I, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&Fints);
	global_dpd_->buf4_close(&G);

	global_dpd_->file2_close(&I);
      }
      else if(params.ref == 1) { /** ROHF **/

	/* I'IJ <-- sum_ABC <IA||BC> G(JA,BC) + 2 sum_AbC <aI|bC> G(aJ,bC) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'IJ");

	global_dpd_->buf4_init(&Fints, PSIF_CC_FINTS, 0, 10, 7, 10, 5, 1, "F <ia|bc>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 11, 7, 11, 7, 0, "GCIAB");
	global_dpd_->buf4_sort(&G, PSIF_CC_TMP0, qprs, 10, 7, "GICAB");
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_init(&G, PSIF_CC_TMP0, 0, 10, 7, 10, 7, 0, "GICAB");
	global_dpd_->contract442(&Fints, &G, &I, 0, 0, -2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Fints);

	global_dpd_->buf4_init(&Fints, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 11, 5, 11, 5, 0, "GcIaB");
	global_dpd_->buf4_sort(&G, PSIF_CC_TMP0, qpsr, 10, 5, "GIcBa");
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_init(&G, PSIF_CC_TMP0, 0, 10, 5, 10, 5, 0, "GIcBa");
	global_dpd_->contract442(&Fints, &G, &I, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Fints);

	global_dpd_->file2_close(&I);

	/* I'ij <-- sum_abc <ia||bcC> G(ja,bc) + 2 sum_AbC <Ai|Bc> G(Aj,Bc) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'ij");

	global_dpd_->buf4_init(&Fints, PSIF_CC_FINTS, 0, 10, 7, 10, 5, 1, "F <ia|bc>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 11, 7, 11, 7, 0, "Gciab");
	global_dpd_->buf4_sort(&G, PSIF_CC_TMP0, qprs, 10, 7, "Gicab");
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_init(&G, PSIF_CC_TMP0, 0, 10, 7, 10, 7, 0, "Gicab");
	global_dpd_->contract442(&Fints, &G, &I, 0, 0, -2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Fints);

	global_dpd_->buf4_init(&Fints, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 11, 5, 11, 5, 0, "GCiAb");
	global_dpd_->buf4_sort(&G, PSIF_CC_TMP0, qpsr, 10, 5, "GiCbA");
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_init(&G, PSIF_CC_TMP0, 0, 10, 5, 10, 5, 0, "GiCbA");
	global_dpd_->contract442(&Fints, &G, &I, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Fints);

	global_dpd_->file2_close(&I);

      }
      else if(params.ref == 2) { /** UHF **/

	/* I'IJ <-- sum_ABC <IA||BC> G(JA,BC) + 2 sum_AbC <aI|bC> G(aJ,bC) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'IJ");

	global_dpd_->buf4_init(&Fints, PSIF_CC_FINTS, 0, 20, 7, 20, 5, 1, "F <IA|BC>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 21, 7, 21, 7, 0, "GCIAB");
	global_dpd_->buf4_sort(&G, PSIF_CC_TMP0, qprs, 20, 7, "GICAB");
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_init(&G, PSIF_CC_TMP0, 0, 20, 7, 20, 7, 0, "GICAB");
	global_dpd_->contract442(&Fints, &G, &I, 0, 0, -2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Fints);

	global_dpd_->buf4_init(&Fints, PSIF_CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 25, 29, 25, 29, 0, "GcIaB");
	global_dpd_->buf4_sort(&G, PSIF_CC_TMP0, qpsr, 24, 28, "GIcBa");
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_init(&G, PSIF_CC_TMP0, 0, 24, 28, 24, 28, 0, "GIcBa");
	global_dpd_->contract442(&Fints, &G, &I, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Fints);

	global_dpd_->file2_close(&I);

	/* I'ij <-- sum_abc <ia||bc> G(ja,bc) + 2 sum_AbC <Ai|Bc> G(Aj,Bc) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 2, 2, "I'ij");

	global_dpd_->buf4_init(&Fints, PSIF_CC_FINTS, 0, 30, 17, 30, 15, 1, "F <ia|bc>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 31, 17, 31, 17, 0, "Gciab");
	global_dpd_->buf4_sort(&G, PSIF_CC_TMP0, qprs, 30, 17, "Gicab");
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_init(&G, PSIF_CC_TMP0, 0, 30, 17, 30, 17, 0, "Gicab");
	global_dpd_->contract442(&Fints, &G, &I, 0, 0, -2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Fints);

	global_dpd_->buf4_init(&Fints, PSIF_CC_FINTS, 0, 27, 29, 27, 29, 0, "F <iA|bC>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 26, 28, 26, 28, 0, "GCiAb");
	global_dpd_->buf4_sort(&G, PSIF_CC_TMP0, qpsr, 27, 29, "GiCbA");
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_init(&G, PSIF_CC_TMP0, 0, 27, 29, 27, 29, 0, "GiCbA");
	global_dpd_->contract442(&Fints, &G, &I, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Fints);

	global_dpd_->file2_close(&I);

      }

      if(params.ref == 0) { /** RHF **/
	/* I'IJ <-- sum_KAB <IK||AB> G(JK,AB) + 2 sum_kAb <Ik|Ab> G(Jk,Ab) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'IJ");

	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
	global_dpd_->buf4_scmcopy(&G, PSIF_CC_GAMMA, "2 Gijab - Gijba", 2);
	global_dpd_->buf4_sort_axpy(&G, PSIF_CC_GAMMA, pqsr, 0, 5, "2 Gijab - Gijba", -1);
	global_dpd_->buf4_close(&G);

	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 5, 0, 5, 0, "2 Gijab - Gijba");
	global_dpd_->buf4_init(&Dints, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
	global_dpd_->contract442(&Dints, &G, &I, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&Dints);
	global_dpd_->buf4_close(&G);

	global_dpd_->file2_close(&I);
      }
      else if(params.ref == 1) { /** ROHF **/

	/* I'IJ <-- sum_KAB <IK||AB> G(JK,AB) + 2 sum_kAb <Ik|Ab> G(Jk,Ab) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'IJ");

	global_dpd_->buf4_init(&Dints, PSIF_CC_DINTS, 0, 0, 7, 0, 7, 0, "D <ij||ab> (ij,a>b)");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 7, 2, 7, 0, "GIJAB");
	global_dpd_->contract442(&Dints, &G, &I, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Dints);

	global_dpd_->buf4_init(&Dints, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
	global_dpd_->contract442(&Dints, &G, &I, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Dints);

	global_dpd_->file2_close(&I);

	/* I'ij <-- sum_kab <ik||ab> G(jk,ab) + 2 sum_KaB <iK|aB> G(jK,aB) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'ij");

	global_dpd_->buf4_init(&Dints, PSIF_CC_DINTS, 0, 0, 7, 0, 7, 0, "D <ij||ab> (ij,a>b)");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 7, 2, 7, 0, "Gijab");
	global_dpd_->contract442(&Dints, &G, &I, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Dints);

	global_dpd_->buf4_init(&Dints, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
	global_dpd_->contract442(&Dints, &G, &I, 1, 1, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Dints);

	global_dpd_->file2_close(&I);

      }
      else if(params.ref == 2) { /** UHF **/

	/* I'IJ <-- sum_KAB <IK||AB> G(JK,AB) + 2 sum_kAb <Ik|Ab> G(Jk,Ab) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'IJ");

	global_dpd_->buf4_init(&Dints, PSIF_CC_DINTS, 0, 0, 7, 0, 7, 0, "D <IJ||AB> (IJ,A>B)");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 7, 2, 7, 0, "GIJAB");
	global_dpd_->contract442(&Dints, &G, &I, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Dints);

	global_dpd_->buf4_init(&Dints, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 22, 28, 22, 28, 0, "GIjAb");
	global_dpd_->contract442(&Dints, &G, &I, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Dints);

	global_dpd_->file2_close(&I);

	/* I'ij <-- sum_kab <ik||ab> G(jk,ab) + 2 sum_KaB <iK|aB> G(jK,aB) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 2, 2, "I'ij");

	global_dpd_->buf4_init(&Dints, PSIF_CC_DINTS, 0, 10, 17, 10, 17, 0, "D <ij||ab> (ij,a>b)");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 17, 12, 17, 0, "Gijab");
	global_dpd_->contract442(&Dints, &G, &I, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Dints);

	global_dpd_->buf4_init(&Dints, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 22, 28, 22, 28, 0, "GIjAb");
	global_dpd_->contract442(&Dints, &G, &I, 1, 1, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Dints);

	global_dpd_->file2_close(&I);

      }

      if(params.ref == 0) { /** RHF **/
	/* I'IJ <-- 2 sum_AKB <IA||KB> G(JA,KB) + 2 sum_aKb <Ia|Kb> G(Ja,Kb) -
	   2 sum_akB <Ik|Ba> GJakB(Jk,Ba) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'IJ");

	global_dpd_->buf4_init(&Cints, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIBJA");
	global_dpd_->contract442(&Cints, &G, &I, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Cints);

	global_dpd_->buf4_init(&Cints, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbJa");
	global_dpd_->contract442(&Cints, &G, &I, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Cints);

	global_dpd_->buf4_init(&Dints, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbjA");
	global_dpd_->buf4_sort(&G, PSIF_CC_TMP0, prsq, 0, 5, "GIbjA (Ij,Ab)");
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_init(&G, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "GIbjA (Ij,Ab)");
	global_dpd_->contract442(&Dints, &G, &I, 0, 0, -2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Dints);

	global_dpd_->file2_close(&I);
      }
      else if(params.ref == 1) { /** ROHF **/

	/* I'IJ <-- 2 sum_AKB <IA||KB> G(JA,KB) + 2 sum_aKb <Ia|Kb> G(Ja,Kb) -
	   2 sum_akB <Ik|Ba> GJakB(Jk,Ba) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'IJ");

	global_dpd_->buf4_init(&Cints, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIBJA");
	global_dpd_->contract442(&Cints, &G, &I, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Cints);

	global_dpd_->buf4_init(&Cints, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbJa");
	global_dpd_->contract442(&Cints, &G, &I, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Cints);

	global_dpd_->buf4_init(&Dints, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbjA");
	global_dpd_->buf4_sort(&G, PSIF_CC_TMP0, prsq, 0, 5, "GIbjA (Ij,Ab)");
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_init(&G, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "GIbjA (Ij,Ab)");
	global_dpd_->contract442(&Dints, &G, &I, 0, 0, -2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Dints);

	global_dpd_->file2_close(&I);

	/* I'ij <-- 2 sum_akb <ia||kb> G(ja,kb) + 2 sum_AkB <iA|kB> G(jA,kB) +
	   2 sum_AKb <iK|bA> GjAKb(jK,bA) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'ij");

	global_dpd_->buf4_init(&Cints, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "Gibja");
	global_dpd_->contract442(&Cints, &G, &I, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Cints);

	global_dpd_->buf4_init(&Cints, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GiBjA");
	global_dpd_->contract442(&Cints, &G, &I, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Cints);

	global_dpd_->buf4_init(&Dints, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GiBJa");
	global_dpd_->buf4_sort(&G, PSIF_CC_TMP0, prsq, 0, 5, "GiBJa (iJ,aB)");
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_init(&G, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "GiBJa (iJ,aB)");
	global_dpd_->contract442(&Dints, &G, &I, 0, 0, -2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Dints);

	global_dpd_->file2_close(&I);

      }
      else if(params.ref == 2) { /** UHF **/

	/* I'IJ <-- 2 sum_AKB <IA||KB> G(JA,KB) + 2 sum_aKb <Ia|Kb> G(Ja,Kb) -
	   2 sum_akB <Ik|Ba> GJakB(Jk,Ba) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'IJ");

	global_dpd_->buf4_init(&Cints, PSIF_CC_CINTS, 0, 20, 20, 20, 20, 0, "C <IA||JB>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 20, 20, 20, 20, 0, "GIBJA");
	global_dpd_->contract442(&Cints, &G, &I, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Cints);

	global_dpd_->buf4_init(&Cints, PSIF_CC_CINTS, 0, 24, 24, 24, 24, 0, "C <Ia|Jb>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 24, 24, 24, 24, 0, "GIbJa");
	global_dpd_->contract442(&Cints, &G, &I, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Cints);

	global_dpd_->buf4_init(&Dints, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 24, 27, 24, 27, 0, "GIbjA");
	global_dpd_->buf4_sort(&G, PSIF_CC_TMP0, prsq, 22, 28, "GIbjA (Ij,Ab)");
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_init(&G, PSIF_CC_TMP0, 0, 22, 28, 22, 28, 0, "GIbjA (Ij,Ab)");
	global_dpd_->contract442(&Dints, &G, &I, 0, 0, -2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Dints);

	global_dpd_->file2_close(&I);

	/* I'ij <-- 2 sum_akb <ia||kb> G(ja,kb) + 2 sum_AkB <iA|kB> G(jA,kB) +
	   2 sum_AKb <iK|bA> GjAKb(jK,bA) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 2, 2, "I'ij");

	global_dpd_->buf4_init(&Cints, PSIF_CC_CINTS, 0, 30, 30, 30, 30, 0, "C <ia||jb>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 30, 30, 30, 30, 0, "Gibja");
	global_dpd_->contract442(&Cints, &G, &I, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Cints);

	global_dpd_->buf4_init(&Cints, PSIF_CC_CINTS, 0, 27, 27, 27, 27, 0, "C <iA|jB>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 27, 27, 27, 27, 0, "GiBjA");
	global_dpd_->contract442(&Cints, &G, &I, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Cints);

	global_dpd_->buf4_init(&Dints, PSIF_CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 27, 24, 27, 24, 0, "GiBJa");
	global_dpd_->buf4_sort(&G, PSIF_CC_TMP0, prsq, 23, 29, "GiBJa (iJ,aB)");
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_init(&G, PSIF_CC_TMP0, 0, 23, 29, 23, 29, 0, "GiBJa (iJ,aB)");
	global_dpd_->contract442(&Dints, &G, &I, 0, 0, -2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Dints);

	global_dpd_->file2_close(&I);

      }

      if(params.ref == 0) { /** RHF **/
	/* I'IJ <-- 2 sum_KLA <IK||LA> G(JK,LA) + 2 sum_kLa <Ik|La> G(Jk,La)
	   + 2 sum_kAl <kI|lA> G(kJ,lA) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'IJ");

	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 10, 0, 10, 0, "GIjKa");
	global_dpd_->buf4_scmcopy(&G, PSIF_CC_GAMMA, "2 Gijka - Gjika", 2);
	global_dpd_->buf4_sort_axpy(&G, PSIF_CC_GAMMA, qprs, 0, 10, "2 Gijka - Gjika", -1);
	global_dpd_->buf4_close(&G);

	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 10, 0, 10, 0, "2 Gijka - Gjika");
	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
	global_dpd_->contract442(&Eints, &G, &I, 0, 0, 2.0, 1.0);
	global_dpd_->contract442(&Eints, &G, &I, 1, 1, 2.0, 1.0);
	global_dpd_->buf4_close(&Eints);
	global_dpd_->buf4_close(&G);

	global_dpd_->file2_close(&I);
      }
      else if(params.ref == 1) { /** ROHF **/

	/* I'IJ <-- 2 sum_KLA <IK||LA> G(JK,LA) + 2 sum_kLa <Ik|La> G(Jk,La)
	   + 2 sum_kAl <kI|lA> G(kJ,lA) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'IJ");

	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 0, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 10, 2, 10, 0, "GIJKA");
	global_dpd_->contract442(&Eints, &G, &I, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Eints);

	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 10, 0, 10, 0, "GIjKa");
	global_dpd_->contract442(&Eints, &G, &I, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Eints);

	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 10, 0, 10, 0, "GiJkA");
	global_dpd_->contract442(&Eints, &G, &I, 1, 1, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Eints);

	global_dpd_->file2_close(&I);

	/* I'ij <-- 2 sum_kla <ik||la> G(jk,la) + 2 sum_KlA <iK|lA> G(jK,lA)
	   + 2 sum_KaL <Ki|La> G(Kj,La) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'ij");

	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 0, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 10, 2, 10, 0, "Gijka");
	global_dpd_->contract442(&Eints, &G, &I, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Eints);

	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 10, 0, 10, 0, "GiJkA");
	global_dpd_->contract442(&Eints, &G, &I, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Eints);

	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 10, 0, 10, 0, "GIjKa");
	global_dpd_->contract442(&Eints, &G, &I, 1, 1, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Eints);

	global_dpd_->file2_close(&I);
      }
      else if(params.ref == 2) { /** UHF **/

	/* I'IJ <-- 2 sum_KLA <IK||LA> G(JK,LA) + 2 sum_kLa <Ik|La> G(Jk,La)
	   + 2 sum_kAl <kI|lA> G(kJ,lA) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'IJ");

	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 0, 20, 2, 20, 0, "E <IJ||KA> (I>J,KA)");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 20, 2, 20, 0, "GIJKA");
	global_dpd_->contract442(&Eints, &G, &I, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Eints);

	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 22, 24, 22, 24, 0, "GIjKa");
	global_dpd_->contract442(&Eints, &G, &I, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Eints);

	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 23, 27, 23, 27, 0, "E <iJ|kA>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 23, 27, 23, 27, 0, "GiJkA");
	global_dpd_->contract442(&Eints, &G, &I, 1, 1, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Eints);

	global_dpd_->file2_close(&I);

	/* I'ij <-- 2 sum_kla <ik||la> G(jk,la) + 2 sum_KlA <iK|lA> G(jK,lA)
	   + 2 sum_KaL <Ki|La> G(Kj,La) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 2, 2, "I'ij");

	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 10, 30, 12, 30, 0, "E <ij||ka> (i>j,ka)");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 30, 12, 30, 0, "Gijka");
	global_dpd_->contract442(&Eints, &G, &I, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Eints);

	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 23, 27, 23, 27, 0, "E <iJ|kA>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 23, 27, 23, 27, 0, "GiJkA");
	global_dpd_->contract442(&Eints, &G, &I, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Eints);

	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 22, 24, 22, 24, 0, "GIjKa");
	global_dpd_->contract442(&Eints, &G, &I, 1, 1, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Eints);

	global_dpd_->file2_close(&I);

      }

      if(params.ref == 0) { /** RHF **/
	/* I'IJ <-- 2 sum_AKL <K>L||IA> G(K>L,JA) + 2 sum_aKl <Kl|Ia> G(Kl,Ja) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'IJ");

	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 10, 0, 10, 0, "GIjKa");
	global_dpd_->buf4_scmcopy(&G, PSIF_CC_GAMMA, "2 Gijka - Gjika", 2);
	global_dpd_->buf4_sort_axpy(&G, PSIF_CC_GAMMA, qprs, 0, 10, "2 Gijka - Gjika", -1);
	global_dpd_->buf4_close(&G);

	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 10, 0, 10, 0, "2 Gijka - Gjika");
	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
	global_dpd_->contract442(&Eints, &G, &I, 2, 2, 2.0, 1.0);
	global_dpd_->buf4_close(&Eints);
	global_dpd_->buf4_close(&G);

	global_dpd_->file2_close(&I);
      }
      else if(params.ref == 1) { /** ROHF **/

	/* I'IJ <-- 2 sum_AKL <K>L||IA> G(K>L,JA) + 2 sum_aKl <Kl|Ia> G(Kl,Ja) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'IJ");

	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 2, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 2, 10, 2, 10, 0, "GIJKA");
	global_dpd_->contract442(&Eints, &G, &I, 2, 2, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Eints);

	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 10, 0, 10, 0, "GIjKa");
	global_dpd_->contract442(&Eints, &G, &I, 2, 2, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Eints);

	global_dpd_->file2_close(&I);

	/* I'ij <-- 2 sum_akl <k>l||ia> G(k>l,ja) + 2 sum_AkL <kL|iA> G(kL,jA) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'ij");

	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 2, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 2, 10, 2, 10, 0, "Gijka");
	global_dpd_->contract442(&Eints, &G, &I, 2, 2, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Eints);

	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 10, 0, 10, 0, "GiJkA");
	global_dpd_->contract442(&Eints, &G, &I, 2, 2, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Eints);

	global_dpd_->file2_close(&I);

      }
      else if(params.ref == 2) { /** UHF **/

	/* I'IJ <-- 2 sum_AKL <K>L||IA> G(K>L,JA) + 2 sum_aKl <Kl|Ia> G(Kl,Ja) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'IJ");

	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 2, 20, 2, 20, 0, "E <IJ||KA> (I>J,KA)");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 2, 20, 2, 20, 0, "GIJKA");
	global_dpd_->contract442(&Eints, &G, &I, 2, 2, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Eints);

	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 22, 24, 22, 24, 0, "GIjKa");
	global_dpd_->contract442(&Eints, &G, &I, 2, 2, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Eints);

	global_dpd_->file2_close(&I);

	/* I'ij <-- 2 sum_akl <k>l||ia> G(k>l,ja) + 2 sum_AkL <kL|iA> G(kL,jA) */
	global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 2, 2, "I'ij");

	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 12, 30, 12, 30, 0, "E <ij||ka> (i>j,ka)");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 12, 30, 12, 30, 0, "Gijka");
	global_dpd_->contract442(&Eints, &G, &I, 2, 2, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Eints);

	global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 23, 27, 23, 27, 0, "E <iJ|kA>");
	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 23, 27, 23, 27, 0, "GiJkA");
	global_dpd_->contract442(&Eints, &G, &I, 2, 2, 2.0, 1.0);
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_close(&Eints);

	global_dpd_->file2_close(&I);

      }
    }

  }} // namespace psi::ccdensity
