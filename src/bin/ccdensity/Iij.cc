/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here 
*/
#include <libdpd/dpd.h>
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
	dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'IJ");

	dpd_file2_init(&F, CC_OEI, 0, 0, 0, "fIJ");
	dpd_file2_init(&D, CC_OEI, 0, 0, 0, rho_params.DIJ_lbl);
	dpd_contract222(&F, &D, &I, 0, 0, 1.0, 0.0);
	dpd_contract222(&F, &D, &I, 0, 1, 1.0, 1.0);
	dpd_file2_close(&D);

	/* Add reference contribution: I'IJ <-- 2 fIJ */
	dpd_file2_axpy(&F, &I, 2.0, 0);
	dpd_file2_close(&F);

	dpd_file2_init(&F, CC_OEI, 0, 0, 1, "fIA");
	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
	dpd_contract222(&F, &D, &I, 0, 0, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
	dpd_contract222(&F, &D, &I, 0, 0, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_file2_close(&F);

	dpd_file2_close(&I);
      }
      else if(params.ref == 1) { /** ROHF **/

	/* I'IJ <-- sum_K fIK (DJK + DKJ) + sum_A fIA (DJA + DAJ) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'IJ");

	dpd_file2_init(&F, CC_OEI, 0, 0, 0, "fIJ");
	dpd_file2_init(&D, CC_OEI, 0, 0, 0, rho_params.DIJ_lbl);
	dpd_contract222(&F, &D, &I, 0, 0, 1.0, 0.0);
	dpd_contract222(&F, &D, &I, 0, 1, 1.0, 1.0);
	dpd_file2_close(&D);

	/* Add reference contribution: I'IJ <-- 2 fIJ */
	dpd_file2_axpy(&F, &I, 2.0, 0);
	dpd_file2_close(&F);

	dpd_file2_init(&F, CC_OEI, 0, 0, 1, "fIA");
	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
	dpd_contract222(&F, &D, &I, 0, 0, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
	dpd_contract222(&F, &D, &I, 0, 0, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_file2_close(&F);

	dpd_file2_close(&I);

	/* I'ij <-- sum_k fik (Djk + Dkj)  + sum_a fia (Dja + Daj) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'ij");

	dpd_file2_init(&F, CC_OEI, 0, 0, 0, "fij");
	dpd_file2_init(&D, CC_OEI, 0, 0, 0, rho_params.Dij_lbl);
	dpd_contract222(&F, &D, &I, 0, 0, 1.0, 0.0);
	dpd_contract222(&F, &D, &I, 0, 1, 1.0, 1.0);
	dpd_file2_close(&D);

	/* Add reference contribution: I'ij <-- 2 fij */
	dpd_file2_axpy(&F, &I, 2.0, 0);
	dpd_file2_close(&F);

	dpd_file2_init(&F, CC_OEI, 0, 0, 1, "fia");
	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.Dia_lbl);
	dpd_contract222(&F, &D, &I, 0, 0, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.Dai_lbl);
	dpd_contract222(&F, &D, &I, 0, 0, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_file2_close(&F);

	dpd_file2_close(&I);
      }
      else if(params.ref == 2) { /** UHF **/

	/* I'IJ <-- sum_K fIK (DJK + DKJ) + sum_A fIA (DJA + DAJ) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'IJ");

	dpd_file2_init(&F, CC_OEI, 0, 0, 0, "fIJ");
	dpd_file2_init(&D, CC_OEI, 0, 0, 0, rho_params.DIJ_lbl);
	dpd_contract222(&F, &D, &I, 0, 0, 1.0, 0.0);
	dpd_contract222(&F, &D, &I, 0, 1, 1.0, 1.0);
	dpd_file2_close(&D);

	/* Add reference contribution: I'IJ <-- 2 fIJ */
	dpd_file2_axpy(&F, &I, 2.0, 0);
	dpd_file2_close(&F);

	dpd_file2_init(&F, CC_OEI, 0, 0, 1, "fIA");
	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
	dpd_contract222(&F, &D, &I, 0, 0, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
	dpd_contract222(&F, &D, &I, 0, 0, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_file2_close(&F);

	dpd_file2_close(&I);

	/* I'ij <-- sum_k fik (Djk + Dkj)  + sum_a fia (Dja + Daj) */
	dpd_file2_init(&I, CC_OEI, 0, 2, 2, "I'ij");

	dpd_file2_init(&F, CC_OEI, 0, 2, 2, "fij");
	dpd_file2_init(&D, CC_OEI, 0, 2, 2, rho_params.Dij_lbl);
	dpd_contract222(&F, &D, &I, 0, 0, 1.0, 0.0);
	dpd_contract222(&F, &D, &I, 0, 1, 1.0, 1.0);
	dpd_file2_close(&D);

	/* Add reference contribution: I'ij <-- 2 fij */
	dpd_file2_axpy(&F, &I, 2.0, 0);
	dpd_file2_close(&F);

	dpd_file2_init(&F, CC_OEI, 0, 2, 3, "fia");
	dpd_file2_init(&D, CC_OEI, 0, 2, 3, rho_params.Dia_lbl);
	dpd_contract222(&F, &D, &I, 0, 0, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_file2_init(&D, CC_OEI, 0, 2, 3, rho_params.Dai_lbl);
	dpd_contract222(&F, &D, &I, 0, 0, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_file2_close(&F);

	dpd_file2_close(&I);
      }

      if(params.ref == 0) { /** RHF **/
	/* I'IJ <-- sum_KL <IK||JL> (D_KL + D_LK) + sum_kl <Ik|Jl> (D_kl + D_lk) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'IJ");

	dpd_file2_init(&D, CC_OEI, 0, 0, 0, rho_params.DIJ_lbl);
	dpd_buf4_init(&Aints, CC_AINTS, 0, 0, 0, 0, 0, 1, "A <ij|kl>");
	dpd_dot24(&D, &Aints, &I, 0, 0, 1.0, 1.0);
	dpd_dot24(&D, &Aints, &I, 1, 0, 1.0, 1.0);
	dpd_buf4_close(&Aints);
	dpd_file2_close(&D);

	dpd_file2_init(&D, CC_OEI, 0, 0, 0, rho_params.Dij_lbl);
	dpd_buf4_init(&Aints, CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
	dpd_dot24(&D, &Aints, &I, 0, 0, 1.0, 1.0);
	dpd_dot24(&D, &Aints, &I, 1, 0, 1.0, 1.0);
	dpd_buf4_close(&Aints);
	dpd_file2_close(&D);
  
	dpd_file2_close(&I);
      }
      else if(params.ref == 1) { /** ROHF **/

	/* I'IJ <-- sum_KL <IK||JL> (D_KL + D_LK) + sum_kl <Ik|Jl> (D_kl + D_lk) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'IJ");

	dpd_file2_init(&D, CC_OEI, 0, 0, 0, rho_params.DIJ_lbl);
	dpd_buf4_init(&Aints, CC_AINTS, 0, 0, 0, 0, 0, 1, "A <ij|kl>");
	dpd_dot24(&D, &Aints, &I, 0, 0, 1.0, 1.0);
	dpd_dot24(&D, &Aints, &I, 1, 0, 1.0, 1.0);
	dpd_buf4_close(&Aints);
	dpd_file2_close(&D);

	dpd_file2_init(&D, CC_OEI, 0, 0, 0, rho_params.Dij_lbl);
	dpd_buf4_init(&Aints, CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
	dpd_dot24(&D, &Aints, &I, 0, 0, 1.0, 1.0);
	dpd_dot24(&D, &Aints, &I, 1, 0, 1.0, 1.0);
	dpd_buf4_close(&Aints);
	dpd_file2_close(&D);
  
	dpd_file2_close(&I);

	/* I'ij <-- sum_kl <ik||jl> (D_kl + D_lk) + sum_KL <iK|jL> (D_KL + D_LK) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'ij");

	dpd_file2_init(&D, CC_OEI, 0, 0, 0, rho_params.Dij_lbl);
	dpd_buf4_init(&Aints, CC_AINTS, 0, 0, 0, 0, 0, 1, "A <ij|kl>");
	dpd_dot24(&D, &Aints, &I, 0, 0, 1.0, 1.0);
	dpd_dot24(&D, &Aints, &I, 1, 0, 1.0, 1.0);
	dpd_buf4_close(&Aints);
	dpd_file2_close(&D);

	dpd_file2_init(&D, CC_OEI, 0, 0, 0, rho_params.DIJ_lbl);
	dpd_buf4_init(&Aints, CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
	dpd_dot24(&D, &Aints, &I, 0, 0, 1.0, 1.0);
	dpd_dot24(&D, &Aints, &I, 1, 0, 1.0, 1.0);
	dpd_buf4_close(&Aints);
	dpd_file2_close(&D);
  
	dpd_file2_close(&I);
      }
      else if(params.ref == 2) { /** UHF **/

	/* I'IJ <-- sum_KL <IK||JL> (D_KL + D_LK) + sum_kl <Ik|Jl> (D_kl + D_lk) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'IJ");

	dpd_file2_init(&D, CC_OEI, 0, 0, 0, rho_params.DIJ_lbl);
	dpd_buf4_init(&Aints, CC_AINTS, 0, 0, 0, 0, 0, 1, "A <IJ|KL>");
	dpd_dot24(&D, &Aints, &I, 0, 0, 1.0, 1.0);
	dpd_dot24(&D, &Aints, &I, 1, 0, 1.0, 1.0);
	dpd_buf4_close(&Aints);
	dpd_file2_close(&D);

	dpd_file2_init(&D, CC_OEI, 0, 2, 2, rho_params.Dij_lbl);
	dpd_buf4_init(&Aints, CC_AINTS, 0, 22, 22, 22, 22, 0, "A <Ij|Kl>");
	dpd_dot24(&D, &Aints, &I, 0, 0, 1.0, 1.0);
	dpd_dot24(&D, &Aints, &I, 1, 0, 1.0, 1.0);
	dpd_buf4_close(&Aints);
	dpd_file2_close(&D);
  
	dpd_file2_close(&I);

	/* I'ij <-- sum_kl <ik||jl> (D_kl + D_lk) + sum_KL <iK|jL> (D_KL + D_LK) */
	dpd_file2_init(&I, CC_OEI, 0, 2, 2, "I'ij");

	dpd_file2_init(&D, CC_OEI, 0, 2, 2, rho_params.Dij_lbl);
	dpd_buf4_init(&Aints, CC_AINTS, 0, 10, 10, 10, 10, 1, "A <ij|kl>");
	dpd_dot24(&D, &Aints, &I, 0, 0, 1.0, 1.0);
	dpd_dot24(&D, &Aints, &I, 1, 0, 1.0, 1.0);
	dpd_buf4_close(&Aints);
	dpd_file2_close(&D);

	dpd_file2_init(&D, CC_OEI, 0, 0, 0, rho_params.DIJ_lbl);
	dpd_buf4_init(&Aints, CC_AINTS, 0, 22, 22, 22, 22, 0, "A <Ij|Kl>");
	dpd_dot13(&D, &Aints, &I, 0, 0, 1.0, 1.0);
	dpd_dot13(&D, &Aints, &I, 1, 0, 1.0, 1.0);
	dpd_buf4_close(&Aints);
	dpd_file2_close(&D);
  
	dpd_file2_close(&I);

      }

      if(params.ref == 0) { /** RHF **/
	/* I'IJ <-- sum_KA <IK||JA> (D_KA + D_AK) + sum_ka <Ik|Ja> (D_ka + D_ak) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'IJ");

	dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
	dpd_dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
	dpd_dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_buf4_close(&Eints);

	dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.Dia_lbl);
	dpd_dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.Dai_lbl);
	dpd_dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_buf4_close(&Eints);
  
	dpd_file2_close(&I);
      }
      else if(params.ref == 1) { /** ROHF **/

	/* I'IJ <-- sum_KA <IK||JA> (D_KA + D_AK) + sum_ka <Ik|Ja> (D_ka + D_ak) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'IJ");

	dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
	dpd_dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
	dpd_dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_buf4_close(&Eints);

	dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.Dia_lbl);
	dpd_dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.Dai_lbl);
	dpd_dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_buf4_close(&Eints);
  
	dpd_file2_close(&I);

	/* I'ij <-- sum_ka <ik||ja> (D_ka + D_ak) + sum_KA <iK|jA> (D_KA + D_AK) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'ij");

	dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.Dia_lbl);
	dpd_dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.Dai_lbl);
	dpd_dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_buf4_close(&Eints);

	dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
	dpd_dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
	dpd_dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_buf4_close(&Eints);
  
	dpd_file2_close(&I);
      }
      else if(params.ref == 2) { /** UHF **/

	/* I'IJ <-- sum_KA <IK||JA> (D_KA + D_AK) + sum_ka <Ik|Ja> (D_ka + D_ak) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'IJ");

	dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 20, 2, 20, 0, "E <IJ||KA> (I>J,KA)");
	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
	dpd_dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
	dpd_dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_buf4_close(&Eints);

	dpd_buf4_init(&Eints, CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
	dpd_file2_init(&D, CC_OEI, 0, 2, 3, rho_params.Dia_lbl);
	dpd_dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_file2_init(&D, CC_OEI, 0, 2, 3, rho_params.Dai_lbl);
	dpd_dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_buf4_close(&Eints);
  
	dpd_file2_close(&I);

	/* I'ij <-- sum_ka <ik||ja> (D_ka + D_ak) + sum_KA <iK|jA> (D_KA + D_AK) */
	dpd_file2_init(&I, CC_OEI, 0, 2, 2, "I'ij");

	dpd_buf4_init(&Eints, CC_EINTS, 0, 10, 30, 12, 30, 0, "E <ij||ka> (i>j,ka)");
	dpd_file2_init(&D, CC_OEI, 0, 2, 3, rho_params.Dia_lbl);
	dpd_dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_file2_init(&D, CC_OEI, 0, 2, 3, rho_params.Dai_lbl);
	dpd_dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_buf4_close(&Eints);

	dpd_buf4_init(&Eints, CC_EINTS, 0, 23, 27, 23, 27, 0, "E <iJ|kA>");
	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
	dpd_dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
	dpd_dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_buf4_close(&Eints);
  
	dpd_file2_close(&I);
      }

      if(params.ref == 0) { /** RHF **/
	/* I'IJ <-- sum_AK <JK||IA> (D_AK + D_KA) + sum_ak <Jk|Ia> (D_ak + D_ka) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'IJ");

	dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
	dpd_dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
	dpd_dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_buf4_close(&Eints);

	dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.Dia_lbl);
	dpd_dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.Dai_lbl);
	dpd_dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_buf4_close(&Eints);
  
	dpd_file2_close(&I);
      }
      else if(params.ref == 1) { /** ROHF **/

	/* I'IJ <-- sum_AK <JK||IA> (D_AK + D_KA) + sum_ak <Jk|Ia> (D_ak + D_ka) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'IJ");

	dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
	dpd_dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
	dpd_dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_buf4_close(&Eints);

	dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.Dia_lbl);
	dpd_dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.Dai_lbl);
	dpd_dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_buf4_close(&Eints);
  
	dpd_file2_close(&I);

	/* I'ij <-- sum_ak <jk||ia> (D_ak + D_ka) + sum_AK <jK|iA> (D_AK + D_KA) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'ij");

	dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.Dia_lbl);
	dpd_dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.Dai_lbl);
	dpd_dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_buf4_close(&Eints);

	dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
	dpd_dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
	dpd_dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_buf4_close(&Eints);
  
	dpd_file2_close(&I);
      }
      else if(params.ref == 2) { /** UHF **/

	/* I'IJ <-- sum_AK <JK||IA> (D_AK + D_KA) + sum_ak <Jk|Ia> (D_ak + D_ka) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'IJ");

	dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 20, 2, 20, 0, "E <IJ||KA> (I>J,KA)");
	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
	dpd_dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
	dpd_dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_buf4_close(&Eints);

	dpd_buf4_init(&Eints, CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
	dpd_file2_init(&D, CC_OEI, 0, 2, 3, rho_params.Dia_lbl);
	dpd_dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_file2_init(&D, CC_OEI, 0, 2, 3, rho_params.Dai_lbl);
	dpd_dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_buf4_close(&Eints);
  
	dpd_file2_close(&I);

	/* I'ij <-- sum_ak <jk||ia> (D_ak + D_ka) + sum_AK <jK|iA> (D_AK + D_KA) */
	dpd_file2_init(&I, CC_OEI, 0, 2, 2, "I'ij");

	dpd_buf4_init(&Eints, CC_EINTS, 0, 10, 30, 12, 30, 0, "E <ij||ka> (i>j,ka)");
	dpd_file2_init(&D, CC_OEI, 0, 2, 3, rho_params.Dia_lbl);
	dpd_dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_file2_init(&D, CC_OEI, 0, 2, 3, rho_params.Dai_lbl);
	dpd_dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_buf4_close(&Eints);

	dpd_buf4_init(&Eints, CC_EINTS, 0, 23, 27, 23, 27, 0, "E <iJ|kA>");
	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
	dpd_dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
	dpd_dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_buf4_close(&Eints);
  
	dpd_file2_close(&I);

      }

      if(params.ref == 0) { /** RHF **/
	/* I'IJ <-- sum_AB <IA||JB> (D_AB + D_BA) + sum_ab <Ia|Jb> (D_ab + D_ba) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'IJ");

	dpd_file2_init(&D, CC_OEI, 0, 1, 1, rho_params.DAB_lbl);
	dpd_buf4_init(&Cints, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
	dpd_dot24(&D, &Cints, &I, 0, 0, 1.0, 1.0);
	dpd_dot24(&D, &Cints, &I, 1, 0, 1.0, 1.0);
	dpd_buf4_close(&Cints);
	dpd_file2_close(&D);

	dpd_file2_init(&D, CC_OEI, 0, 1, 1, rho_params.Dab_lbl);
	dpd_buf4_init(&Cints, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
	dpd_dot24(&D, &Cints, &I, 0, 0, 1.0, 1.0);
	dpd_dot24(&D, &Cints, &I, 1, 0, 1.0, 1.0);
	dpd_buf4_close(&Cints);
	dpd_file2_close(&D);

	dpd_file2_close(&I);
      }
      else if(params.ref == 1) { /** ROHF **/

	/* I'IJ <-- sum_AB <IA||JB> (D_AB + D_BA) + sum_ab <Ia|Jb> (D_ab + D_ba) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'IJ");

	dpd_file2_init(&D, CC_OEI, 0, 1, 1, rho_params.DAB_lbl);
	dpd_buf4_init(&Cints, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
	dpd_dot24(&D, &Cints, &I, 0, 0, 1.0, 1.0);
	dpd_dot24(&D, &Cints, &I, 1, 0, 1.0, 1.0);
	dpd_buf4_close(&Cints);
	dpd_file2_close(&D);

	dpd_file2_init(&D, CC_OEI, 0, 1, 1, rho_params.Dab_lbl);
	dpd_buf4_init(&Cints, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
	dpd_dot24(&D, &Cints, &I, 0, 0, 1.0, 1.0);
	dpd_dot24(&D, &Cints, &I, 1, 0, 1.0, 1.0);
	dpd_buf4_close(&Cints);
	dpd_file2_close(&D);

	dpd_file2_close(&I);

	/* I'ij <-- sum_ab <ia||jb> (D_ab + D_ba) + sum_AB <iA|jB> (D_AB + D_BA) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'ij");

	dpd_file2_init(&D, CC_OEI, 0, 1, 1, rho_params.Dab_lbl);
	dpd_buf4_init(&Cints, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
	dpd_dot24(&D, &Cints, &I, 0, 0, 1.0, 1.0);
	dpd_dot24(&D, &Cints, &I, 1, 0, 1.0, 1.0);
	dpd_buf4_close(&Cints);
	dpd_file2_close(&D);

	dpd_file2_init(&D, CC_OEI, 0, 1, 1, rho_params.DAB_lbl);
	dpd_buf4_init(&Cints, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
	dpd_dot24(&D, &Cints, &I, 0, 0, 1.0, 1.0);
	dpd_dot24(&D, &Cints, &I, 1, 0, 1.0, 1.0);
	dpd_buf4_close(&Cints);
	dpd_file2_close(&D);

	dpd_file2_close(&I);
      }
      else if(params.ref == 2) { /** UHF **/

	/* I'IJ <-- sum_AB <IA||JB> (D_AB + D_BA) + sum_ab <Ia|Jb> (D_ab + D_ba) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'IJ");

	dpd_file2_init(&D, CC_OEI, 0, 1, 1, rho_params.DAB_lbl);
	dpd_buf4_init(&Cints, CC_CINTS, 0, 20, 20, 20, 20, 0, "C <IA||JB>");
	dpd_dot24(&D, &Cints, &I, 0, 0, 1.0, 1.0);
	dpd_dot24(&D, &Cints, &I, 1, 0, 1.0, 1.0);
	dpd_buf4_close(&Cints);
	dpd_file2_close(&D);

	dpd_file2_init(&D, CC_OEI, 0, 3, 3, rho_params.Dab_lbl);
	dpd_buf4_init(&Cints, CC_CINTS, 0, 24, 24, 24, 24, 0, "C <Ia|Jb>");
	dpd_dot24(&D, &Cints, &I, 0, 0, 1.0, 1.0);
	dpd_dot24(&D, &Cints, &I, 1, 0, 1.0, 1.0);
	dpd_buf4_close(&Cints);
	dpd_file2_close(&D);

	dpd_file2_close(&I);

	/* I'ij <-- sum_ab <ia||jb> (D_ab + D_ba) + sum_AB <iA|jB> (D_AB + D_BA) */
	dpd_file2_init(&I, CC_OEI, 0, 2, 2, "I'ij");

	dpd_file2_init(&D, CC_OEI, 0, 3, 3, rho_params.Dab_lbl);
	dpd_buf4_init(&Cints, CC_CINTS, 0, 30, 30, 30, 30, 0, "C <ia||jb>");
	dpd_dot24(&D, &Cints, &I, 0, 0, 1.0, 1.0);
	dpd_dot24(&D, &Cints, &I, 1, 0, 1.0, 1.0);
	dpd_buf4_close(&Cints);
	dpd_file2_close(&D);

	dpd_file2_init(&D, CC_OEI, 0, 1, 1, rho_params.DAB_lbl);
	dpd_buf4_init(&Cints, CC_CINTS, 0, 27, 27, 27, 27, 0, "C <iA|jB>");
	dpd_dot24(&D, &Cints, &I, 0, 0, 1.0, 1.0);
	dpd_dot24(&D, &Cints, &I, 1, 0, 1.0, 1.0);
	dpd_buf4_close(&Cints);
	dpd_file2_close(&D);

	dpd_file2_close(&I);

      }

      if(params.ref == 0) { /** RHF **/
	/* I'IJ <-- 2 [ 2 (Gjklm - Gjkml) <ik|lm>] */
	dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'IJ");

	dpd_buf4_init(&G, CC_GAMMA, 0, 0, 0, 0, 0, 0, "GIjKl");
	dpd_buf4_scmcopy(&G, CC_GAMMA, "2 Gijkl - Gijlk", 2);
	dpd_buf4_sort_axpy(&G, CC_GAMMA, pqsr, 0, 0, "2 Gijkl - Gijlk", -1);
	dpd_buf4_close(&G);
	dpd_buf4_init(&G, CC_GAMMA, 0, 0, 0, 0, 0, 0, "2 Gijkl - Gijlk");
	dpd_buf4_init(&Aints, CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
	dpd_contract442(&Aints, &G, &I, 0, 0, 2.0, 1.0); 
	dpd_buf4_close(&Aints);
	dpd_buf4_close(&G);
	dpd_file2_close(&I);
      }
      else if(params.ref == 1) { /** ROHF **/

	/* I'IJ <-- sum_KLM <IK||LM> G(JK,LM) + 2 sum_kLm <Ik|Lm> G(Jk,Lm) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'IJ");

	dpd_buf4_init(&Aints, CC_AINTS, 0, 0, 2, 0, 0, 1, "A <ij|kl>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 0, 2, 2, 2, 0, "GIJKL");
	dpd_contract442(&Aints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Aints);

	dpd_buf4_init(&Aints, CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 0, 0, 0, 0, 0, "GIjKl");
	dpd_contract442(&Aints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Aints);

	dpd_file2_close(&I);

	/* I'ij <-- sum_klm <ik||lm> G(jk,lm) + 2 sum_KlM <Ki|Ml> G(Kj,Ml) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'ij");

	dpd_buf4_init(&Aints, CC_AINTS, 0, 0, 2, 0, 0, 1, "A <ij|kl>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 0, 2, 2, 2, 0, "Gijkl");
	dpd_contract442(&Aints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Aints);

	dpd_buf4_init(&Aints, CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 0, 0, 0, 0, 0, "GIjKl");
	dpd_contract442(&Aints, &G, &I, 1, 1, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Aints);

	dpd_file2_close(&I);

      }
      else if(params.ref == 2) { /** UHF **/

	/* I'IJ <-- sum_KLM <IK||LM> G(JK,LM) + 2 sum_kLm <Ik|Lm> G(Jk,Lm) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'IJ");

	dpd_buf4_init(&Aints, CC_AINTS, 0, 0, 2, 0, 0, 1, "A <IJ|KL>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 0, 2, 2, 2, 0, "GIJKL");
	dpd_contract442(&Aints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Aints);

	dpd_buf4_init(&Aints, CC_AINTS, 0, 22, 22, 22, 22, 0, "A <Ij|Kl>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 22, 22, 22, 22, 0, "GIjKl");
	dpd_contract442(&Aints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Aints);

	dpd_file2_close(&I);

	/* I'ij <-- sum_klm <ik||lm> G(jk,lm) + 2 sum_KlM <Ki|Ml> G(Kj,Ml) */
	dpd_file2_init(&I, CC_OEI, 0, 2, 2, "I'ij");

	dpd_buf4_init(&Aints, CC_AINTS, 0, 10, 12, 10, 10, 1, "A <ij|kl>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 10, 12, 12, 12, 0, "Gijkl");
	dpd_contract442(&Aints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Aints);

	dpd_buf4_init(&Aints, CC_AINTS, 0, 22, 22, 22, 22, 0, "A <Ij|Kl>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 22, 22, 22, 22, 0, "GIjKl");
	dpd_contract442(&Aints, &G, &I, 1, 1, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Aints);

	dpd_file2_close(&I);

      }

      if(params.ref == 0) { /** RHF **/
	/* I'IJ <-- sum_ABC <IA||BC> G(JA,BC) + 2 sum_AbC <aI|bC> G(aJ,bC) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'IJ");

	dpd_buf4_init(&G, CC_GAMMA, 0, 11, 5, 11, 5, 0, "GCiAb");
	dpd_buf4_scmcopy(&G, CC_GAMMA, "2 Gciab - Gciba", 2);
	dpd_buf4_sort_axpy(&G, CC_GAMMA, pqsr, 11, 5, "2 Gciab - Gciba", -1);
	dpd_buf4_close(&G);

	dpd_buf4_init(&G, CC_GAMMA, 0, 11, 5, 11, 5, 0, "2 Gciab - Gciba");
	dpd_buf4_sort(&G, CC_GAMMA, qpsr, 10, 5, "2 Giabc - Giacb");
	dpd_buf4_close(&G);

	dpd_buf4_init(&G, CC_GAMMA, 0, 10, 5, 10, 5, 0, "2 Giabc - Giacb");
	dpd_buf4_init(&Fints, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
	dpd_contract442(&Fints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&Fints);
	dpd_buf4_close(&G);

	dpd_file2_close(&I);
      }
      else if(params.ref == 1) { /** ROHF **/

	/* I'IJ <-- sum_ABC <IA||BC> G(JA,BC) + 2 sum_AbC <aI|bC> G(aJ,bC) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'IJ");

	dpd_buf4_init(&Fints, CC_FINTS, 0, 10, 7, 10, 5, 1, "F <ia|bc>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 11, 7, 11, 7, 0, "GCIAB");
	dpd_buf4_sort(&G, CC_TMP0, qprs, 10, 7, "GICAB");
	dpd_buf4_close(&G);
	dpd_buf4_init(&G, CC_TMP0, 0, 10, 7, 10, 7, 0, "GICAB");
	dpd_contract442(&Fints, &G, &I, 0, 0, -2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Fints);

	dpd_buf4_init(&Fints, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 11, 5, 11, 5, 0, "GcIaB");
	dpd_buf4_sort(&G, CC_TMP0, qpsr, 10, 5, "GIcBa");
	dpd_buf4_close(&G);
	dpd_buf4_init(&G, CC_TMP0, 0, 10, 5, 10, 5, 0, "GIcBa");
	dpd_contract442(&Fints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Fints);

	dpd_file2_close(&I);

	/* I'ij <-- sum_abc <ia||bcC> G(ja,bc) + 2 sum_AbC <Ai|Bc> G(Aj,Bc) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'ij");

	dpd_buf4_init(&Fints, CC_FINTS, 0, 10, 7, 10, 5, 1, "F <ia|bc>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 11, 7, 11, 7, 0, "Gciab");
	dpd_buf4_sort(&G, CC_TMP0, qprs, 10, 7, "Gicab");
	dpd_buf4_close(&G);
	dpd_buf4_init(&G, CC_TMP0, 0, 10, 7, 10, 7, 0, "Gicab");
	dpd_contract442(&Fints, &G, &I, 0, 0, -2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Fints);

	dpd_buf4_init(&Fints, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 11, 5, 11, 5, 0, "GCiAb");
	dpd_buf4_sort(&G, CC_TMP0, qpsr, 10, 5, "GiCbA");
	dpd_buf4_close(&G);
	dpd_buf4_init(&G, CC_TMP0, 0, 10, 5, 10, 5, 0, "GiCbA");
	dpd_contract442(&Fints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Fints);

	dpd_file2_close(&I);

      }
      else if(params.ref == 2) { /** UHF **/

	/* I'IJ <-- sum_ABC <IA||BC> G(JA,BC) + 2 sum_AbC <aI|bC> G(aJ,bC) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'IJ");

	dpd_buf4_init(&Fints, CC_FINTS, 0, 20, 7, 20, 5, 1, "F <IA|BC>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 21, 7, 21, 7, 0, "GCIAB");
	dpd_buf4_sort(&G, CC_TMP0, qprs, 20, 7, "GICAB");
	dpd_buf4_close(&G);
	dpd_buf4_init(&G, CC_TMP0, 0, 20, 7, 20, 7, 0, "GICAB");
	dpd_contract442(&Fints, &G, &I, 0, 0, -2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Fints);

	dpd_buf4_init(&Fints, CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 25, 29, 25, 29, 0, "GcIaB");
	dpd_buf4_sort(&G, CC_TMP0, qpsr, 24, 28, "GIcBa");
	dpd_buf4_close(&G);
	dpd_buf4_init(&G, CC_TMP0, 0, 24, 28, 24, 28, 0, "GIcBa");
	dpd_contract442(&Fints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Fints);

	dpd_file2_close(&I);

	/* I'ij <-- sum_abc <ia||bc> G(ja,bc) + 2 sum_AbC <Ai|Bc> G(Aj,Bc) */
	dpd_file2_init(&I, CC_OEI, 0, 2, 2, "I'ij");

	dpd_buf4_init(&Fints, CC_FINTS, 0, 30, 17, 30, 15, 1, "F <ia|bc>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 31, 17, 31, 17, 0, "Gciab");
	dpd_buf4_sort(&G, CC_TMP0, qprs, 30, 17, "Gicab");
	dpd_buf4_close(&G);
	dpd_buf4_init(&G, CC_TMP0, 0, 30, 17, 30, 17, 0, "Gicab");
	dpd_contract442(&Fints, &G, &I, 0, 0, -2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Fints);

	dpd_buf4_init(&Fints, CC_FINTS, 0, 27, 29, 27, 29, 0, "F <iA|bC>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 26, 28, 26, 28, 0, "GCiAb");
	dpd_buf4_sort(&G, CC_TMP0, qpsr, 27, 29, "GiCbA");
	dpd_buf4_close(&G);
	dpd_buf4_init(&G, CC_TMP0, 0, 27, 29, 27, 29, 0, "GiCbA");
	dpd_contract442(&Fints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Fints);

	dpd_file2_close(&I);

      }

      if(params.ref == 0) { /** RHF **/
	/* I'IJ <-- sum_KAB <IK||AB> G(JK,AB) + 2 sum_kAb <Ik|Ab> G(Jk,Ab) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'IJ");

	dpd_buf4_init(&G, CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
	dpd_buf4_scmcopy(&G, CC_GAMMA, "2 Gijab - Gijba", 2);
	dpd_buf4_sort_axpy(&G, CC_GAMMA, pqsr, 0, 5, "2 Gijab - Gijba", -1);
	dpd_buf4_close(&G);

	dpd_buf4_init(&G, CC_GAMMA, 0, 0, 5, 0, 5, 0, "2 Gijab - Gijba");
	dpd_buf4_init(&Dints, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
	dpd_contract442(&Dints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&Dints);
	dpd_buf4_close(&G);

	dpd_file2_close(&I);
      }
      else if(params.ref == 1) { /** ROHF **/

	/* I'IJ <-- sum_KAB <IK||AB> G(JK,AB) + 2 sum_kAb <Ik|Ab> G(Jk,Ab) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'IJ");

	dpd_buf4_init(&Dints, CC_DINTS, 0, 0, 7, 0, 7, 0, "D <ij||ab> (ij,a>b)");
	dpd_buf4_init(&G, CC_GAMMA, 0, 0, 7, 2, 7, 0, "GIJAB");
	dpd_contract442(&Dints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Dints);

	dpd_buf4_init(&Dints, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
	dpd_contract442(&Dints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Dints);

	dpd_file2_close(&I);

	/* I'ij <-- sum_kab <ik||ab> G(jk,ab) + 2 sum_KaB <iK|aB> G(jK,aB) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'ij");

	dpd_buf4_init(&Dints, CC_DINTS, 0, 0, 7, 0, 7, 0, "D <ij||ab> (ij,a>b)");
	dpd_buf4_init(&G, CC_GAMMA, 0, 0, 7, 2, 7, 0, "Gijab");
	dpd_contract442(&Dints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Dints);

	dpd_buf4_init(&Dints, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
	dpd_contract442(&Dints, &G, &I, 1, 1, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Dints);

	dpd_file2_close(&I);

      }
      else if(params.ref == 2) { /** UHF **/

	/* I'IJ <-- sum_KAB <IK||AB> G(JK,AB) + 2 sum_kAb <Ik|Ab> G(Jk,Ab) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'IJ");

	dpd_buf4_init(&Dints, CC_DINTS, 0, 0, 7, 0, 7, 0, "D <IJ||AB> (IJ,A>B)");
	dpd_buf4_init(&G, CC_GAMMA, 0, 0, 7, 2, 7, 0, "GIJAB");
	dpd_contract442(&Dints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Dints);

	dpd_buf4_init(&Dints, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 22, 28, 22, 28, 0, "GIjAb");
	dpd_contract442(&Dints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Dints);

	dpd_file2_close(&I);

	/* I'ij <-- sum_kab <ik||ab> G(jk,ab) + 2 sum_KaB <iK|aB> G(jK,aB) */
	dpd_file2_init(&I, CC_OEI, 0, 2, 2, "I'ij");

	dpd_buf4_init(&Dints, CC_DINTS, 0, 10, 17, 10, 17, 0, "D <ij||ab> (ij,a>b)");
	dpd_buf4_init(&G, CC_GAMMA, 0, 10, 17, 12, 17, 0, "Gijab");
	dpd_contract442(&Dints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Dints);

	dpd_buf4_init(&Dints, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 22, 28, 22, 28, 0, "GIjAb");
	dpd_contract442(&Dints, &G, &I, 1, 1, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Dints);

	dpd_file2_close(&I);

      }

      if(params.ref == 0) { /** RHF **/
	/* I'IJ <-- 2 sum_AKB <IA||KB> G(JA,KB) + 2 sum_aKb <Ia|Kb> G(Ja,Kb) -
	   2 sum_akB <Ik|Ba> GJakB(Jk,Ba) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'IJ");

	dpd_buf4_init(&Cints, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIBJA");
	dpd_contract442(&Cints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Cints);

	dpd_buf4_init(&Cints, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbJa");
	dpd_contract442(&Cints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Cints);

	dpd_buf4_init(&Dints, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbjA");
	dpd_buf4_sort(&G, CC_TMP0, prsq, 0, 5, "GIbjA (Ij,Ab)");
	dpd_buf4_close(&G);
	dpd_buf4_init(&G, CC_TMP0, 0, 0, 5, 0, 5, 0, "GIbjA (Ij,Ab)");
	dpd_contract442(&Dints, &G, &I, 0, 0, -2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Dints);

	dpd_file2_close(&I);
      }
      else if(params.ref == 1) { /** ROHF **/

	/* I'IJ <-- 2 sum_AKB <IA||KB> G(JA,KB) + 2 sum_aKb <Ia|Kb> G(Ja,Kb) -
	   2 sum_akB <Ik|Ba> GJakB(Jk,Ba) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'IJ");

	dpd_buf4_init(&Cints, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIBJA");
	dpd_contract442(&Cints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Cints);

	dpd_buf4_init(&Cints, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbJa");
	dpd_contract442(&Cints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Cints);

	dpd_buf4_init(&Dints, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbjA");
	dpd_buf4_sort(&G, CC_TMP0, prsq, 0, 5, "GIbjA (Ij,Ab)");
	dpd_buf4_close(&G);
	dpd_buf4_init(&G, CC_TMP0, 0, 0, 5, 0, 5, 0, "GIbjA (Ij,Ab)");
	dpd_contract442(&Dints, &G, &I, 0, 0, -2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Dints);

	dpd_file2_close(&I);

	/* I'ij <-- 2 sum_akb <ia||kb> G(ja,kb) + 2 sum_AkB <iA|kB> G(jA,kB) +
	   2 sum_AKb <iK|bA> GjAKb(jK,bA) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'ij");

	dpd_buf4_init(&Cints, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "Gibja");
	dpd_contract442(&Cints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Cints);

	dpd_buf4_init(&Cints, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GiBjA");
	dpd_contract442(&Cints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Cints);

	dpd_buf4_init(&Dints, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GiBJa");
	dpd_buf4_sort(&G, CC_TMP0, prsq, 0, 5, "GiBJa (iJ,aB)");
	dpd_buf4_close(&G);
	dpd_buf4_init(&G, CC_TMP0, 0, 0, 5, 0, 5, 0, "GiBJa (iJ,aB)");
	dpd_contract442(&Dints, &G, &I, 0, 0, -2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Dints);

	dpd_file2_close(&I);

      }
      else if(params.ref == 2) { /** UHF **/

	/* I'IJ <-- 2 sum_AKB <IA||KB> G(JA,KB) + 2 sum_aKb <Ia|Kb> G(Ja,Kb) -
	   2 sum_akB <Ik|Ba> GJakB(Jk,Ba) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'IJ");

	dpd_buf4_init(&Cints, CC_CINTS, 0, 20, 20, 20, 20, 0, "C <IA||JB>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 20, 20, 20, 20, 0, "GIBJA");
	dpd_contract442(&Cints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Cints);

	dpd_buf4_init(&Cints, CC_CINTS, 0, 24, 24, 24, 24, 0, "C <Ia|Jb>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 24, 24, 24, 24, 0, "GIbJa");
	dpd_contract442(&Cints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Cints);

	dpd_buf4_init(&Dints, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 24, 27, 24, 27, 0, "GIbjA");
	dpd_buf4_sort(&G, CC_TMP0, prsq, 22, 28, "GIbjA (Ij,Ab)");
	dpd_buf4_close(&G);
	dpd_buf4_init(&G, CC_TMP0, 0, 22, 28, 22, 28, 0, "GIbjA (Ij,Ab)");
	dpd_contract442(&Dints, &G, &I, 0, 0, -2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Dints);

	dpd_file2_close(&I);

	/* I'ij <-- 2 sum_akb <ia||kb> G(ja,kb) + 2 sum_AkB <iA|kB> G(jA,kB) +
	   2 sum_AKb <iK|bA> GjAKb(jK,bA) */
	dpd_file2_init(&I, CC_OEI, 0, 2, 2, "I'ij");

	dpd_buf4_init(&Cints, CC_CINTS, 0, 30, 30, 30, 30, 0, "C <ia||jb>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 30, 30, 30, 30, 0, "Gibja");
	dpd_contract442(&Cints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Cints);

	dpd_buf4_init(&Cints, CC_CINTS, 0, 27, 27, 27, 27, 0, "C <iA|jB>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 27, 27, 27, 27, 0, "GiBjA");
	dpd_contract442(&Cints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Cints);

	dpd_buf4_init(&Dints, CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 27, 24, 27, 24, 0, "GiBJa");
	dpd_buf4_sort(&G, CC_TMP0, prsq, 23, 29, "GiBJa (iJ,aB)");
	dpd_buf4_close(&G);
	dpd_buf4_init(&G, CC_TMP0, 0, 23, 29, 23, 29, 0, "GiBJa (iJ,aB)");
	dpd_contract442(&Dints, &G, &I, 0, 0, -2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Dints);

	dpd_file2_close(&I);

      }

      if(params.ref == 0) { /** RHF **/
	/* I'IJ <-- 2 sum_KLA <IK||LA> G(JK,LA) + 2 sum_kLa <Ik|La> G(Jk,La)
	   + 2 sum_kAl <kI|lA> G(kJ,lA) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'IJ");

	dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 0, 10, 0, "GIjKa");
	dpd_buf4_scmcopy(&G, CC_GAMMA, "2 Gijka - Gjika", 2);
	dpd_buf4_sort_axpy(&G, CC_GAMMA, qprs, 0, 10, "2 Gijka - Gjika", -1);
	dpd_buf4_close(&G);

	dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 0, 10, 0, "2 Gijka - Gjika");
	dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
	dpd_contract442(&Eints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_contract442(&Eints, &G, &I, 1, 1, 2.0, 1.0);
	dpd_buf4_close(&Eints);
	dpd_buf4_close(&G);

	dpd_file2_close(&I);
      }
      else if(params.ref == 1) { /** ROHF **/

	/* I'IJ <-- 2 sum_KLA <IK||LA> G(JK,LA) + 2 sum_kLa <Ik|La> G(Jk,La)
	   + 2 sum_kAl <kI|lA> G(kJ,lA) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'IJ");

	dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
	dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 2, 10, 0, "GIJKA");
	dpd_contract442(&Eints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Eints);

	dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 0, 10, 0, "GIjKa");
	dpd_contract442(&Eints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Eints);

	dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 0, 10, 0, "GiJkA");
	dpd_contract442(&Eints, &G, &I, 1, 1, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Eints);

	dpd_file2_close(&I);

	/* I'ij <-- 2 sum_kla <ik||la> G(jk,la) + 2 sum_KlA <iK|lA> G(jK,lA)
	   + 2 sum_KaL <Ki|La> G(Kj,La) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'ij");

	dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
	dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 2, 10, 0, "Gijka");
	dpd_contract442(&Eints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Eints);

	dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 0, 10, 0, "GiJkA");
	dpd_contract442(&Eints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Eints);

	dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 0, 10, 0, "GIjKa");
	dpd_contract442(&Eints, &G, &I, 1, 1, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Eints);

	dpd_file2_close(&I);
      }
      else if(params.ref == 2) { /** UHF **/

	/* I'IJ <-- 2 sum_KLA <IK||LA> G(JK,LA) + 2 sum_kLa <Ik|La> G(Jk,La)
	   + 2 sum_kAl <kI|lA> G(kJ,lA) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'IJ");

	dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 20, 2, 20, 0, "E <IJ||KA> (I>J,KA)");
	dpd_buf4_init(&G, CC_GAMMA, 0, 0, 20, 2, 20, 0, "GIJKA");
	dpd_contract442(&Eints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Eints);

	dpd_buf4_init(&Eints, CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 22, 24, 22, 24, 0, "GIjKa");
	dpd_contract442(&Eints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Eints);

	dpd_buf4_init(&Eints, CC_EINTS, 0, 23, 27, 23, 27, 0, "E <iJ|kA>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 23, 27, 23, 27, 0, "GiJkA");
	dpd_contract442(&Eints, &G, &I, 1, 1, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Eints);

	dpd_file2_close(&I);

	/* I'ij <-- 2 sum_kla <ik||la> G(jk,la) + 2 sum_KlA <iK|lA> G(jK,lA)
	   + 2 sum_KaL <Ki|La> G(Kj,La) */
	dpd_file2_init(&I, CC_OEI, 0, 2, 2, "I'ij");

	dpd_buf4_init(&Eints, CC_EINTS, 0, 10, 30, 12, 30, 0, "E <ij||ka> (i>j,ka)");
	dpd_buf4_init(&G, CC_GAMMA, 0, 10, 30, 12, 30, 0, "Gijka");
	dpd_contract442(&Eints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Eints);

	dpd_buf4_init(&Eints, CC_EINTS, 0, 23, 27, 23, 27, 0, "E <iJ|kA>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 23, 27, 23, 27, 0, "GiJkA");
	dpd_contract442(&Eints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Eints);

	dpd_buf4_init(&Eints, CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 22, 24, 22, 24, 0, "GIjKa");
	dpd_contract442(&Eints, &G, &I, 1, 1, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Eints);

	dpd_file2_close(&I);

      }

      if(params.ref == 0) { /** RHF **/
	/* I'IJ <-- 2 sum_AKL <K>L||IA> G(K>L,JA) + 2 sum_aKl <Kl|Ia> G(Kl,Ja) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'IJ");

	dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 0, 10, 0, "GIjKa");
	dpd_buf4_scmcopy(&G, CC_GAMMA, "2 Gijka - Gjika", 2);
	dpd_buf4_sort_axpy(&G, CC_GAMMA, qprs, 0, 10, "2 Gijka - Gjika", -1);
	dpd_buf4_close(&G);

	dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 0, 10, 0, "2 Gijka - Gjika");
	dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
	dpd_contract442(&Eints, &G, &I, 2, 2, 2.0, 1.0);
	dpd_buf4_close(&Eints);
	dpd_buf4_close(&G);

	dpd_file2_close(&I);
      }
      else if(params.ref == 1) { /** ROHF **/

	/* I'IJ <-- 2 sum_AKL <K>L||IA> G(K>L,JA) + 2 sum_aKl <Kl|Ia> G(Kl,Ja) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'IJ");

	dpd_buf4_init(&Eints, CC_EINTS, 0, 2, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
	dpd_buf4_init(&G, CC_GAMMA, 0, 2, 10, 2, 10, 0, "GIJKA");
	dpd_contract442(&Eints, &G, &I, 2, 2, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Eints);

	dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 0, 10, 0, "GIjKa");
	dpd_contract442(&Eints, &G, &I, 2, 2, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Eints);

	dpd_file2_close(&I);

	/* I'ij <-- 2 sum_akl <k>l||ia> G(k>l,ja) + 2 sum_AkL <kL|iA> G(kL,jA) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'ij");

	dpd_buf4_init(&Eints, CC_EINTS, 0, 2, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
	dpd_buf4_init(&G, CC_GAMMA, 0, 2, 10, 2, 10, 0, "Gijka");
	dpd_contract442(&Eints, &G, &I, 2, 2, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Eints);

	dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 0, 10, 0, "GiJkA");
	dpd_contract442(&Eints, &G, &I, 2, 2, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Eints);

	dpd_file2_close(&I);

      }
      else if(params.ref == 2) { /** UHF **/

	/* I'IJ <-- 2 sum_AKL <K>L||IA> G(K>L,JA) + 2 sum_aKl <Kl|Ia> G(Kl,Ja) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'IJ");

	dpd_buf4_init(&Eints, CC_EINTS, 0, 2, 20, 2, 20, 0, "E <IJ||KA> (I>J,KA)");
	dpd_buf4_init(&G, CC_GAMMA, 0, 2, 20, 2, 20, 0, "GIJKA");
	dpd_contract442(&Eints, &G, &I, 2, 2, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Eints);

	dpd_buf4_init(&Eints, CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 22, 24, 22, 24, 0, "GIjKa");
	dpd_contract442(&Eints, &G, &I, 2, 2, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Eints);

	dpd_file2_close(&I);

	/* I'ij <-- 2 sum_akl <k>l||ia> G(k>l,ja) + 2 sum_AkL <kL|iA> G(kL,jA) */
	dpd_file2_init(&I, CC_OEI, 0, 2, 2, "I'ij");

	dpd_buf4_init(&Eints, CC_EINTS, 0, 12, 30, 12, 30, 0, "E <ij||ka> (i>j,ka)");
	dpd_buf4_init(&G, CC_GAMMA, 0, 12, 30, 12, 30, 0, "Gijka");
	dpd_contract442(&Eints, &G, &I, 2, 2, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Eints);

	dpd_buf4_init(&Eints, CC_EINTS, 0, 23, 27, 23, 27, 0, "E <iJ|kA>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 23, 27, 23, 27, 0, "GiJkA");
	dpd_contract442(&Eints, &G, &I, 2, 2, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Eints);

	dpd_file2_close(&I);

      }
    }

  }} // namespace psi::ccdensity
