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

    /* Iia(): Build the occupied-virtual block of the orbital Lagrangian
    ** using the expression given in lag.c.
    **
    ** Note that the code currently produces only the unique I_IJ terms,
    ** but the actual contractions still need to be spin-adapted for
    ** greater efficiency.
    **
    ** TDC, 2/2008
    ** */

    void Iia(struct RHO_Params rho_params)
    {
      dpdfile2 F, D, I;
      dpdbuf4 G, Aints, Fints, Eints, Dints, Cints;

      if(params.ref == 0) { /** RHF **/
	/* I'IA <-- sum_J fIJ (DAJ + DJA) + sum_B fIB (DAB + DBA) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I'IA");

	dpd_file2_init(&F, CC_OEI, 0, 0, 0, "fIJ");
	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
	dpd_contract222(&F, &D, &I, 0, 1, 1.0, 0.0);
	dpd_file2_close(&D);
	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
	dpd_contract222(&F, &D, &I, 0, 1, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_file2_close(&F);

	dpd_file2_init(&F, CC_OEI, 0, 0, 1, "fIA");
	dpd_file2_init(&D, CC_OEI, 0, 1, 1, rho_params.DAB_lbl);
	dpd_contract222(&F, &D, &I, 0, 0, 1.0, 1.0);
	dpd_contract222(&F, &D, &I, 0, 1, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_file2_close(&F);

	dpd_file2_close(&I);
      }
      else if(params.ref == 1) { /** ROHF **/

	/* I'IA <-- sum_J fIJ (DAJ + DJA) + sum_B fIB (DAB + DBA) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I'IA");

	dpd_file2_init(&F, CC_OEI, 0, 0, 0, "fIJ");
	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
	dpd_contract222(&F, &D, &I, 0, 1, 1.0, 0.0);
	dpd_file2_close(&D);
	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
	dpd_contract222(&F, &D, &I, 0, 1, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_file2_close(&F);

	dpd_file2_init(&F, CC_OEI, 0, 0, 1, "fIA");
	dpd_file2_init(&D, CC_OEI, 0, 1, 1, rho_params.DAB_lbl);
	dpd_contract222(&F, &D, &I, 0, 0, 1.0, 1.0);
	dpd_contract222(&F, &D, &I, 0, 1, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_file2_close(&F);

	dpd_file2_close(&I);

	/* I'ia <-- sum_j fij (Daj + Dja) + sum_b fib (Dab + Dba) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I'ia");

	dpd_file2_init(&F, CC_OEI, 0, 0, 0, "fij");
	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.Dai_lbl);
	dpd_contract222(&F, &D, &I, 0, 1, 1.0, 0.0);
	dpd_file2_close(&D);
	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.Dia_lbl);
	dpd_contract222(&F, &D, &I, 0, 1, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_file2_close(&F);

	dpd_file2_init(&F, CC_OEI, 0, 0, 1, "fia");
	dpd_file2_init(&D, CC_OEI, 0, 1, 1, rho_params.Dab_lbl);
	dpd_contract222(&F, &D, &I, 0, 0, 1.0, 1.0);
	dpd_contract222(&F, &D, &I, 0, 1, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_file2_close(&F);

	dpd_file2_close(&I);
      }
      else if(params.ref == 2) { /** UHF **/

	/* I'IA <-- sum_J fIJ (DAJ + DJA) + sum_B fIB (DAB + DBA) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I'IA");

	dpd_file2_init(&F, CC_OEI, 0, 0, 0, "fIJ");
	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
	dpd_contract222(&F, &D, &I, 0, 1, 1.0, 0.0);
	dpd_file2_close(&D);
	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
	dpd_contract222(&F, &D, &I, 0, 1, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_file2_close(&F);

	dpd_file2_init(&F, CC_OEI, 0, 0, 1, "fIA");
	dpd_file2_init(&D, CC_OEI, 0, 1, 1, rho_params.DAB_lbl);
	dpd_contract222(&F, &D, &I, 0, 0, 1.0, 1.0);
	dpd_contract222(&F, &D, &I, 0, 1, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_file2_close(&F);

	dpd_file2_close(&I);

	/* I'ia <-- sum_j fij (Daj + Dja) + sum_b fib (Dab + Dba) */
	dpd_file2_init(&I, CC_OEI, 0, 2, 3, "I'ia");

	dpd_file2_init(&F, CC_OEI, 0, 2, 2, "fij");
	dpd_file2_init(&D, CC_OEI, 0, 2, 3, rho_params.Dai_lbl);
	dpd_contract222(&F, &D, &I, 0, 1, 1.0, 0.0);
	dpd_file2_close(&D);
	dpd_file2_init(&D, CC_OEI, 0, 2, 3, rho_params.Dia_lbl);
	dpd_contract222(&F, &D, &I, 0, 1, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_file2_close(&F);

	dpd_file2_init(&F, CC_OEI, 0, 2, 3, "fia");
	dpd_file2_init(&D, CC_OEI, 0, 3, 3, rho_params.Dab_lbl);
	dpd_contract222(&F, &D, &I, 0, 0, 1.0, 1.0);
	dpd_contract222(&F, &D, &I, 0, 1, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_file2_close(&F);

	dpd_file2_close(&I);

      }

  

      if(params.ref == 0) { /** RHF **/
	/* I'IA <-- sum_JKL <LK||JI> G(LK,JA) + 2 sum_jKl <lK|jI> G(lK,jA) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I'IA");

	dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 0, 10, 0, "GIjKa");
	dpd_buf4_scmcopy(&G, CC_GAMMA, "2 Gijka - Gjika", 2);
	dpd_buf4_sort_axpy(&G, CC_GAMMA, qprs, 0, 10, "2 Gijka - Gjika", -1);
	dpd_buf4_close(&G);

	dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 0, 10, 0, "2 Gijka - Gjika");
	dpd_buf4_init(&Aints, CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
	dpd_contract442(&Aints, &G, &I, 3, 3, 2.0, 1.0);
	dpd_buf4_close(&Aints);
	dpd_buf4_close(&G);

	dpd_file2_close(&I);
      }
      else if(params.ref == 1) { /** ROHF **/

	/* I'IA <-- sum_JKL <LK||JI> G(LK,JA) + 2 sum_jKl <lK|jI> G(lK,jA) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I'IA");

	dpd_buf4_init(&Aints, CC_AINTS, 0, 2, 0, 0, 0, 1, "A <ij|kl>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 2, 10, 2, 10, 0, "GIJKA");
	dpd_contract442(&Aints, &G, &I, 3, 3, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Aints);

	dpd_buf4_init(&Aints, CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 0, 10, 0, "GiJkA");
	dpd_contract442(&Aints, &G, &I, 3, 3, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Aints);

	dpd_file2_close(&I);

	/* I'ia <-- sum_jkl <lk||ji> G(lk,ja) + 2 sum_JkL <Lk|Ji> G(Lk,Ja) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I'ia");

	dpd_buf4_init(&Aints, CC_AINTS, 0, 2, 0, 0, 0, 1, "A <ij|kl>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 2, 10, 2, 10, 0, "Gijka");
	dpd_contract442(&Aints, &G, &I, 3, 3, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Aints);

	dpd_buf4_init(&Aints, CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 0, 10, 0, "GIjKa");
	dpd_contract442(&Aints, &G, &I, 3, 3, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Aints);

	dpd_file2_close(&I);
      }
      else if(params.ref == 2) { /** UHF **/

	/* I'IA <-- sum_JKL <LK||JI> G(LK,JA) + 2 sum_jKl <lK|jI> G(lK,jA) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I'IA");

	dpd_buf4_init(&Aints, CC_AINTS, 0, 2, 0, 0, 0, 1, "A <IJ|KL>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 2, 20, 2, 20, 0, "GIJKA");
	dpd_contract442(&Aints, &G, &I, 3, 3, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Aints);

	dpd_buf4_init(&Aints, CC_AINTS, 0, 22, 22, 22, 22, 0, "A <Ij|Kl>");
	dpd_buf4_sort(&Aints, CC_AINTS, qpsr, 23, 23, "A <iJ|kL>");
	dpd_buf4_close(&Aints);
	dpd_buf4_init(&Aints, CC_AINTS, 0, 23, 23, 23, 23, 0, "A <iJ|kL>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 23, 27, 23, 27, 0, "GiJkA");
	dpd_contract442(&Aints, &G, &I, 3, 3, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Aints);

	dpd_file2_close(&I);

	/* I'ia <-- sum_jkl <lk||ji> G(lk,ja) + 2 sum_JkL <Lk|Ji> G(Lk,Ja) */
	dpd_file2_init(&I, CC_OEI, 0, 2, 3, "I'ia");

	dpd_buf4_init(&Aints, CC_AINTS, 0, 12, 10, 10, 10, 1, "A <ij|kl>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 12, 30, 12, 30, 0, "Gijka");
	dpd_contract442(&Aints, &G, &I, 3, 3, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Aints);

	dpd_buf4_init(&Aints, CC_AINTS, 0, 22, 22, 22, 22, 0, "A <Ij|Kl>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 22, 24, 22, 24, 0, "GIjKa");
	dpd_contract442(&Aints, &G, &I, 3, 3, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Aints);

	dpd_file2_close(&I);

      }

  

      if(params.ref == 0) { /** RHF **/
	/* I'IA <-- sum_BCD <IB||CD> G(AB,CD) + 2 sum_bCd <Ib|Cd> G(Ab,Cd) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I'IA");

	dpd_buf4_init(&G, CC_GAMMA, 0, 5, 5, 5, 5, 0, "GAbCd");
	dpd_buf4_scmcopy(&G, CC_GAMMA, "2 Gabcd - Gabdc", 2);
	dpd_buf4_sort_axpy(&G, CC_GAMMA, pqsr, 5, 5, "2 Gabcd - Gabdc", -1);
	dpd_buf4_close(&G);

	dpd_buf4_init(&G, CC_GAMMA, 0, 5, 5, 5, 5, 0, "2 Gabcd - Gabdc");
	dpd_buf4_init(&Fints, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
	dpd_contract442(&Fints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&Fints);
	dpd_buf4_close(&G);


	dpd_file2_close(&I);
      }
      else if(params.ref == 1) { /** ROHF **/

	/* I'IA <-- sum_BCD <IB||CD> G(AB,CD) + 2 sum_bCd <Ib|Cd> G(Ab,Cd) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I'IA");

	dpd_buf4_init(&Fints, CC_FINTS, 0, 10, 7, 10, 5, 1, "F <ia|bc>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 5, 7, 7, 7, 0, "GABCD");
	dpd_contract442(&Fints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Fints);

	dpd_buf4_init(&Fints, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 5, 5, 5, 5, 0, "GAbCd");
	dpd_contract442(&Fints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Fints);

	dpd_file2_close(&I);

	/* I'ia <-- sum_bcd <ib||cd> G(ab,cd) + 2 sum_BcD <Dc|Bi> G(Dc,Ba) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I'ia");

	dpd_buf4_init(&Fints, CC_FINTS, 0, 10, 7, 10, 5, 1, "F <ia|bc>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 5, 7, 7, 7, 0, "Gabcd");
	dpd_contract442(&Fints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Fints);

	dpd_buf4_init(&Fints, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
	dpd_buf4_sort(&Fints, CC_TMP0, srqp, 5, 11, "F <cb|ai>");
	dpd_buf4_close(&Fints);
	dpd_buf4_init(&Fints, CC_TMP0, 0, 5, 11, 5, 11, 0, "F <cb|ai>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 5, 5, 5, 5, 0, "GAbCd");
	dpd_contract442(&Fints, &G, &I, 3, 3, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Fints);

	dpd_file2_close(&I);
      }
      else if(params.ref == 2) { /** UHF **/

	/* I'IA <-- sum_BCD <IB||CD> G(AB,CD) + 2 sum_bCd <Ib|Cd> G(Ab,Cd) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I'IA");

	dpd_buf4_init(&Fints, CC_FINTS, 0, 20, 7, 20, 5, 1, "F <IA|BC>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 5, 7, 7, 7, 0, "GABCD");
	dpd_contract442(&Fints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Fints);

	dpd_buf4_init(&Fints, CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 28, 28, 28, 28, 0, "GAbCd");
	dpd_contract442(&Fints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Fints);

	dpd_file2_close(&I);

	/* I'ia <-- sum_bcd <ib||cd> G(ab,cd) + 2 sum_BcD <Dc|Bi> G(Dc,Ba) */
	dpd_file2_init(&I, CC_OEI, 0, 2, 3, "I'ia");

	dpd_buf4_init(&Fints, CC_FINTS, 0, 30, 17, 30, 15, 1, "F <ia|bc>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 15, 17, 17, 17, 0, "Gabcd");
	dpd_contract442(&Fints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Fints);

	dpd_buf4_init(&Fints, CC_FINTS, 0, 27, 29, 27, 29, 0, "F <iA|bC>");
	dpd_buf4_sort(&Fints, CC_TMP0, srqp, 28, 26, "F <Cb|Ai>");
	dpd_buf4_close(&Fints);
	dpd_buf4_init(&Fints, CC_TMP0, 0, 28, 26, 28, 26, 0, "F <Cb|Ai>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 28, 28, 28, 28, 0, "GAbCd");
	dpd_contract442(&Fints, &G, &I, 3, 3, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Fints);

	dpd_file2_close(&I);

      }

  

      if(params.ref == 0) { /** RHF **/
	/* I'IA <-- 2 sum_JKB <JI||KB> G(JA,KB) + 2 sum_jKb <Ij|Kb> G(Aj,Kb) +
	   2 sum_jkB <jI|kB> G(jA,kB) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I'IA");

	dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIBJA");
	dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
	dpd_contract442(&Eints, &G, &I, 1, 1, 2.0, 1.0);
	dpd_buf4_close(&Eints);
	dpd_buf4_close(&G);

	dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GiBJa");
	dpd_buf4_sort(&G, CC_TMP0, qprs, 11, 10, "GiBJa (Bi,Ja)");
	dpd_buf4_close(&G);
	dpd_buf4_init(&G, CC_TMP0, 0, 11, 10, 11, 10, 0, "GiBJa (Bi,Ja)");
	dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
	dpd_contract442(&Eints, &G, &I, 0, 0, -2.0, 1.0);
	dpd_buf4_close(&Eints);
	dpd_buf4_close(&G);

	dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GiBjA");
	dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
	dpd_contract442(&Eints, &G, &I, 1, 1, 2.0, 1.0);
	dpd_buf4_close(&Eints);
	dpd_buf4_close(&G);

	dpd_file2_close(&I);
      }
      else if(params.ref == 1) { /** ROHF **/

	/* I'IA <-- 2 sum_JKB <JI||KB> G(JA,KB) + 2 sum_jKb <Ij|Kb> G(Aj,Kb) +
	   2 sum_jkB <jI|kB> G(jA,kB) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I'IA");

	dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIBJA");
	dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
	dpd_contract442(&Eints, &G, &I, 1, 1, 2.0, 1.0);
	dpd_buf4_close(&Eints);
	dpd_buf4_close(&G);

	dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GiBJa");
	dpd_buf4_sort(&G, CC_TMP0, qprs, 11, 10, "GiBJa (Bi,Ja)");
	dpd_buf4_close(&G);
	dpd_buf4_init(&G, CC_TMP0, 0, 11, 10, 11, 10, 0, "GiBJa (Bi,Ja)");
	dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
	dpd_contract442(&Eints, &G, &I, 0, 0, -2.0, 1.0);
	dpd_buf4_close(&Eints);
	dpd_buf4_close(&G);

	dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GiBjA");
	dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
	dpd_contract442(&Eints, &G, &I, 1, 1, 2.0, 1.0);
	dpd_buf4_close(&Eints);
	dpd_buf4_close(&G);

	dpd_file2_close(&I);

	/* I'ia <-- 2 sum_jkb <ji||kb> G(ja,kb) + 2 sum_JkB <iJ|kB> G(aJ,kB) +
	   2 sum_JKb <Ji|kB> G(Ja,Kb) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I'ia");

	dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "Gibja");
	dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
	dpd_contract442(&Eints, &G, &I, 1, 1, 2.0, 1.0);
	dpd_buf4_close(&Eints);
	dpd_buf4_close(&G);

	dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbjA");
	dpd_buf4_sort(&G, CC_TMP0, qprs, 11, 10, "GIbjA (bI,jA)");
	dpd_buf4_close(&G);
	dpd_buf4_init(&G, CC_TMP0, 0, 11, 10, 11, 10, 0, "GIbjA (bI,jA)");
	dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
	dpd_contract442(&Eints, &G, &I, 0, 0, -2.0, 1.0);
	dpd_buf4_close(&Eints);
	dpd_buf4_close(&G);

	dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbJa");
	dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
	dpd_contract442(&Eints, &G, &I, 1, 1, 2.0, 1.0);
	dpd_buf4_close(&Eints);
	dpd_buf4_close(&G);

	dpd_file2_close(&I);
      }
      else if(params.ref == 2) { /** UHF **/

	/* I'IA <-- 2 sum_JKB <JI||KB> G(JA,KB) + 2 sum_jKb <Ij|Kb> G(Aj,Kb) +
	   2 sum_jkB <jI|kB> G(jA,kB) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I'IA");

	dpd_buf4_init(&G, CC_GAMMA, 0, 20, 20, 20, 20, 0, "GIBJA");
	dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 20, 2, 20, 0, "E <IJ||KA> (I>J,KA)");
	dpd_contract442(&Eints, &G, &I, 1, 1, 2.0, 1.0);
	dpd_buf4_close(&Eints);
	dpd_buf4_close(&G);

	dpd_buf4_init(&G, CC_GAMMA, 0, 27, 24, 27, 24, 0, "GiBJa");
	dpd_buf4_sort(&G, CC_TMP0, qprs, 26, 24, "GiBJa (Bi,Ja)");
	dpd_buf4_close(&G);
	dpd_buf4_init(&G, CC_TMP0, 0, 26, 24, 26, 24, 0, "GiBJa (Bi,Ja)");
	dpd_buf4_init(&Eints, CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
	dpd_contract442(&Eints, &G, &I, 0, 0, -2.0, 1.0);
	dpd_buf4_close(&Eints);
	dpd_buf4_close(&G);

	dpd_buf4_init(&G, CC_GAMMA, 0, 27, 27, 27, 27, 0, "GiBjA");
	dpd_buf4_init(&Eints, CC_EINTS, 0, 23, 27, 23, 27, 0, "E <iJ|kA>");
	dpd_contract442(&Eints, &G, &I, 1, 1, 2.0, 1.0);
	dpd_buf4_close(&Eints);
	dpd_buf4_close(&G);

	dpd_file2_close(&I);

	/* I'ia <-- 2 sum_jkb <ji||kb> G(ja,kb) + 2 sum_JkB <iJ|kB> G(aJ,kB) +
	   2 sum_JKb <Ji|Kb> G(Ja,Kb) */
	dpd_file2_init(&I, CC_OEI, 0, 2, 3, "I'ia");

	dpd_buf4_init(&G, CC_GAMMA, 0, 30, 30, 30, 30, 0, "Gibja");
	dpd_buf4_init(&Eints, CC_EINTS, 0, 10, 30, 12, 30, 0, "E <ij||ka> (i>j,ka)");
	dpd_contract442(&Eints, &G, &I, 1, 1, 2.0, 1.0);
	dpd_buf4_close(&Eints);
	dpd_buf4_close(&G);

	dpd_buf4_init(&G, CC_GAMMA, 0, 24, 27, 24, 27, 0, "GIbjA");
	dpd_buf4_sort(&G, CC_TMP0, qprs, 25, 27, "GIbjA (bI,jA)");
	dpd_buf4_close(&G);
	dpd_buf4_init(&G, CC_TMP0, 0, 25, 27, 25, 27, 0, "GIbjA (bI,jA)");
	dpd_buf4_init(&Eints, CC_EINTS, 0, 23, 27, 23, 27, 0, "E <iJ|kA>");
	dpd_contract442(&Eints, &G, &I, 0, 0, -2.0, 1.0);
	dpd_buf4_close(&Eints);
	dpd_buf4_close(&G);

	dpd_buf4_init(&G, CC_GAMMA, 0, 24, 24, 24, 24, 0, "GIbJa");
	dpd_buf4_init(&Eints, CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
	dpd_contract442(&Eints, &G, &I, 1, 1, 2.0, 1.0);
	dpd_buf4_close(&Eints);
	dpd_buf4_close(&G);

	dpd_file2_close(&I);

      }

  

      if(params.ref == 0) { /** RHF **/
	/* I'IA <-- sum_BJK <JK||IB> G(JK,AB) + 2 sum_bJk <Jk|Ib> G(Jk,Ab) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I'IA");

	dpd_buf4_init(&G, CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
	dpd_buf4_scmcopy(&G, CC_GAMMA, "2 Gijab - Gijba", 2);
	dpd_buf4_sort_axpy(&G, CC_GAMMA, pqsr, 0, 5, "2 Gijab - Gijba", -1);
	dpd_buf4_close(&G);

	dpd_buf4_init(&G, CC_GAMMA, 0, 0, 5, 0, 5, 0, "2 Gijab - Gijba");
	dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
	dpd_contract442(&Eints, &G, &I, 2, 2, 2.0, 1.0);
	dpd_buf4_close(&Eints);
	dpd_buf4_close(&G);

	dpd_file2_close(&I);
      }
      else if(params.ref == 1) { /** ROHF **/

	/* I'IA <-- sum_BJK <JK||IB> G(JK,AB) + 2 sum_bJk <Jk|Ib> G(Jk,Ab) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I'IA");

	dpd_buf4_init(&Eints, CC_EINTS, 0, 2, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
	dpd_buf4_init(&G, CC_GAMMA, 0, 2, 5, 2, 7, 0, "GIJAB");
	dpd_contract442(&Eints, &G, &I, 2, 2, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Eints);

	dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
	dpd_contract442(&Eints, &G, &I, 2, 2, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Eints);

	dpd_file2_close(&I);

	/* I'ia <-- sum_bjk <jk||ib> G(jk,ab) + 2 sum_BjK <Kj|Bi> G(Kj,Ba) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I'ia");

	dpd_buf4_init(&Eints, CC_EINTS, 0, 2, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
	dpd_buf4_init(&G, CC_GAMMA, 0, 2, 5, 2, 7, 0, "Gijab");
	dpd_contract442(&Eints, &G, &I, 2, 2, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Eints);

	dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
	dpd_buf4_sort(&G, CC_TMP0, qpsr, 0, 5, "GjIbA");
	dpd_buf4_close(&G);
	dpd_buf4_init(&G, CC_TMP0, 0, 0, 5, 0, 5, 0, "GjIbA");
	dpd_contract442(&Eints, &G, &I, 2, 2, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Eints);

	dpd_file2_close(&I);
      }
      else if(params.ref == 2) { /** UHF **/

	/* I'IA <-- sum_BJK <JK||IB> G(JK,AB) + 2 sum_bJk <Jk|Ib> G(Jk,Ab) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I'IA");

	dpd_buf4_init(&Eints, CC_EINTS, 0, 2, 20, 2, 20, 0, "E <IJ||KA> (I>J,KA)");
	dpd_buf4_init(&G, CC_GAMMA, 0, 2, 5, 2, 7, 0, "GIJAB");
	dpd_contract442(&Eints, &G, &I, 2, 2, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Eints);

	dpd_buf4_init(&Eints, CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 22, 28, 22, 28, 0, "GIjAb");
	dpd_contract442(&Eints, &G, &I, 2, 2, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Eints);

	dpd_file2_close(&I);

	/* I'ia <-- sum_bjk <jk||ib> G(jk,ab) + 2 sum_BjK <Kj|Bi> G(Kj,Ba) */
	dpd_file2_init(&I, CC_OEI, 0, 2, 3, "I'ia");

	dpd_buf4_init(&Eints, CC_EINTS, 0, 12, 30, 12, 30, 0, "E <ij||ka> (i>j,ka)");
	dpd_buf4_init(&G, CC_GAMMA, 0, 12, 15, 12, 17, 0, "Gijab");
	dpd_contract442(&Eints, &G, &I, 2, 2, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Eints);

	dpd_buf4_init(&Eints, CC_EINTS, 0, 22, 26, 22, 26, 0, "E <Ij|Ak>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 22, 28, 22, 28, 0, "GIjAb");
	dpd_contract442(&Eints, &G, &I, 3, 3, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Eints);

	dpd_file2_close(&I);

      }

  

      if(params.ref == 0) { /** RHF **/
	/* I'IA <-- sum_JBC <IJ||BC> G(AJ,BC) + 2 sum_jBc <Ij|Bc> G(Aj,Bc) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I'IA");

	dpd_buf4_init(&G, CC_GAMMA, 0, 11, 5, 11, 5, 0, "GCiAb");
	dpd_buf4_scmcopy(&G, CC_GAMMA, "2 Gciab - Gciba", 2);
	dpd_buf4_sort_axpy(&G, CC_GAMMA, pqsr, 11, 5, "2 Gciab - Gciba", -1);
	dpd_buf4_close(&G);

	dpd_buf4_init(&G, CC_GAMMA, 0, 11, 5, 11, 5, 0, "2 Gciab - Gciba");
	dpd_buf4_init(&Dints, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
	dpd_contract442(&Dints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&Dints);
	dpd_buf4_close(&G);

	dpd_file2_close(&I);
      }
      else if(params.ref == 1) { /** ROHF **/

	/* I'IA <-- sum_JBC <IJ||BC> G(AJ,BC) + 2 sum_jBc <Ij|Bc> G(Aj,Bc) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I'IA");

	dpd_buf4_init(&Dints, CC_DINTS, 0, 0, 7, 0, 7, 0, "D <ij||ab> (ij,a>b)");
	dpd_buf4_init(&G, CC_GAMMA, 0, 11, 7, 11, 7, 0, "GCIAB");
	dpd_contract442(&Dints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Dints);

	dpd_buf4_init(&Dints, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 11, 5, 11, 5, 0, "GCiAb");
	dpd_contract442(&Dints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Dints);

	dpd_file2_close(&I);

	/* I'ia <-- sum_jbc <ij||bc> G(aj,bc) + 2 sum_JbC <iJ|bC> G(aJ,bC) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I'ia");

	dpd_buf4_init(&Dints, CC_DINTS, 0, 0, 7, 0, 7, 0, "D <ij||ab> (ij,a>b)");
	dpd_buf4_init(&G, CC_GAMMA, 0, 11, 7, 11, 7, 0, "Gciab");
	dpd_contract442(&Dints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Dints);

	dpd_buf4_init(&Dints, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 11, 5, 11, 5, 0, "GcIaB");
	dpd_contract442(&Dints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Dints);

	dpd_file2_close(&I);
      }
      else if(params.ref == 2) { /** UHF **/

	/* I'IA <-- sum_JBC <IJ||BC> G(AJ,BC) + 2 sum_jBc <Ij|Bc> G(Aj,Bc) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I'IA");

	dpd_buf4_init(&Dints, CC_DINTS, 0, 0, 7, 0, 7, 0, "D <IJ||AB> (IJ,A>B)");
	dpd_buf4_init(&G, CC_GAMMA, 0, 21, 7, 21, 7, 0, "GCIAB");
	dpd_contract442(&Dints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Dints);

	dpd_buf4_init(&Dints, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 26, 28, 26, 28, 0, "GCiAb");
	dpd_contract442(&Dints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Dints);

	dpd_file2_close(&I);

	/* I'ia <-- sum_jbc <ij||bc> G(aj,bc) + 2 sum_JbC <iJ|bC> G(aJ,bC) */
	dpd_file2_init(&I, CC_OEI, 0, 2, 3, "I'ia");

	dpd_buf4_init(&Dints, CC_DINTS, 0, 10, 17, 10, 17, 0, "D <ij||ab> (ij,a>b)");
	dpd_buf4_init(&G, CC_GAMMA, 0, 31, 17, 31, 17, 0, "Gciab");
	dpd_contract442(&Dints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Dints);

	dpd_buf4_init(&Dints, CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 25, 29, 25, 29, 0, "GcIaB");
	dpd_contract442(&Dints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Dints);

	dpd_file2_close(&I);

      }

  

      if(params.ref == 0) { /** RHF **/
	/* I'IA <-- 2 sum_BJC <JC||IB> G(JC,AB) + 2 sum_bJc <Jc|Ib> G(Jc,Ab) +
	   2 sum_bjC <Cj|Ib> G(Cj,Ab) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I'IA");

	dpd_buf4_init(&G, CC_GAMMA, 0, 11, 5, 11, 5, 0, "GCiAb");
	dpd_buf4_scmcopy(&G, CC_GAMMA, "2 Gciab - Gciba", 2);
	dpd_buf4_sort_axpy(&G, CC_GAMMA, pqsr, 11, 5, "2 Gciab - Gciba", -1);
	dpd_buf4_close(&G);

	dpd_buf4_init(&G, CC_GAMMA, 0, 11, 5, 11, 5, 0, "2 Gciab - Gciba");
	dpd_buf4_init(&Cints, CC_CINTS, 0, 11, 11, 11, 11, 0, "C <ai|bj>");
	dpd_contract442(&Cints, &G, &I, 3, 3, 2.0, 1.0);
	dpd_buf4_close(&Cints);
	dpd_buf4_init(&Dints, CC_DINTS, 0, 11, 10, 11, 10, 0, "D <ij|ab> (aj,ib)");
	dpd_contract442(&Dints, &G, &I, 2, 2, 2.0, 1.0);
	dpd_buf4_close(&Dints);
	dpd_buf4_close(&G);

	dpd_file2_close(&I);
      }
      else if(params.ref == 1) { /** ROHF **/

	/* I'IA <-- 2 sum_BJC <JC||IB> G(JC,AB) + 2 sum_bJc <Jc|Ib> G(Jc,Ab) +
	   2 sum_bjC <Cj|Ib> G(Cj,Ab) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I'IA");

	dpd_buf4_init(&G, CC_GAMMA, 0, 11, 5, 11, 7, 0, "GCIAB");
	dpd_buf4_sort(&G, CC_TMP0, qpsr, 10, 5, "GICBA");
	dpd_buf4_close(&G);
	dpd_buf4_init(&G, CC_TMP0, 0, 10, 5, 10, 5, 0, "GICBA");
	dpd_buf4_init(&Cints, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
	dpd_contract442(&Cints, &G, &I, 2, 2, 2.0, 1.0);
	dpd_buf4_close(&Cints);
	dpd_buf4_close(&G);

	dpd_buf4_init(&G, CC_GAMMA, 0, 11, 5, 11, 5, 0, "GcIaB");
	dpd_buf4_sort(&G, CC_TMP0, qpsr, 10, 5, "GIcBa");
	dpd_buf4_close(&G);
	dpd_buf4_init(&G, CC_TMP0, 0, 10, 5, 10, 5, 0, "GIcBa");
	dpd_buf4_init(&Cints, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
	dpd_contract442(&Cints, &G, &I, 2, 2, 2.0, 1.0);
	dpd_buf4_close(&Cints);
	dpd_buf4_close(&G);

	dpd_buf4_init(&Dints, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
	dpd_buf4_sort(&Dints, CC_TMP0, spqr, 11, 10, "D <ij|ab> (bi,ja)");
	dpd_buf4_close(&Dints);
	dpd_buf4_init(&Dints, CC_TMP0, 0, 11, 10, 11, 10, 0, "D <ij|ab> (bi,ja)");
	dpd_buf4_init(&G, CC_GAMMA, 0, 11, 5, 11, 5, 0, "GCiAb");
	dpd_contract442(&Dints, &G, &I, 2, 2, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Dints);

	dpd_file2_close(&I);

	/* I'ia <-- 2 sum_bjc <jc||ib> G(jc,ab) + 2 sum_BjC <jC|iB> G(jC,aB) +
	   2 sum_BJc <cJ|iB> G(cJ,aB) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I'ia");

	dpd_buf4_init(&G, CC_GAMMA, 0, 11, 5, 11, 7, 0, "Gciab");
	dpd_buf4_sort(&G, CC_TMP0, qpsr, 10, 5, "Gicba");
	dpd_buf4_close(&G);
	dpd_buf4_init(&G, CC_TMP0, 0, 10, 5, 10, 5, 0, "Gicba");
	dpd_buf4_init(&Cints, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
	dpd_contract442(&Cints, &G, &I, 2, 2, 2.0, 1.0);
	dpd_buf4_close(&Cints);
	dpd_buf4_close(&G);

	dpd_buf4_init(&G, CC_GAMMA, 0, 11, 5, 11, 5, 0, "GCiAb");
	dpd_buf4_sort(&G, CC_TMP0, qpsr, 10, 5, "GiCbA");
	dpd_buf4_close(&G);
	dpd_buf4_init(&G, CC_TMP0, 0, 10, 5, 10, 5, 0, "GiCbA");
	dpd_buf4_init(&Cints, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
	dpd_contract442(&Cints, &G, &I, 2, 2, 2.0, 1.0);
	dpd_buf4_close(&Cints);
	dpd_buf4_close(&G);

	/* This set of sorted D-integrals is generated in the previous code block */
	dpd_buf4_init(&Dints, CC_TMP0, 0, 11, 10, 11, 10, 0, "D <ij|ab> (bi,ja)");
	dpd_buf4_init(&G, CC_GAMMA, 0, 11, 5, 11, 5, 0, "GcIaB");
	dpd_contract442(&Dints, &G, &I, 2, 2, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Dints);

	dpd_file2_close(&I);
      }
      else if(params.ref == 2) { /** UHF **/

	/* I'IA <-- 2 sum_BJC <JC||IB> G(JC,AB) + 2 sum_bJc <Jc|Ib> G(Jc,Ab) +
	   2 sum_bjC <Cj|Ib> G(Cj,Ab) */
	dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I'IA");

	dpd_buf4_init(&G, CC_GAMMA, 0, 21, 5, 21, 7, 0, "GCIAB");
	dpd_buf4_sort(&G, CC_TMP0, qpsr, 20, 5, "GICBA");
	dpd_buf4_close(&G);
	dpd_buf4_init(&G, CC_TMP0, 0, 20, 5, 20, 5, 0, "GICBA");
	dpd_buf4_init(&Cints, CC_CINTS, 0, 20, 20, 20, 20, 0, "C <IA||JB>");
	dpd_contract442(&Cints, &G, &I, 2, 2, 2.0, 1.0);
	dpd_buf4_close(&Cints);
	dpd_buf4_close(&G);

	dpd_buf4_init(&G, CC_GAMMA, 0, 25, 29, 25, 29, 0, "GcIaB");
	dpd_buf4_sort(&G, CC_TMP0, qpsr, 24, 28, "GIcBa");
	dpd_buf4_close(&G);
	dpd_buf4_init(&G, CC_TMP0, 0, 24, 28, 24, 28, 0, "GIcBa");
	dpd_buf4_init(&Cints, CC_CINTS, 0, 24, 24, 24, 24, 0, "C <Ia|Jb>");
	dpd_contract442(&Cints, &G, &I, 2, 2, 2.0, 1.0);
	dpd_buf4_close(&Cints);
	dpd_buf4_close(&G);

	dpd_buf4_init(&Dints, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
	dpd_buf4_sort(&Dints, CC_TMP0, rqps, 26, 24, "D <Ij|Ab> (Aj,Ib)");
	dpd_buf4_close(&Dints);
	dpd_buf4_init(&Dints, CC_TMP0, 0, 26, 24, 26, 24, 0, "D <Ij|Ab> (Aj,Ib)");
	dpd_buf4_init(&G, CC_GAMMA, 0, 26, 28, 26, 28, 0, "GCiAb");
	dpd_contract442(&Dints, &G, &I, 2, 2, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Dints);

	dpd_file2_close(&I);

	/* I'ia <-- 2 sum_bjc <jc||ib> G(jc,ab) + 2 sum_BjC <jC|iB> G(jC,aB) +
	   2 sum_BJc <cJ|iB> G(cJ,aB) */
	dpd_file2_init(&I, CC_OEI, 0, 2, 3, "I'ia");

	dpd_buf4_init(&G, CC_GAMMA, 0, 31, 15, 31, 17, 0, "Gciab");
	dpd_buf4_sort(&G, CC_TMP0, qpsr, 30, 15, "Gicba");
	dpd_buf4_close(&G);
	dpd_buf4_init(&G, CC_TMP0, 0, 30, 15, 30, 15, 0, "Gicba");
	dpd_buf4_init(&Cints, CC_CINTS, 0, 30, 30, 30, 30, 0, "C <ia||jb>");
	dpd_contract442(&Cints, &G, &I, 2, 2, 2.0, 1.0);
	dpd_buf4_close(&Cints);
	dpd_buf4_close(&G);

	dpd_buf4_init(&G, CC_GAMMA, 0, 26, 28, 26, 28, 0, "GCiAb");
	dpd_buf4_sort(&G, CC_TMP0, qpsr, 27, 29, "GiCbA");
	dpd_buf4_close(&G);
	dpd_buf4_init(&G, CC_TMP0, 0, 27, 29, 27, 29, 0, "GiCbA");
	dpd_buf4_init(&Cints, CC_CINTS, 0, 27, 27, 27, 27, 0, "C <iA|jB>");
	dpd_contract442(&Cints, &G, &I, 2, 2, 2.0, 1.0);
	dpd_buf4_close(&Cints);
	dpd_buf4_close(&G);

	dpd_buf4_init(&Dints, CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");
	dpd_buf4_sort(&Dints, CC_TMP0, rqps, 25, 27, "D <iJ|aB> (aJ,iB)");
	dpd_buf4_close(&Dints);
	dpd_buf4_init(&Dints, CC_TMP0, 0, 25, 27, 25, 27, 0, "D <iJ|aB> (aJ,iB)");
	dpd_buf4_init(&G, CC_GAMMA, 0, 25, 29, 25, 29, 0, "GcIaB");
	dpd_contract442(&Dints, &G, &I, 2, 2, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Dints);

	dpd_file2_close(&I);

      }
  
    }

  }} // namespace psi::ccdensity
