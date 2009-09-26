/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <strings.h>
#include <string.h>
#include <libdpd/dpd.h>
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
  dpdfile2 D, T1, L1, Z;
  dpdbuf4 T2, L2;
  double trace=0.0, dot_AI, dot_IA, dot_ai, dot_ia;
  double factor=0.0;

  if(params.ref == 0 || params.ref == 1) { /** RHF/ROHF **/
    dpd_file2_init(&D, CC_OEI, 0, 0, 0, rho_params.DIJ_lbl);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 7, 2, 7, 0, "LIJAB");
    dpd_contract442(&T2, &L2, &D, 0, 0, -1.0, 0.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2); 
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 5, 0, 5, 0, "LIjAb");
    dpd_contract442(&T2, &L2, &D, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2); 
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
    dpd_contract222(&T1, &L1, &D, 0, 0, -1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_file2_close(&T1);
    trace += dpd_file2_trace(&D);
    dpd_file2_close(&D);

    dpd_file2_init(&D, CC_OEI, 0, 0, 0, rho_params.Dij_lbl);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tijab");
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 7, 2, 7, 0, "Lijab");
    dpd_contract442(&T2, &L2, &D, 0, 0, -1.0, 0.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 5, 0, 5, 0, "LiJaB");
    dpd_contract442(&T2, &L2, &D, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "Lia");
    dpd_contract222(&T1, &L1, &D, 0, 0, -1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_file2_close(&T1);
    trace += dpd_file2_trace(&D); 
    dpd_file2_close(&D);

    dpd_file2_init(&D, CC_OEI, 0, 1, 1, rho_params.DAB_lbl);
    dpd_buf4_init(&L2, CC_GLG, 0, 2, 5, 2, 7, 0, "LIJAB");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    dpd_contract442(&L2, &T2, &D, 3, 3, 1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 5, 0, 5, 0, "LiJaB");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
    dpd_contract442(&L2, &T2, &D, 3, 3, 1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract222(&L1, &T1, &D, 1, 1, 1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_file2_close(&T1);
    trace += dpd_file2_trace(&D); 
    dpd_file2_close(&D);

    dpd_file2_init(&D, CC_OEI, 0, 1, 1, rho_params.Dab_lbl);
    dpd_buf4_init(&L2, CC_GLG, 0, 2, 5, 2, 7, 0, "Lijab");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tijab");
    dpd_contract442(&L2, &T2, &D, 3, 3, 1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 5, 0, 5, 0, "LIjAb");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_contract442(&L2, &T2, &D, 3, 3, 1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "Lia");
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
    dpd_contract222(&L1, &T1, &D, 1, 1, 1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_file2_close(&T1);
    trace += dpd_file2_trace(&D); 
    dpd_file2_close(&D);

    /*fprintf(outfile, "\n\tTrace of onepdm = %20.15f\n", trace);*/

    /* This term is * L0 = 0 for excited states */
    if (rho_params.L_ground) {
      dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
      dpd_file2_copy(&T1, CC_OEI, rho_params.DIA_lbl);
      dpd_file2_close(&T1);
    }
    dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.DIA_lbl);

    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
    dpd_dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_buf4_close(&T2);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "Lia");
    dpd_dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_buf4_close(&T2);
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&Z, CC_TMP0, 0, 0, 0, "Z(I,M)");
    dpd_contract222(&T1, &L1, &Z, 0, 0, 1.0, 0.0);
    dpd_file2_close(&L1);
    dpd_contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    dpd_file2_close(&T1);
    dpd_file2_close(&Z);
    /*  D(I,A) << L2(MN,EF) T2(IN,EF) T(M,A) + L2(Mn,Ef) T2(In,Ef) T(M,A) */
    dpd_file2_init(&Z, CC_TMP0, 0, 0, 0, "Z(I,M)");
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 7, 2, 7, 0, "LIJAB");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
    dpd_contract442(&T2, &L2, &Z, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 5, 0, 5, 0, "LIjAb");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_contract442(&T2, &L2, &Z, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    dpd_file2_close(&Z);
    dpd_file2_close(&T1);
    /* T2(MN,AF) L2(MN,EF) T(I,E) + T2(Mn,Af) L2(Mn,Ef) T(I,E) */
    dpd_file2_init(&Z, CC_TMP0, 0, 1, 1, "Z(A,E)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    dpd_buf4_init(&L2, CC_GLG, 0, 2, 5, 2, 7, 0, "LIJAB");
    dpd_contract442(&T2, &L2, &Z, 2, 2, 1.0, 0.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 5, 0, 5, 0, "LIjAb");
    dpd_contract442(&T2, &L2, &Z, 2, 2, 1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract222(&T1, &Z, &D, 0, 0, -1.0, 1.0);
    dpd_file2_close(&T1);
    dpd_file2_close(&Z);
    dpd_file2_close(&D);

    /* This term is * L0 = 0 for excited states */
    if (rho_params.L_ground) {
      dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
      dpd_file2_copy(&T1, CC_OEI, rho_params.Dia_lbl);
      dpd_file2_close(&T1);
    }
    dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.Dia_lbl);

    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tijab");
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "Lia");
    dpd_dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_buf4_close(&T2);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
    dpd_dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_buf4_close(&T2);
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "Lia");
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
    dpd_file2_init(&Z, CC_TMP0, 0, 0, 0, "Z(i,m)");
    dpd_contract222(&T1, &L1, &Z, 0, 0, 1.0, 0.0);
    dpd_file2_close(&L1);
    dpd_contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    dpd_file2_close(&T1);
    dpd_file2_close(&Z);
    dpd_file2_init(&Z, CC_TMP0, 0, 0, 0, "Z(i,m)");
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 7, 2, 7, 0, "Lijab");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tijab");
    dpd_contract442(&T2, &L2, &Z, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 5, 0, 5, 0, "LiJaB");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
    dpd_contract442(&T2, &L2, &Z, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
    dpd_contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    dpd_file2_close(&Z);
    dpd_file2_close(&T1);
    dpd_file2_init(&Z, CC_TMP0, 0, 1, 1, "Z(a,e)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tijab");
    dpd_buf4_init(&L2, CC_GLG, 0, 2, 5, 2, 7, 0, "Lijab");
    dpd_contract442(&T2, &L2, &Z, 2, 2, 1.0, 0.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 5, 0, 5, 0, "LiJaB");
    dpd_contract442(&T2, &L2, &Z, 2, 2, 1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
    dpd_contract222(&T1, &Z, &D, 0, 0, -1.0, 1.0);
    dpd_file2_close(&T1);
    dpd_file2_close(&Z);
    dpd_file2_close(&D);

    /* Note that these blocks are still stored occ/vir */
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
    dpd_file2_copy(&L1, CC_OEI, rho_params.DAI_lbl);
    dpd_file2_close(&L1);

    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "Lia");
    dpd_file2_copy(&L1, CC_OEI, rho_params.Dai_lbl);
    dpd_file2_close(&L1);

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
    fprintf(outfile,"\tOverlaps of onepdm after ground-state parts added.\n");
    fprintf(outfile,"\t<DIA|DIA> = %15.10lf     <Dia|Dia> = %15.10lf\n", dot_IA, dot_ia);
    fprintf(outfile,"\t<DAI|DAI> = %15.10lf     <Dai|Dai> = %15.10lf\n", dot_AI, dot_ai);
    fprintf(outfile,"\t<Dpq|Dqp> = %15.10lf\n", dot_IA+dot_ia+dot_AI+dot_ai);
    */
  }
  else if(params.ref == 2) { /** UHF **/

    if(params.wfn == "CCSD_T" && params.dertype == 1) {
      /* For CCSD(T) gradients, some density contributions are
	 calculated in cctriples */
      factor = 1.0;
      dpd_file2_init(&D, CC_OEI, 0, 0, 0, "DIJ");
      dpd_file2_copy(&D, CC_OEI, rho_params.DIJ_lbl);
      dpd_file2_close(&D);
      dpd_file2_init(&D, CC_OEI, 0, 2, 2, "Dij");
      dpd_file2_copy(&D, CC_OEI, rho_params.Dij_lbl);
      dpd_file2_close(&D);
      dpd_file2_init(&D, CC_OEI, 0, 1, 1, "DAB");
      dpd_file2_copy(&D, CC_OEI, rho_params.DAB_lbl);
      dpd_file2_close(&D);
      dpd_file2_init(&D, CC_OEI, 0, 3, 3, "Dab");
      dpd_file2_copy(&D, CC_OEI, rho_params.Dab_lbl);
      dpd_file2_close(&D);
    }
    else factor = 0.0;

    dpd_file2_init(&D, CC_OEI, 0, 0, 0, rho_params.DIJ_lbl);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 7, 2, 7, 0, "LIJAB");
    dpd_contract442(&T2, &L2, &D, 0, 0, -1.0, factor);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2); 
    dpd_buf4_init(&T2, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_buf4_init(&L2, CC_GLG, 0, 22, 28, 22, 28, 0, "LIjAb");
    dpd_contract442(&T2, &L2, &D, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2); 
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
    dpd_contract222(&T1, &L1, &D, 0, 0, -1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_file2_close(&T1);
    trace += dpd_file2_trace(&D);
    dpd_file2_close(&D);

    dpd_file2_init(&D, CC_OEI, 0, 2, 2, rho_params.Dij_lbl);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 17, 12, 17, 0, "tijab");
    dpd_buf4_init(&L2, CC_GLG, 0, 10, 17, 12, 17, 0, "Lijab");
    dpd_contract442(&T2, &L2, &D, 0, 0, -1.0, factor);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    dpd_buf4_init(&L2, CC_GLG, 0, 23, 29, 23, 29, 0, "LiJaB");
    dpd_contract442(&T2, &L2, &D, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
    dpd_file2_init(&L1, CC_GLG, 0, 2, 3, "Lia");
    dpd_contract222(&T1, &L1, &D, 0, 0, -1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_file2_close(&T1);
    trace += dpd_file2_trace(&D); 
    dpd_file2_close(&D);

    dpd_file2_init(&D, CC_OEI, 0, 1, 1, rho_params.DAB_lbl);
    dpd_buf4_init(&L2, CC_GLG, 0, 2, 5, 2, 7, 0, "LIJAB");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    dpd_contract442(&L2, &T2, &D, 3, 3, 1.0, factor);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GLG, 0, 23, 29, 23, 29, 0, "LiJaB");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    dpd_contract442(&L2, &T2, &D, 3, 3, 1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract222(&L1, &T1, &D, 1, 1, 1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_file2_close(&T1);
    trace += dpd_file2_trace(&D); 
    dpd_file2_close(&D);

    dpd_file2_init(&D, CC_OEI, 0, 3, 3, rho_params.Dab_lbl);
    dpd_buf4_init(&L2, CC_GLG, 0, 12, 15, 12, 17, 0, "Lijab");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 12, 15, 12, 17, 0, "tijab");
    dpd_contract442(&L2, &T2, &D, 3, 3, 1.0, factor);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GLG, 0, 22, 28, 22, 28, 0, "LIjAb");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_contract442(&L2, &T2, &D, 3, 3, 1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_file2_init(&L1, CC_GLG, 0, 2, 3, "Lia");
    dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
    dpd_contract222(&L1, &T1, &D, 1, 1, 1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_file2_close(&T1);
    trace += dpd_file2_trace(&D); 
    dpd_file2_close(&D);

    /*fprintf(outfile, "\n\tTrace of onepdm = %20.15f\n", trace);*/

    /* This term is * L0 = 0 for excited states */
    if (rho_params.L_ground) {
      dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
      dpd_file2_copy(&T1, CC_OEI, rho_params.DIA_lbl);
      dpd_file2_close(&T1);
    }
    dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.DIA_lbl);

    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
    dpd_dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_buf4_close(&T2);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_file2_init(&L1, CC_GLG, 0, 2, 3, "Lia");
    dpd_dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_buf4_close(&T2);
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&Z, CC_TMP0, 0, 0, 0, "Z(I,M)");
    dpd_contract222(&T1, &L1, &Z, 0, 0, 1.0, 0.0);
    dpd_file2_close(&L1);
    dpd_contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    dpd_file2_close(&T1);
    dpd_file2_close(&Z);
    /*  D(I,A) << L2(MN,EF) T2(IN,EF) T(M,A) + L2(Mn,Ef) T2(In,Ef) T(M,A) */
    dpd_file2_init(&Z, CC_TMP0, 0, 0, 0, "Z(I,M)");
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 7, 2, 7, 0, "LIJAB");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
    dpd_contract442(&T2, &L2, &Z, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GLG, 0, 22, 28, 22, 28, 0, "LIjAb");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_contract442(&T2, &L2, &Z, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    dpd_file2_close(&Z);
    dpd_file2_close(&T1);
    /* T2(MN,AF) L2(MN,EF) T(I,E) + T2(Mn,Af) L2(Mn,Ef) T(I,E) */
    dpd_file2_init(&Z, CC_TMP0, 0, 1, 1, "Z(A,E)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    dpd_buf4_init(&L2, CC_GLG, 0, 2, 5, 2, 7, 0, "LIJAB");
    dpd_contract442(&T2, &L2, &Z, 2, 2, 1.0, 0.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_buf4_init(&L2, CC_GLG, 0, 22, 28, 22, 28, 0, "LIjAb");
    dpd_contract442(&T2, &L2, &Z, 2, 2, 1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract222(&T1, &Z, &D, 0, 0, -1.0, 1.0);
    dpd_file2_close(&T1);
    dpd_file2_close(&Z);
    dpd_file2_close(&D);

    /* This term is * L0 = 0 for excited states */
    if (rho_params.L_ground) {
      dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
      dpd_file2_copy(&T1, CC_OEI, rho_params.Dia_lbl);
      dpd_file2_close(&T1);
    }
    dpd_file2_init(&D, CC_OEI, 0, 2, 3, rho_params.Dia_lbl);

    dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 15, 12, 17, 0, "tijab");
    dpd_file2_init(&L1, CC_GLG, 0, 2, 3, "Lia");
    dpd_dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_buf4_close(&T2);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
    dpd_dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_buf4_close(&T2);
    dpd_file2_init(&L1, CC_GLG, 0, 2, 3, "Lia");
    dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
    dpd_file2_init(&Z, CC_TMP0, 0, 2, 2, "Z(i,m)");
    dpd_contract222(&T1, &L1, &Z, 0, 0, 1.0, 0.0);
    dpd_file2_close(&L1);
    dpd_contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    dpd_file2_close(&T1);
    dpd_file2_close(&Z);
    dpd_file2_init(&Z, CC_TMP0, 0, 2, 2, "Z(i,m)");
    dpd_buf4_init(&L2, CC_GLG, 0, 10, 17, 12, 17, 0, "Lijab");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 17, 12, 17, 0, "tijab");
    dpd_contract442(&T2, &L2, &Z, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GLG, 0, 23, 29, 23, 29, 0, "LiJaB");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    dpd_contract442(&T2, &L2, &Z, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
    dpd_contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    dpd_file2_close(&Z);
    dpd_file2_close(&T1);
    dpd_file2_init(&Z, CC_TMP0, 0, 3, 3, "Z(a,e)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 12, 15, 12, 17, 0, "tijab");
    dpd_buf4_init(&L2, CC_GLG, 0, 12, 15, 12, 17, 0, "Lijab");
    dpd_contract442(&T2, &L2, &Z, 2, 2, 1.0, 0.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    dpd_buf4_init(&L2, CC_GLG, 0, 23, 29, 23, 29, 0, "LiJaB");
    dpd_contract442(&T2, &L2, &Z, 2, 2, 1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
    dpd_contract222(&T1, &Z, &D, 0, 0, -1.0, 1.0);
    dpd_file2_close(&T1);
    dpd_file2_close(&Z);
    dpd_file2_close(&D);

    /* Note that these blocks are still stored occ/vir */
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
    dpd_file2_copy(&L1, CC_OEI, rho_params.DAI_lbl);
    dpd_file2_close(&L1);

    dpd_file2_init(&L1, CC_GLG, 0, 2, 3, "Lia");
    dpd_file2_copy(&L1, CC_OEI, rho_params.Dai_lbl);
    dpd_file2_close(&L1);

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
    fprintf(outfile,"\tOverlaps of onepdm after ground-state parts added.\n");
    fprintf(outfile,"\t<DIA|DIA> = %15.10lf     <Dia|Dia> = %15.10lf\n", dot_IA, dot_ia);
    fprintf(outfile,"\t<DAI|DAI> = %15.10lf     <Dai|Dai> = %15.10lf\n", dot_AI, dot_ai);
    fprintf(outfile,"\t<Dpq|Dqp> = %15.10lf\n", dot_IA+dot_ia+dot_AI+dot_ai);
    */
  }
}

}} // namespace psi::ccdensity
