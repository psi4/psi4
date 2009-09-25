/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

void x_onepdm_rohf(struct RHO_Params rho_params);
extern void x_onepdm_uhf(struct RHO_Params rho_params);

void x_onepdm(struct RHO_Params rho_params) {
  if (params.ref == 0 || params.ref == 1)
    x_onepdm_rohf(rho_params);
  else
    x_onepdm_uhf(rho_params);
  return;
}

/* onepdm(): Computes the non-R0 parts of the unrelaxed EOM 1pdm
* intermediates are defined in x_oe_intermediates.c
* 
*   D[i][j] = -LR_oo[j][i] - t1[i][f] * L2R1_ov[j][f]
*
*   D[a][b] = +LR_vv[a][b] + t1[n][b] * L2R1_ov[n][a]
*
*   D[a][i] = +L2R1_ov[i][a]
*
*   D[i][a] = + L1R2_ov[i][a] 
*            - t1[m][a] * LR_oo[m][i] 
*            - t1[i][e] * LR_vv[e][a]
*            - r1[m][a] * LT2_oo[m][i]
*            - r1[i][e] * LT2_vv[e][a]
*            + L2R1_ov[M][E] * (t2[i][m][a][e] - t1[i][e] * t1[m][a])
*
* RAK 2003
*/

void x_onepdm_rohf(struct RHO_Params rho_params)
{
  dpdfile2 DAI, Dai, DIA, Dia, DIJ, DAB, Dij, Dab, TIA, Tia;
  dpdfile2 LIA, Lia, RIA, Ria, I, XIJ, Xij;
  dpdbuf4 T2, L2, R2, I2;
  int L_irr, R_irr, G_irr;
  double dot_IA, dot_ia, dot_AI, dot_ai;
  double dot_IJ;
  L_irr = rho_params.L_irr;
  R_irr = rho_params.R_irr;
  G_irr = rho_params.G_irr;

  dpd_file2_init(&TIA, CC_OEI, 0, 0, 1, "tIA");
  dpd_file2_init(&Tia, CC_OEI, 0, 0, 1, "tia");
  dpd_file2_init(&RIA, CC_GR, R_irr, 0, 1, "RIA");
  dpd_file2_init(&Ria, CC_GR, R_irr, 0, 1, "Ria");
  dpd_file2_init(&LIA, CC_GL, L_irr, 0, 1, "LIA");
  dpd_file2_init(&Lia, CC_GL, L_irr, 0, 1, "Lia");

  /* D[i][j] = -LR_oo[j][i] - t1[i][f] * L2R1_ov[j][f] */
  dpd_file2_init(&DIJ, CC_OEI, G_irr, 0, 0, rho_params.DIJ_lbl);
  dpd_file2_init(&I, EOM_TMP, G_irr, 0, 0, "LR_OO");
  dpd_file2_axpy(&I, &DIJ, -1.0, 1);
  dpd_file2_close(&I);

  if (!params.connect_xi) {
    dpd_file2_init(&I, EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    dpd_contract222(&TIA, &I, &DIJ, 0, 0, -1.0, 1.0);
    dpd_file2_close(&I);
  }
  dpd_file2_close(&DIJ);

  dpd_file2_init(&Dij, CC_OEI, G_irr, 0, 0, rho_params.Dij_lbl);
  dpd_file2_init(&I, EOM_TMP, G_irr, 0, 0, "LR_oo");
  dpd_file2_axpy(&I, &Dij, -1.0, 1);
  dpd_file2_close(&I);

  if (!params.connect_xi) {
    dpd_file2_init(&I, EOM_TMP, G_irr, 0, 1, "L2R1_ov");
    dpd_contract222(&Tia, &I, &Dij, 0, 0, -1.0, 1.0);
    dpd_file2_close(&I);
  }
  dpd_file2_close(&Dij);

  /* D[a][b] = +LR_vv[a][b] + L2R1_ov[n][a] * t1[n][b] */
  dpd_file2_init(&DAB, CC_OEI, G_irr, 1, 1, rho_params.DAB_lbl);
  dpd_file2_init(&I, EOM_TMP, G_irr, 1, 1, "LR_VV");
  dpd_file2_axpy(&I, &DAB, 1.0, 0);
  dpd_file2_close(&I);

  if (!params.connect_xi) {
    dpd_file2_init(&I, EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    dpd_contract222(&I, &TIA, &DAB, 1, 1, 1.0, 1.0);
    dpd_file2_close(&I);
  }
  dpd_file2_close(&DAB);

  dpd_file2_init(&Dab, CC_OEI, G_irr, 1, 1, rho_params.Dab_lbl);
  dpd_file2_init(&I, EOM_TMP, G_irr, 1, 1, "LR_vv");
  dpd_file2_axpy(&I, &Dab, 1.0, 0);
  dpd_file2_close(&I);

  if (!params.connect_xi) {
    dpd_file2_init(&I, EOM_TMP, G_irr, 0, 1, "L2R1_ov");
    dpd_contract222(&I, &Tia, &Dab, 1, 1, 1.0, 1.0);
    dpd_file2_close(&I);
  }
  dpd_file2_close(&Dab);

  /* D[a][i] = +L2R1_ov[i][a] */

  if (!params.connect_xi) {
    dpd_file2_init(&DAI, CC_OEI, G_irr, 0, 1, rho_params.DAI_lbl);
    dpd_file2_init(&I, EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    dpd_file2_axpy(&I, &DAI, 1.0, 0);
    dpd_file2_close(&I);
    dpd_file2_close(&DAI);

    dpd_file2_init(&Dai, CC_OEI, G_irr, 0, 1, rho_params.Dai_lbl);
    dpd_file2_init(&I, EOM_TMP, G_irr, 0, 1, "L2R1_ov");
    dpd_file2_axpy(&I, &Dai, 1.0, 0);
    dpd_file2_close(&I);
    dpd_file2_close(&Dai);
  }

  /*
     D[I][A] = (1-R0)*tIA + L1R2_OV[I][A] 
               - LR_OO[M][I] * t1[M][A]
               - t1[I][E] * LR_vv[E][A]
               - LT2_OO[M][I] * r1[M][A]
               - r1[I][E] * LT2_vv[E][A]
               + L2R1_ov[M][E] * (t2[i][m][a][e] - t1[i][e] * t1[m][a])
  */

  /* (1-R0) * tIA */
  dpd_file2_init(&DIA, CC_OEI, G_irr, 0, 1, rho_params.DIA_lbl);

  if ( (G_irr == 0) /* && (!params.connect_xi)*/ ) {
    dpd_file2_init(&I, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_axpy(&I, &DIA, 1.0, 0);
    dpd_file2_close(&I);
  }

  dpd_file2_init(&Dia, CC_OEI, G_irr, 0, 1, rho_params.Dia_lbl);
  if ( (G_irr == 0) /* && (!params.connect_xi)*/ ) {
    dpd_file2_init(&I, CC_OEI, 0, 0, 1, "tia");
    dpd_file2_axpy(&I, &Dia, 1.0, 0);
    dpd_file2_close(&I);
  }

  /* D[i][a] = L1R2_ov[i][a] */
  dpd_file2_init(&I, EOM_TMP, G_irr, 0, 1, "L1R2_OV");
  dpd_file2_axpy(&I, &DIA, 1.0, 0);
  dpd_file2_close(&I);

  dpd_file2_init(&I, EOM_TMP, G_irr, 0, 1, "L1R2_ov");
  dpd_file2_axpy(&I, &Dia, 1.0, 0);
  dpd_file2_close(&I);

  /* - LR_OO[M][I] * t1[M][A] */

  dpd_file2_init(&I, EOM_TMP, G_irr, 0, 0, "LR_OO");
  dpd_contract222(&I, &TIA, &DIA, 1, 1, -1.0, 1.0);
  dpd_file2_close(&I);

  dpd_file2_init(&I, EOM_TMP, G_irr, 0, 0, "LR_oo");
  dpd_contract222(&I, &Tia, &Dia, 1, 1, -1.0, 1.0);
  dpd_file2_close(&I);

  /* - t1[I][E] * LR_vv[E][A] */

  dpd_file2_init(&I, EOM_TMP, G_irr, 1, 1, "LR_VV");
  dpd_contract222(&TIA, &I, &DIA, 0, 1, -1.0, 1.0);
  dpd_file2_close(&I);

  dpd_file2_init(&I, EOM_TMP, G_irr, 1, 1, "LR_vv");
  dpd_contract222(&Tia, &I, &Dia, 0, 1, -1.0, 1.0);
  dpd_file2_close(&I);

  /* - LT2_OO[M][I] * r1[M][A] */

  dpd_file2_init(&I, EOM_TMP, L_irr, 0, 0, "LT2_OO");
  dpd_contract222(&I, &RIA, &DIA, 1, 1, -1.0, 1.0);
  dpd_file2_close(&I);

  dpd_file2_init(&I, EOM_TMP, L_irr, 0, 0, "LT2_oo");
  dpd_contract222(&I, &Ria, &Dia, 1, 1, -1.0, 1.0);
  dpd_file2_close(&I);

  /* - r1[I][E] * LT2_vv[E][A] */

  dpd_file2_init(&I, EOM_TMP, L_irr, 1, 1, "LT2_VV");
  dpd_contract222(&RIA, &I, &DIA, 0, 1, -1.0, 1.0);
  dpd_file2_close(&I);

  dpd_file2_init(&I, EOM_TMP, L_irr, 1, 1, "LT2_vv");
  dpd_contract222(&Ria, &I, &Dia, 0, 1, -1.0, 1.0);
  dpd_file2_close(&I);

  /* term 6, + L2R1_ov[M][E] * t2[i][m][a][e] */

  if (!params.connect_xi) {
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB"); 
    dpd_file2_init(&I, EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    dpd_dot24(&I, &T2, &DIA, 0, 0, 1.0, 1.0);
    dpd_file2_close(&I);
    dpd_buf4_close(&T2);

    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb"); 
    dpd_file2_init(&I, EOM_TMP, G_irr, 0, 1, "L2R1_ov");
    dpd_dot24(&I, &T2, &DIA, 0, 0, 1.0, 1.0);
    dpd_file2_close(&I);
    dpd_buf4_close(&T2);

    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tijab"); 
    dpd_file2_init(&I, EOM_TMP, G_irr, 0, 1, "L2R1_ov");
    dpd_dot24(&I, &T2, &Dia, 0, 0, 1.0, 1.0);
    dpd_file2_close(&I);
    dpd_buf4_close(&T2);

    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB"); 
    dpd_file2_init(&I, EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    dpd_dot24(&I, &T2, &Dia, 0, 0, 1.0, 1.0);
    dpd_file2_close(&I);
    dpd_buf4_close(&T2);
  }
    
  /* term 7, - (t1[i][e] * L2R1_ov[M][E]) * t1[m][a] */

  if (!params.connect_xi) {
    dpd_file2_init(&I, EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    dpd_file2_init(&XIJ, EOM_TMP, G_irr, 0, 0, "XIJ");
    dpd_contract222(&TIA, &I, &XIJ, 0, 0, 1.0, 0.0);
    dpd_file2_close(&I);

    dpd_file2_init(&XIJ, EOM_TMP, G_irr, 0, 0, "XIJ");
    dpd_contract222(&XIJ, &TIA, &DIA, 0, 1, -1.0, 1.0);
    dpd_file2_close(&XIJ);

    dpd_file2_init(&I, EOM_TMP, G_irr, 0, 1, "L2R1_ov");
    dpd_file2_init(&Xij, EOM_TMP, G_irr, 0, 0, "Xij");
    dpd_contract222(&Tia, &I, &Xij, 0, 0, 1.0, 0.0);
    dpd_file2_close(&I);

    dpd_file2_init(&Xij, EOM_TMP, G_irr, 0, 0, "Xij");
    dpd_contract222(&Xij, &Tia, &Dia, 0, 1, -1.0, 1.0);
    dpd_file2_close(&Xij);
  }

  dpd_file2_close(&DIA);
  dpd_file2_close(&Dia);

  dpd_file2_close(&TIA);
  dpd_file2_close(&Tia);
  dpd_file2_close(&RIA);
  dpd_file2_close(&Ria);
  dpd_file2_close(&LIA);
  dpd_file2_close(&Lia);

  /* check overlaps */
  /*
  dpd_file2_init(&DIA, CC_OEI, G_irr, 0, 1, rho_params.DIA_lbl);
  dot_IA = dpd_file2_dot_self(&DIA);
  dpd_file2_close(&DIA);
  dpd_file2_init(&Dia, CC_OEI, G_irr, 0, 1, rho_params.Dia_lbl);
  dot_ia = dpd_file2_dot_self(&Dia);
  dpd_file2_close(&Dia);
  dpd_file2_init(&DAI, CC_OEI, G_irr, 0, 1, rho_params.DAI_lbl);
  dot_AI = dpd_file2_dot_self(&DAI);
  dpd_file2_close(&DAI);
  dpd_file2_init(&Dai, CC_OEI, G_irr, 0, 1, rho_params.Dai_lbl);
  dot_ai = dpd_file2_dot_self(&Dai);
  dpd_file2_close(&Dai);
  dpd_file2_init(&DIJ, CC_OEI, G_irr, 0, 0, rho_params.DIJ_lbl);
  dot_IJ = dpd_file2_dot_self(&DIJ);
  dpd_file2_close(&DIJ);
  dpd_file2_init(&DIJ, CC_OEI, G_irr, 0, 0, rho_params.Dij_lbl);
  dot_IJ += dpd_file2_dot_self(&DIJ);
  dpd_file2_close(&DIJ); 
  fprintf(outfile,"\tOverlaps of onepdm after excited-state parts added.\n");
  fprintf(outfile,"\t<DIA|DIA> = %15.10lf     <Dia|Dia> = %15.10lf\n", dot_IA, dot_ia);
  fprintf(outfile,"\t<DAI|DAI> = %15.10lf     <Dai|Dai> = %15.10lf\n", dot_AI, dot_ai);
  fprintf(outfile,"\t<Dpq|Dqp> = %15.10lf\n", dot_IA+dot_ia+dot_AI+dot_ai);
  */

  return;
}

}} // namespace psi::ccdensity
