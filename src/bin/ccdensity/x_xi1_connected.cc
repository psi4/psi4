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
    dpd_file2_init(&HIA, EOM_TMP0, G_irr, 0, 1, "HIA");

    dpd_buf4_init(&L2, CC_GL, L_irr, 0, 5, 0, 5, 0, "2LIjAb - LIjbA");
    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
    dpd_dot24(&R1, &L2, &HIA, 0, 0, 1.0, 0.0);
    dpd_file2_close(&R1);
    dpd_buf4_close(&L2);

    tval =  params.cceom_energy;
    /* fprintf(outfile,"\nenergy: %15.10lf\n",tval); */
    dpd_file2_scm(&HIA, tval);

    /* -= (Fme Rme) Lia */
    if (R_irr == 0) {
      dpd_file2_init(&F1, CC_OEI, 0, 0, 1, "FME");
      dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
      dot = 2.0 * dpd_file2_dot(&F1, &R1);
      dpd_file2_close(&R1);
      dpd_file2_close(&F1);

      dpd_file2_init(&L1, CC_GL, L_irr, 0, 1, "LIA");
      dpd_file2_axpy(&L1, &HIA, -dot, 0);
      dpd_file2_close(&L1);
    }

    /* -= - (Rme Lmnea) Fin */
    dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    dpd_file2_init(&F1, CC_OEI, 0, 0, 0, "FMI");
    dpd_contract222(&F1, &I1, &HIA, 0, 1, 1.0, 1.0);
    dpd_file2_close(&F1);
    dpd_file2_close(&I1);

    /* -= Rme Lmief Ffa */
    dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    dpd_file2_init(&F1, CC_OEI, 0, 1, 1, "FAE");
    dpd_contract222(&I1, &F1, &HIA, 0, 1, -1.0, 1.0);
    dpd_file2_close(&F1);
    dpd_file2_close(&I1);

    /* -= Rme Lmnef Wifan */
    dpd_buf4_init(&H2, CC_HBAR, 0, 10, 10, 10, 10, 0, "2 W(ME,jb) + W(Me,Jb)");
    dpd_buf4_sort(&H2, EOM_TMP_XI, prqs, 0, 5, "2 W(ME,jb) + W(Me,Jb) (Mj,Eb)");
    dpd_buf4_close(&H2);

    dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    dpd_buf4_init(&H2, EOM_TMP_XI, 0, 0, 5, 0, 5, 0, "2 W(ME,jb) + W(Me,Jb) (Mj,Eb)");
    dpd_dot24(&I1, &H2, &HIA, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&H2);
    dpd_file2_close(&I1);

    /* -= Limae ( Rmf Fef - Rne Fnm + Rnf Wnefm ) */
    dpd_file2_init(&IME, EOM_TMP_XI, R_irr, 0, 1, "IME");

    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
    dpd_file2_init(&F1, CC_OEI, 0, 1, 1, "FAE");
    dpd_contract222(&R1, &F1, &IME, 0, 0, 1.0, 0.0);
    dpd_file2_close(&F1);
    dpd_file2_close(&R1);

    dpd_file2_init(&F1, CC_OEI, 0, 0, 0, "FMI");
    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
    dpd_contract222(&F1, &R1, &IME, 1, 1, -1.0, 1.0);
    dpd_file2_close(&R1);
    dpd_file2_close(&F1);

    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
    dpd_buf4_init(&H2, EOM_TMP_XI, 0, 0, 5, 0, 5, 0, "2 W(ME,jb) + W(Me,Jb) (Mj,Eb)");
    dpd_dot13(&R1, &H2, &IME, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&H2);
    dpd_file2_close(&R1);

    /* HIA -= LIAME IME + LIAme Ime */
    dpd_buf4_init(&L2, CC_GL, L_irr, 0, 5, 0, 5, 0, "2LIjAb - LIjbA");
    dpd_dot24(&IME, &L2, &HIA, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&L2);

    dpd_file2_close(&IME);

    /* add to Xi1 */
    /* aug_xi_check(&HIA, &Hia); */

    dpd_file2_init(&XIA, EOM_XI, G_irr, 0, 1, "XIA");
    dpd_file2_axpy(&HIA, &XIA, 1.0, 0);
    dpd_file2_close(&XIA);
    dpd_file2_close(&HIA);
  }


  else if (params.ref == 1) { /* ROHF */
    dpd_file2_init(&HIA, EOM_TMP0, G_irr, 0, 1, "HIA");
    dpd_file2_init(&Hia, EOM_TMP0, G_irr, 0, 1, "Hia");

    dpd_buf4_init(&L2, CC_GL, L_irr, 0, 5, 2, 7, 0, "LIJAB");
    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
    dpd_dot24(&R1, &L2, &HIA, 0, 0, 1.0, 0.0);
    dpd_file2_close(&R1);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GL, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "Ria");
    dpd_dot24(&R1, &L2, &HIA, 0, 0, 1.0, 1.0);
    dpd_file2_close(&R1);
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, CC_GL, L_irr, 0, 5, 2, 7, 0, "Lijab");
    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "Ria");
    dpd_dot24(&R1, &L2, &Hia, 0, 0, 1.0, 0.0);
    dpd_file2_close(&R1);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GL, L_irr, 0, 5, 0, 5, 0, "LiJaB");
    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
    dpd_dot24(&R1, &L2, &Hia, 0, 0, 1.0, 1.0);
    dpd_file2_close(&R1);
    dpd_buf4_close(&L2);

    tval =  params.cceom_energy;
    fprintf(outfile,"\nenergy: %15.10lf\n",tval);
    dpd_file2_scm(&HIA, tval);
    dpd_file2_scm(&Hia, tval);

    /* -= (Fme Rme) Lia */
    if (R_irr == 0) {
    dpd_file2_init(&F1, CC_OEI, 0, 0, 1, "FME");
    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
    dot = dpd_file2_dot(&F1, &R1);
    dpd_file2_close(&R1);
    dpd_file2_close(&F1);
    dpd_file2_init(&F1, CC_OEI, 0, 0, 1, "Fme");
    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "Ria");
    dot += dpd_file2_dot(&F1, &R1);
    dpd_file2_close(&R1);
    dpd_file2_close(&F1);

    dpd_file2_init(&L1, CC_GL, L_irr, 0, 1, "LIA");
    dpd_file2_axpy(&L1, &HIA, -dot, 0);
    dpd_file2_close(&L1);
    dpd_file2_init(&L1, CC_GL, L_irr, 0, 1, "Lia");
    dpd_file2_axpy(&L1, &Hia, -dot, 0);
    dpd_file2_close(&L1);
    }

    /* -= - (Rme Lmnea) Fin */
    dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    dpd_file2_init(&F1, CC_OEI, 0, 0, 0, "FMI");
    dpd_contract222(&F1, &I1, &HIA, 0, 1, 1.0, 1.0);
    dpd_file2_close(&F1);
    dpd_file2_close(&I1);

    dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 1, "L2R1_ov");
    dpd_file2_init(&F1, CC_OEI, 0, 0, 0, "Fmi");
    dpd_contract222(&F1, &I1, &Hia, 0, 1, 1.0, 1.0);
    dpd_file2_close(&F1);
    dpd_file2_close(&I1);

    /* -= Rme Lmief Ffa */
    dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    dpd_file2_init(&F1, CC_OEI, 0, 1, 1, "FAE");
    dpd_contract222(&I1, &F1, &HIA, 0, 1, -1.0, 1.0);
    dpd_file2_close(&F1);
    dpd_file2_close(&I1);

    dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 1, "L2R1_ov");
    dpd_file2_init(&F1, CC_OEI, 0, 1, 1, "Fae");
    dpd_contract222(&I1, &F1, &Hia, 0, 1, -1.0, 1.0);
    dpd_file2_close(&F1);
    dpd_file2_close(&I1);

    /* -= Rme Lmnef Wifan */
    dpd_buf4_init(&H2, CC_HBAR, 0, 10, 10, 10, 10, 0, "WMBEJ");
    dpd_buf4_sort(&H2, EOM_TMP_XI, prqs, 0, 5, "WMBEJ (MJ,EB)");
    dpd_buf4_close(&H2);
    dpd_buf4_init(&H2, CC_HBAR, 0, 10, 10, 10, 10, 0, "Wmbej");
    dpd_buf4_sort(&H2, EOM_TMP_XI, prqs, 0, 5, "Wmbej (mj,eb)");
    dpd_buf4_close(&H2);
    dpd_buf4_init(&H2, CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbEj");
    dpd_buf4_sort(&H2, EOM_TMP_XI, prqs, 0, 5, "WMbEj (Mj,Eb)");
    dpd_buf4_close(&H2);
    dpd_buf4_init(&H2, CC_HBAR, 0, 10, 10, 10, 10, 0, "WmBeJ");
    dpd_buf4_sort(&H2, EOM_TMP_XI, prqs, 0, 5, "WmBeJ (mJ,eB)");
    dpd_buf4_close(&H2);

    dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    dpd_buf4_init(&H2, EOM_TMP_XI, 0, 0, 5, 0, 5, 0, "WMBEJ (MJ,EB)");
    dpd_dot24(&I1, &H2, &HIA, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&H2);
    dpd_buf4_init(&H2, EOM_TMP_XI, 0, 0, 5, 0, 5, 0, "WmBeJ (mJ,eB)");
    dpd_dot24(&I1, &H2, &Hia, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&H2);
    dpd_file2_close(&I1);

    dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 1, "L2R1_ov");
    dpd_buf4_init(&H2, EOM_TMP_XI, 0, 0, 5, 0, 5, 0, "Wmbej (mj,eb)");
    dpd_dot24(&I1, &H2, &Hia, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&H2);
    dpd_buf4_init(&H2, EOM_TMP_XI, 0, 0, 5, 0, 5, 0, "WMbEj (Mj,Eb)");
    dpd_dot24(&I1, &H2, &HIA, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&H2);
    dpd_file2_close(&I1);

    /* -= Limae ( Rmf Fef - Rne Fnm + Rnf Wnefm ) */
    dpd_file2_init(&IME, EOM_TMP_XI, R_irr, 0, 1, "IME");
    dpd_file2_init(&Ime, EOM_TMP_XI, R_irr, 0, 1, "Ime");

    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
    dpd_file2_init(&F1, CC_OEI, 0, 1, 1, "FAE");
    dpd_contract222(&R1, &F1, &IME, 0, 0, 1.0, 0.0);
    dpd_file2_close(&F1);
    dpd_file2_close(&R1);
    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "Ria");
    dpd_file2_init(&F1, CC_OEI, 0, 1, 1, "Fae");
    dpd_contract222(&R1, &F1, &Ime, 0, 0, 1.0, 0.0);
    dpd_file2_close(&F1);
    dpd_file2_close(&R1);

    dpd_file2_init(&F1, CC_OEI, 0, 0, 0, "FMI");
    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
    dpd_contract222(&F1, &R1, &IME, 1, 1, -1.0, 1.0);
    dpd_file2_close(&R1);
    dpd_file2_close(&F1);
    dpd_file2_init(&F1, CC_OEI, 0, 0, 0, "Fmi");
    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "Ria");
    dpd_contract222(&F1, &R1, &Ime, 1, 1, -1.0, 1.0);
    dpd_file2_close(&R1);
    dpd_file2_close(&F1);

    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
    dpd_buf4_init(&H2, EOM_TMP_XI, 0, 0, 5, 0, 5, 0, "WMBEJ (MJ,EB)");
    dpd_dot13(&R1, &H2, &IME, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&H2);
    dpd_buf4_init(&H2, EOM_TMP_XI, 0, 0, 5, 0, 5, 0, "WMbEj (Mj,Eb)");
    dpd_dot13(&R1, &H2, &Ime, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&H2);
    dpd_file2_close(&R1);

    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "Ria");
    dpd_buf4_init(&H2, EOM_TMP_XI, 0, 0, 5, 0, 5, 0, "Wmbej (mj,eb)");
    dpd_dot13(&R1, &H2, &Ime, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&H2);
    dpd_buf4_init(&H2, EOM_TMP_XI, 0, 0, 5, 0, 5, 0, "WmBeJ (mJ,eB)");
    dpd_dot13(&R1, &H2, &IME, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&H2);
    dpd_file2_close(&R1);

    /* HIA -= LIAME IME + LIAme Ime */
    dpd_buf4_init(&L2, CC_GL, L_irr, 0, 5, 2, 7, 0, "LIJAB");
    dpd_dot24(&IME, &L2, &HIA, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GL, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    dpd_dot24(&Ime, &L2, &HIA, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GL, L_irr, 0, 5, 2, 7, 0, "Lijab");
    dpd_dot24(&Ime, &L2, &Hia, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GL, L_irr, 0, 5, 0, 5, 0, "LiJaB");
    dpd_dot24(&IME, &L2, &Hia, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&L2);

    dpd_file2_close(&IME);
    dpd_file2_close(&Ime);

    /* add to Xi1 */
    /* aug_xi_check(&HIA, &Hia); */

    dpd_file2_init(&XIA, EOM_XI, G_irr, 0, 1, "XIA");
    dpd_file2_init(&Xia, EOM_XI, G_irr, 0, 1, "Xia");
    dpd_file2_axpy(&HIA, &XIA, 1.0, 0);
    dpd_file2_axpy(&Hia, &Xia, 1.0, 0);
    dpd_file2_close(&XIA);
    dpd_file2_close(&Xia);

    dpd_file2_close(&HIA);
    dpd_file2_close(&Hia);
  }


  else { /* UHF */
    dpd_file2_init(&HIA, EOM_TMP0, G_irr, 0, 1, "HIA");
    dpd_file2_init(&Hia, EOM_TMP0, G_irr, 2, 3, "Hia");

    dpd_buf4_init(&L2, CC_GL, L_irr, 0, 5, 2, 7, 0, "LIJAB");
    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
    dpd_dot24(&R1, &L2, &HIA, 0, 0, 1.0, 0.0);
    dpd_file2_close(&R1);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GL, L_irr, 22, 28, 22, 28, 0, "LIjAb");
    dpd_file2_init(&R1, CC_GR, R_irr, 2, 3, "Ria");
    dpd_dot24(&R1, &L2, &HIA, 0, 0, 1.0, 1.0);
    dpd_file2_close(&R1);
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, CC_GL, L_irr, 10, 15, 12, 17, 0, "Lijab");
    dpd_file2_init(&R1, CC_GR, R_irr, 2, 3, "Ria");
    dpd_dot24(&R1, &L2, &Hia, 0, 0, 1.0, 0.0);
    dpd_file2_close(&R1);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GL, L_irr, 23, 29, 23, 29, 0, "LiJaB");
    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
    dpd_dot24(&R1, &L2, &Hia, 0, 0, 1.0, 1.0);
    dpd_file2_close(&R1);
    dpd_buf4_close(&L2);

    tval =  params.cceom_energy;
    fprintf(outfile,"\nenergy: %15.10lf\n",tval);
    dpd_file2_scm(&HIA, tval);
    dpd_file2_scm(&Hia, tval);

    /* -= (Fme Rme) Lia */
    if (R_irr == 0) {
    dpd_file2_init(&F1, CC_OEI, 0, 0, 1, "FME");
    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
    dot = dpd_file2_dot(&F1, &R1);
    dpd_file2_close(&R1);
    dpd_file2_close(&F1);
    dpd_file2_init(&F1, CC_OEI, 0, 2, 3, "Fme");
    dpd_file2_init(&R1, CC_GR, R_irr, 2, 3, "Ria");
    dot += dpd_file2_dot(&F1, &R1);
    dpd_file2_close(&R1);
    dpd_file2_close(&F1);

    dpd_file2_init(&L1, CC_GL, L_irr, 0, 1, "LIA");
    dpd_file2_axpy(&L1, &HIA, -dot, 0);
    dpd_file2_close(&L1);
    dpd_file2_init(&L1, CC_GL, L_irr, 2, 3, "Lia");
    dpd_file2_axpy(&L1, &Hia, -dot, 0);
    dpd_file2_close(&L1);
    }

    /* -= - (Rme Lmnea) Fin */
    dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    dpd_file2_init(&F1, CC_OEI, 0, 0, 0, "FMI");
    dpd_contract222(&F1, &I1, &HIA, 0, 1, 1.0, 1.0);
    dpd_file2_close(&F1);
    dpd_file2_close(&I1);

    dpd_file2_init(&I1, EOM_TMP, G_irr, 2, 3, "L2R1_ov");
    dpd_file2_init(&F1, CC_OEI, 0, 2, 2, "Fmi");
    dpd_contract222(&F1, &I1, &Hia, 0, 1, 1.0, 1.0);
    dpd_file2_close(&F1);
    dpd_file2_close(&I1);

    /* -= Rme Lmief Ffa */
    dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    dpd_file2_init(&F1, CC_OEI, 0, 1, 1, "FAE");
    dpd_contract222(&I1, &F1, &HIA, 0, 1, -1.0, 1.0);
    dpd_file2_close(&F1);
    dpd_file2_close(&I1);

    dpd_file2_init(&I1, EOM_TMP, G_irr, 2, 3, "L2R1_ov");
    dpd_file2_init(&F1, CC_OEI, 0, 3, 3, "Fae");
    dpd_contract222(&I1, &F1, &Hia, 0, 1, -1.0, 1.0);
    dpd_file2_close(&F1);
    dpd_file2_close(&I1);

    /* -= Rme Lmnef Wifan */
    dpd_buf4_init(&H2, CC_HBAR, 0, 20, 20, 20, 20, 0, "WMBEJ");
    dpd_buf4_sort(&H2, EOM_TMP_XI, prqs, 0, 5, "WMBEJ (MJ,EB)");
    dpd_buf4_close(&H2);
    dpd_buf4_init(&H2, CC_HBAR, 0, 30, 30, 30, 30, 0, "Wmbej");
    dpd_buf4_sort(&H2, EOM_TMP_XI, prqs, 10, 15, "Wmbej (mj,eb)");
    dpd_buf4_close(&H2);
    dpd_buf4_init(&H2, CC_HBAR, 0, 20, 30, 20, 30, 0, "WMbEj");
    dpd_buf4_sort(&H2, EOM_TMP_XI, prqs, 22, 28, "WMbEj (Mj,Eb)");
    dpd_buf4_close(&H2);
    dpd_buf4_init(&H2, CC_HBAR, 0, 30, 20, 30, 20, 0, "WmBeJ");
    dpd_buf4_sort(&H2, EOM_TMP_XI, prqs, 23, 29, "WmBeJ (mJ,eB)");
    dpd_buf4_close(&H2);

    dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    dpd_buf4_init(&H2, EOM_TMP_XI, 0, 0, 5, 0, 5, 0, "WMBEJ (MJ,EB)");
    dpd_dot24(&I1, &H2, &HIA, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&H2);
    dpd_buf4_init(&H2, EOM_TMP_XI, 0, 23, 29, 23, 29, 0, "WmBeJ (mJ,eB)");
    dpd_dot24(&I1, &H2, &Hia, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&H2);
    dpd_file2_close(&I1);
    dpd_file2_init(&I1, EOM_TMP, G_irr, 2, 3, "L2R1_ov");
    dpd_buf4_init(&H2, EOM_TMP_XI, 0, 10, 15, 10, 15, 0, "Wmbej (mj,eb)");
    dpd_dot24(&I1, &H2, &Hia, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&H2);
    dpd_buf4_init(&H2, EOM_TMP_XI, 0, 22, 28, 22, 28, 0, "WMbEj (Mj,Eb)");
    dpd_dot24(&I1, &H2, &HIA, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&H2);
    dpd_file2_close(&I1);

    /* -= Limae ( Rmf Fef - Rne Fnm + Rnf Wnefm ) */
    dpd_file2_init(&IME, EOM_TMP_XI, R_irr, 0, 1, "IME");
    dpd_file2_init(&Ime, EOM_TMP_XI, R_irr, 2, 3, "Ime");

    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
    dpd_file2_init(&F1, CC_OEI, 0, 1, 1, "FAE");
    dpd_contract222(&R1, &F1, &IME, 0, 0, 1.0, 0.0);
    dpd_file2_close(&F1);
    dpd_file2_close(&R1);
    dpd_file2_init(&R1, CC_GR, R_irr, 2, 3, "Ria");
    dpd_file2_init(&F1, CC_OEI, 0, 3, 3, "Fae");
    dpd_contract222(&R1, &F1, &Ime, 0, 0, 1.0, 0.0);
    dpd_file2_close(&F1);
    dpd_file2_close(&R1);

    dpd_file2_init(&F1, CC_OEI, 0, 0, 0, "FMI");
    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
    dpd_contract222(&F1, &R1, &IME, 1, 1, -1.0, 1.0);
    dpd_file2_close(&R1);
    dpd_file2_close(&F1);
    dpd_file2_init(&F1, CC_OEI, 0, 2, 2, "Fmi");
    dpd_file2_init(&R1, CC_GR, R_irr, 2, 3, "Ria");
    dpd_contract222(&F1, &R1, &Ime, 1, 1, -1.0, 1.0);
    dpd_file2_close(&R1);
    dpd_file2_close(&F1);

    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
    dpd_buf4_init(&H2, EOM_TMP_XI, 0, 0, 5, 0, 5, 0, "WMBEJ (MJ,EB)");
    dpd_dot13(&R1, &H2, &IME, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&H2);
    dpd_buf4_init(&H2, EOM_TMP_XI, 0, 22, 28, 22, 28, 0, "WMbEj (Mj,Eb)");
    dpd_dot13(&R1, &H2, &Ime, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&H2);
    dpd_file2_close(&R1);
    dpd_file2_init(&R1, CC_GR, R_irr, 2, 3, "Ria");
    dpd_buf4_init(&H2, EOM_TMP_XI, 0, 10, 15, 10, 15, 0, "Wmbej (mj,eb)");
    dpd_dot13(&R1, &H2, &Ime, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&H2);
    dpd_buf4_init(&H2, EOM_TMP_XI, 0, 23, 29, 23, 29, 0, "WmBeJ (mJ,eB)");
    dpd_dot13(&R1, &H2, &IME, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&H2);
    dpd_file2_close(&R1);

    /* HIA -= LIAME IME + LIAme Ime */
    dpd_buf4_init(&L2, CC_GL, L_irr, 0, 5, 2, 7, 0, "LIJAB");
    dpd_dot24(&IME, &L2, &HIA, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GL, L_irr, 22, 28, 22, 28, 0, "LIjAb");
    dpd_dot24(&Ime, &L2, &HIA, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GL, L_irr, 10, 15, 12, 17, 0, "Lijab");
    dpd_dot24(&Ime, &L2, &Hia, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GL, L_irr, 23, 29, 23, 29, 0, "LiJaB");
    dpd_dot24(&IME, &L2, &Hia, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&L2);

    dpd_file2_close(&IME);
    dpd_file2_close(&Ime);

    /* aug_xi_check(&HIA, &Hia); */

    /* add to Xi1 */
    dpd_file2_init(&XIA, EOM_XI, G_irr, 0, 1, "XIA");
    dpd_file2_init(&Xia, EOM_XI, G_irr, 2, 3, "Xia");
    dpd_file2_axpy(&HIA, &XIA, 1.0, 0);
    dpd_file2_axpy(&Hia, &Xia, 1.0, 0);
    dpd_file2_close(&XIA);
    dpd_file2_close(&Xia);

    dpd_file2_close(&HIA);
    dpd_file2_close(&Hia);
  }
}

/*
double aug_xi_check(dpdfile2 *HIA, dpdfile2 *Hia) 
{
  double tvalA, tvalB;
  tvalA = tvalB = 0.0;

  if (params.ref == 0) {
    tvalA = dpd_file2_dot_self(HIA);
    fprintf(outfile, "<HIA|HIA> = %15.10lf\n", tvalA);
  }
  else {
    tvalB = dpd_file2_dot_self(Hia);
    fprintf(outfile, "<HIA|HIA> = %15.10lf\n", tvalA);
  }
  fprintf(outfile, "<H1|H1> = %15.10lf\n", tvalA + tvalB);
  return;
}
*/

}} // namespace psi::ccdensity
