/*! \file
    \ingroup CCEOM
    \brief Enter brief description of file here 
*/
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cceom {

/* compute [ [H,eT1] , C1 ]
 *
 * [H,eT^1] matrix elements are in CC3_HET1
 *
 * C1 has symmetry C_irr and is read from EOM_CME with name "CME i"
 *
 * illustrative code that computes [H,C1] is in cc3_HC1.c
 *
 * output matrix elements are written to CC_HC1ET1
 *
*/
void HC1ET1_Wmbij(int i, int C_irr);
void HC1ET1_Wabei(int i, int C_irr);
void HC1ET1_Wmbij_rhf(int i, int C_irr);
void HC1ET1_Wabei_rhf(int i, int C_irr);

void cc3_HC1ET1 (int i, int C_irr) {

  /* only need new Wmbij and Wabei for CC3 EOM energies */

  if (params.ref == 0) {
    HC1ET1_Wmbij_rhf(i, C_irr);
    HC1ET1_Wabei_rhf(i, C_irr);
  }
  else {
    HC1ET1_Wmbij(i, C_irr);
    HC1ET1_Wabei(i, C_irr);
  }

  return;
}

void HC1ET1_Wmbij(int i, int C_irr)
{
  double dot;
  dpdbuf4 C, D, E, F, Ht, W, W1, X, Z;
  dpdfile2 CME, Cme;
  char CME_lbl[32], Cme_lbl[32];
  sprintf(CME_lbl, "%s %d", "CME", i);
  sprintf(Cme_lbl, "%s %d", "Cme", i);

  if (params.ref == 1) { /** ROHF **/
    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);
    dpd_file2_init(&Cme, EOM_Cme, C_irr, 0, 1, Cme_lbl);

    dpd_file2_close(&CME);
    dpd_file2_close(&Cme);
  }
  else if (params.ref == 2) { /** UHF **/

    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);
    dpd_file2_init(&Cme, EOM_Cme, C_irr, 2, 3, Cme_lbl);

    /**** Term I ****/

    /***** Ht (MB,I>J) <--- -WMNIJ * CNB *****/
    dpd_buf4_init(&Ht, CC3_HC1ET1, C_irr, 20, 2, 20, 2, 0, "Ht_WMBIJ (MB,I>J)");
    dpd_buf4_init(&W, CC3_HET1, 0, 0, 2, 2, 2, 0, "CC3 WMNIJ (M>N,I>J)");
    dpd_contract424(&W, &CME, &Ht, 1, 0, 1, -1.0, 0.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Ht);
 
    /***** Ht (mb,i>j) <--- -Wmnij * Cnb *****/
    dpd_buf4_init(&Ht, CC3_HC1ET1, C_irr, 30, 12, 30, 12, 0, "Ht_Wmbij (mb,i>j)");
    dpd_buf4_init(&W, CC3_HET1, 0, 10, 12, 12, 12, 0, "CC3 Wmnij (m>n,i>j)");
    dpd_contract424(&W, &Cme, &Ht, 1, 0, 1, -1.0, 0.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Ht);

    /***** Ht (Mb,Ij) <--- -WMnIj * Cnb *****/
    dpd_buf4_init(&Ht, CC3_HC1ET1, C_irr, 24, 22, 24, 22, 0, "Ht_WMbIj (Mb,Ij)");
    dpd_buf4_init(&W, CC3_HET1, 0, 22, 22, 22, 22, 0, "CC3 WMnIj (Mn,Ij)");
    dpd_contract424(&W, &Cme, &Ht, 1, 0, 1, -1.0, 0.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Ht);

    /***** Ht (mB,iJ) <--- ZBmJi <--- CNB * WNmJi *****/
    dpd_buf4_init(&Z, CC_TMP0, C_irr, 26, 22, 26, 22, 0, "Z (Bm,Ji)");
    dpd_buf4_init(&W, CC3_HET1, 0, 22, 22, 22, 22, 0, "CC3 WMnIj (Mn,Ij)");
    dpd_contract244(&CME, &W, &Z, 0, 0, 0, -1.0, 0.0);
    dpd_buf4_close(&W);
    dpd_buf4_sort(&Z, CC3_HC1ET1, qpsr, 27, 23, "Ht_WmBiJ (mB,iJ)");
    dpd_buf4_close(&Z);
		
    /**** Term II ****/

    /***** Ht (MB,I>J) <--- -P(I/J) X (MB,JI) <--- WMBEJ(ME,JB) * CIE *****/
    dpd_buf4_init(&Z, CC_TMP0, C_irr, 0, 20, 0, 20, 0, "Z (MI,JB)");
    dpd_buf4_init(&W, CC3_HET1, 0, 20, 20, 20, 20, 0, "CC3 WMBEJ (ME,JB)");
    dpd_contract424(&W, &CME, &Z, 1, 1, 1, 1.0, 0.0);
    dpd_buf4_close(&W);

    dpd_buf4_sort(&Z, CC_TMP0, psrq, 20, 0, "X (MB,JI)");
    dpd_buf4_close(&Z);
    dpd_buf4_init(&X, CC_TMP0, C_irr, 20, 0, 20, 0, 0, "X (MB,JI)");
    dpd_buf4_sort_axpy(&X, CC3_HC1ET1, pqsr, 20, 2, "Ht_WMBIJ (MB,I>J)", 1);

    dpd_buf4_init(&Ht, CC3_HC1ET1, C_irr, 20, 0, 20, 2, 0, "Ht_WMBIJ (MB,I>J)");
    dpd_buf4_axpy(&X, &Ht, -1);
    dpd_buf4_close(&X);
    dpd_buf4_close(&Ht);

    /***** Ht (mb,i>j) <--- -P(i/j) X (mb,ij) <--- Wmbej (me,jb) * Cie *****/

    dpd_buf4_init(&Z, CC_TMP0, C_irr, 10, 30, 10, 30, 0, "Z (mi,jb)");
    dpd_buf4_init(&W, CC3_HET1, 0, 30, 30, 30, 30, 0, "CC3 Wmbej (me,jb)");
    dpd_contract424(&W, &Cme, &Z, 1, 1, 1, 1.0, 0.0);
    dpd_buf4_close(&W);

    dpd_buf4_sort(&Z, CC_TMP0, psrq, 30, 10, "X (mb,ji)");
    dpd_buf4_close(&Z);
    dpd_buf4_init(&X, CC_TMP0, C_irr, 30, 10, 30, 10, 0, "X (mb,ji)");
    dpd_buf4_sort_axpy(&X, CC3_HC1ET1, pqsr, 30, 12, "Ht_Wmbij (mb,i>j)", 1);

    dpd_buf4_init(&Ht, CC3_HC1ET1, C_irr, 30, 10, 30, 12, 0, "Ht_Wmbij (mb,i>j)");
    dpd_buf4_axpy(&X, &Ht, -1);
    dpd_buf4_close(&X);
    dpd_buf4_close(&Ht);

    /***** Ht (Mb,Ij) <--- CIE * WMbEj *****/
    dpd_buf4_init(&Z, CC_TMP0, C_irr, 0, 30, 0, 30, 0, "Z (MI,jb)");
    dpd_buf4_init(&W, CC3_HET1, 0, 20, 30, 20, 30, 0, "CC3 WMbEj (ME,jb)");
    dpd_contract424(&W, &CME, &Z, 1, 1, 1, 1.0, 0.0);
    dpd_buf4_close(&W);
    dpd_buf4_sort_axpy(&Z, CC3_HC1ET1, psqr, 24, 22, "Ht_WMbIj (Mb,Ij)", 1);
    dpd_buf4_close(&Z);

    /***** Ht (Mb,Ij) <--- Cje * WMbeI *****/
    dpd_buf4_init(&Z, CC_TMP0, C_irr, 22, 24, 22, 24, 0, "Z (Mj,Ib)");
    dpd_buf4_init(&W, CC3_HET1, 0, 24, 24, 24, 24, 0, "CC3 WMbeJ (Me,Jb)");
    dpd_contract424(&W, &Cme, &Z, 1, 1, 1, 1.0, 0.0);
    dpd_buf4_close(&W);
    dpd_buf4_sort_axpy(&Z, CC3_HC1ET1, psrq, 24, 22, "Ht_WMbIj (Mb,Ij)", -1);
    dpd_buf4_close(&Z);

    /***** Ht (mB,iJ) <--- Cie * WmBiJ *****/
    dpd_buf4_init(&Z, CC_TMP0, C_irr, 10, 20, 10, 20, 0, "Z (mi,JB)");
    dpd_buf4_init(&W, CC3_HET1, 0, 30, 20, 30, 20, 0, "CC3 WmBeJ (me,JB)");
    dpd_contract424(&W, &Cme, &Z, 1, 1, 1, 1.0, 0.0);
    dpd_buf4_close(&W);
    dpd_buf4_sort_axpy(&Z, CC3_HC1ET1, psqr, 27, 23, "Ht_WmBiJ (mB,iJ)", 1);
    dpd_buf4_close(&Z);

    /***** Ht (mB,iJ) <--- CJE * WmEiB *****/
    dpd_buf4_init(&Z, CC_TMP0, C_irr, 23, 27, 23, 27, 0, "Z (mJ,iB)");
    dpd_buf4_init(&W, CC3_HET1, 0, 27, 27, 27, 27, 0, "CC3 WmBEj (mE,jB)");
    dpd_contract424(&W, &CME, &Z, 1, 1, 1, 1.0, 0.0);
    dpd_buf4_close(&W);
    dpd_buf4_sort_axpy(&Z, CC3_HC1ET1, psrq, 27, 23, "Ht_WmBiJ (mB,iJ)", -1);
    dpd_buf4_close(&Z);


    dpd_buf4_init(&Ht, CC3_HC1ET1, C_irr, 20, 2, 20, 2, 0, "Ht_WMBIJ (MB,I>J)");
    dpd_buf4_sort(&Ht, CC3_HC1ET1, rspq, 2, 20, "Ht_WMBIJ (I>J,MB)");
    dpd_buf4_close(&Ht);

    dpd_buf4_init(&Ht, CC3_HC1ET1, C_irr, 30, 12, 30, 12, 0, "Ht_Wmbij (mb,i>j)");
    dpd_buf4_sort(&Ht, CC3_HC1ET1, rspq, 12, 30, "Ht_Wmbij (i>j,mb)");
    dpd_buf4_close(&Ht);

    dpd_buf4_init(&Ht, CC3_HC1ET1, C_irr, 24, 22, 24, 22, 0, "Ht_WMbIj (Mb,Ij)");
    dpd_buf4_sort(&Ht, CC3_HC1ET1, rspq, 22, 24, "Ht_WMbIj (Ij,Mb)");
    dpd_buf4_close(&Ht);

    dpd_buf4_init(&Ht, CC3_HC1ET1, C_irr, 27, 23, 27, 23, 0, "Ht_WmBiJ (mB,iJ)");
    dpd_buf4_sort(&Ht, CC3_HC1ET1, rspq, 23, 27, "Ht_WmBiJ (iJ,mB)");
    dpd_buf4_close(&Ht);


    /************ TEST *************/

    /*
    dpd_buf4_init(&W, CC3_HC1ET1, 0, 0, 20, 2, 20, 0, "Ht_WMBIJ (I>J,MB)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    fprintf(outfile,"\t<WMBIJ (I>J,MB)|WMBIJ> = %15.10lf\n", dot);

    dpd_buf4_init(&W, CC3_HC1ET1, 0, 10, 30, 12, 30, 0, "Ht_Wmbij (i>j,mb)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    fprintf(outfile,"\t<Wmbij (i>j,mb)|Wmbij> = %15.10lf\n", dot);

    dpd_buf4_init(&W, CC3_HC1ET1, 0, 22, 24, 22, 24, 0, "Ht_WMbIj (Ij,Mb)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    fprintf(outfile,"\t<WMbIj (Ij,Mb)|WMbIj> = %15.10lf\n", dot);

    dpd_buf4_init(&W, CC3_HC1ET1, 0, 23, 27, 23, 27, 0, "Ht_WmBiJ (iJ,mB)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    fprintf(outfile,"\t<WmBiJ (iJ,mB)|WmBiJ> = %15.10lf\n", dot);
    */

    dpd_file2_close(&CME);
    dpd_file2_close(&Cme);
  }

  return;
}

void HC1ET1_Wabei(int i, int C_irr)
{
  double dot;
  dpdfile2 CME, Cme, tIA, tia;
  dpdbuf4 Ht, Z, Z1, Z2, Z3, B, C, D, E, F, W, X;

  char CME_lbl[32], Cme_lbl[32];
  sprintf(CME_lbl, "%s %d", "CME", i);
  sprintf(Cme_lbl, "%s %d", "Cme", i);

  if (params.ref == 1) { /** ROHF **/
    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);
    dpd_file2_init(&Cme, EOM_Cme, C_irr, 0, 1, Cme_lbl);

    dpd_file2_close(&CME);
    dpd_file2_close(&Cme);
  }
  else if (params.ref == 2) { /** UHF **/
    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);
    dpd_file2_init(&Cme, EOM_Cme, C_irr, 2, 3, Cme_lbl);
    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 2, 3, "tia");


    /**** Term I ****/

    /***** Ht_WABEI <--- -P(A/B) CMA * WMBEI *****/

    dpd_buf4_init(&Z, CC_TMP0, C_irr, 5, 21, 5, 21, 0, "Z (AB,EI)");
    dpd_buf4_init(&W, CC3_HET1, 0, 20, 21, 20, 21, 0, "CC3 WMBEJ (MB,EJ)");
    dpd_contract244(&CME, &W, &Z, 0, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&W);

    dpd_buf4_sort(&Z, CC_TMP0, qprs, 7, 21, "Ht_WABEI (A>B,EI)");
    dpd_buf4_init(&Ht, CC_TMP0, C_irr, 5, 21, 7, 21, 0, "Ht_WABEI (A>B,EI)");
    dpd_buf4_axpy(&Z, &Ht, -1);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&Ht);

    /***** Ht_Wabei <--- Xbaei <--- -P(a/b) Zabei <-- Cma * Wmbei *****/

    dpd_buf4_init(&Z, CC_TMP0, C_irr, 15, 31, 15, 31, 0, "Z (ab,ei)");
    dpd_buf4_init(&W, CC3_HET1, 0, 30, 31, 30, 31, 0, "CC3 Wmbej (mb,ej)");
    dpd_contract244(&Cme, &W, &Z, 0, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&W);

    dpd_buf4_sort(&Z, CC_TMP0, qprs, 17, 31, "Ht_Wabei (a>b,ei)");
    dpd_buf4_init(&Ht, CC_TMP0, C_irr, 15, 31, 17, 31, 0, "Ht_Wabei (a>b,ei)");
    dpd_buf4_axpy(&Z, &Ht, -1);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&Ht);

    /***** Ht_WAbEi <--- -CMA * WMbEi *****/
    dpd_buf4_init(&Ht, CC_TMP0, C_irr, 28, 26, 28, 26, 0, "Ht_WAbEi (Ab,Ei)");
    dpd_buf4_init(&W, CC3_HET1, 0, 24, 26, 24, 26, 0, "CC3 WMbEj (Mb,Ej)");
    dpd_contract244(&CME, &W, &Ht, 0, 0, 0, -1.0, 0.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Ht);

    /***** Ht_WAbEi <--- WAmEi * Cmb *****/
    dpd_buf4_init(&Ht, CC_TMP0, C_irr, 26, 28, 26, 28, 0, "Ht_WAbEi (Ei,Ab)");
    dpd_buf4_init(&W, CC3_HET1, 0, 26, 26, 26, 26, 0, "CC3 WmBEj (Bm,Ej)");
    dpd_contract424(&W, &Cme, &Ht, 1, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Ht);

    /***** Ht_WaBeI <--- -Cma * WmBeI *****/
    dpd_buf4_init(&Ht, CC_TMP0, C_irr, 29, 25, 29, 25, 0, "Ht_WaBeI (aB,eI)");
    dpd_buf4_init(&W, CC3_HET1, 0, 27, 25, 27, 25, 0, "CC3 WmBeJ (mB,eJ)");
    dpd_contract244(&Cme, &W, &Ht, 0, 0, 0, -1.0, 0.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Ht);

    /***** Ht_WaBeI <--- WaMeI * CMB *****/
    dpd_buf4_init(&Ht, CC_TMP0, C_irr, 25, 29, 25, 29, 0, "Ht_WaBeI (eI,aB)");
    dpd_buf4_init(&W, CC3_HET1, 0, 25, 25, 25, 25, 0, "CC3 WMbeJ (bM,eJ)");
    dpd_contract424(&W, &CME, &Ht, 1, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Ht);

    /**** Term II ****/

    /***** Ht_WABEI <--- <AB||EF> * CIF *****/
    dpd_buf4_init(&Ht, CC_TMP0, C_irr, 7, 21, 7, 21, 0, "Ht_WABEI (A>B,EI)");
    dpd_buf4_init(&B, CC_BINTS, 0, 7, 5, 5, 5, 1, "B <AB|CD>");
    dpd_contract424(&B, &CME, &Ht, 3, 1, 0, 1.0, 1.0);
    dpd_buf4_close(&B);
    dpd_buf4_close(&Ht);


    /***** Ht_Wabei <--- <ab||ef> * Cif *****/
    dpd_buf4_init(&Ht, CC_TMP0, C_irr, 17, 31, 17, 31, 0, "Ht_Wabei (a>b,ei)");
    dpd_buf4_init(&B, CC_BINTS, 0, 17, 15, 15, 15, 1, "B <ab|cd>");
    dpd_contract424(&B, &Cme, &Ht, 3, 1, 0, 1.0, 1.0);
    dpd_buf4_close(&B);
    dpd_buf4_close(&Ht);

    /***** HAbEi <--- <Ab|Ef> * Cif *****/
    dpd_buf4_init(&Ht, CC_TMP0, C_irr, 28, 26, 28, 26, 0, "Ht_WAbEi (Ab,Ei)");
    dpd_buf4_init(&B, CC_BINTS, 0, 28, 28, 28, 28, 0, "B <Ab|Cd>");
    dpd_contract424(&B, &Cme, &Ht, 3, 1, 0, 1.0, 1.0);
    dpd_buf4_close(&B);
    dpd_buf4_close(&Ht);

    /***** HaBeI <--- C[I][F] * <aB|eF>  *****/
    dpd_buf4_init(&Z, CC_TMP0, C_irr, 24, 28, 24, 28, 0, "Z (Ie,Ba)");
    dpd_buf4_init(&B, CC_BINTS, 0, 28, 28, 28, 28, 0, "B <Ab|Cd>");
    dpd_contract244(&CME, &B, &Z, 1, 0, 0, 1, 0);
    dpd_buf4_close(&B);
    dpd_buf4_sort_axpy(&Z, CC_TMP0, qpsr, 25, 29, "Ht_WaBeI (eI,aB)", 1);
    dpd_buf4_close(&Z);


    /** term III **/

    /** H( A>B,EI) <--- -P(A/B)( <AM||EF>*C[I][F]*t1[M][B] ) **/
    dpd_buf4_init(&Z, CC_TMP0, C_irr, 20, 21, 20, 21, 0, "Z (MA,EI)");
    dpd_buf4_init(&F, CC_FINTS, 0, 20, 5, 20, 5, 1, "F <IA|BC>");
    dpd_contract424(&F, &CME, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&F);

    dpd_buf4_init(&Z1, CC_TMP0, C_irr, 5, 21, 5, 21, 0, "Z1(BA,EI)");
    dpd_contract244(&tIA, &Z, &Z1, 0, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&Z);

    dpd_buf4_sort_axpy(&Z1, CC_TMP0, qprs, 7, 21, "Ht_WABEI (A>B,EI)", 1);
    dpd_buf4_init(&Ht, CC_TMP0, C_irr, 5, 21, 7, 21, 0, "Ht_WABEI (A>B,EI)");
    dpd_buf4_axpy(&Z1, &Ht, -1);
    dpd_buf4_close(&Ht);
    dpd_buf4_close(&Z1);


    /** H(a>b,ei) <--- -P(a/b)( <am||ef>*C[i][f]*t1[m][b] ) **/
    dpd_buf4_init(&Z, CC_TMP0, C_irr, 30, 31, 30, 31, 0, "Z (ma,ei)");
    dpd_buf4_init(&F, CC_FINTS, 0, 30, 15, 30, 15, 1, "F <ia|bc>");
    dpd_contract424(&F, &Cme, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&F);

    dpd_buf4_init(&Z1, CC_TMP0, C_irr, 15, 31, 15, 31, 0, "Z1(ba,ei)");
    dpd_contract244(&tia, &Z, &Z1, 0, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&Z);

    dpd_buf4_sort_axpy(&Z1, CC_TMP0, qprs, 17, 31, "Ht_Wabei (a>b,ei)", 1);
    dpd_buf4_init(&Ht, CC_TMP0, C_irr, 15, 31, 17, 31, 0, "Ht_Wabei (a>b,ei)");
    dpd_buf4_axpy(&Z1, &Ht, -1);
    dpd_buf4_close(&Ht);
    dpd_buf4_close(&Z1);

    /** H(Ab,Ei) <--- -<Am|Ef>*C[i][f]*t1[m][b] **/

    dpd_buf4_init(&Z, CC_TMP0, C_irr, 26, 26, 26, 26, 0, "Z(Am,Ei)");
    dpd_buf4_init(&F, CC_FINTS, 0, 26, 28, 26, 28, 0, "F <Ai|Bc>");
    dpd_contract424(&F, &Cme, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&F);

    dpd_buf4_init(&Ht, CC_TMP0, C_irr, 26, 28, 26, 28, 0, "Ht_WAbEi (Ei,Ab)");
    dpd_contract424(&Z, &tia, &Ht, 1, 0, 0, -1, 1);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&Ht);

    /** Ht(Ab,Ei) <--- -<Mb|Ef>*C[i][f]*t1[M][A] **/

    dpd_buf4_init(&Z, CC_TMP0, C_irr, 24, 26, 24, 26, 0, "Z(Mb,Ei)");
    dpd_buf4_init(&F, CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
    dpd_contract424(&F, &Cme, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&F);

    dpd_buf4_init(&Ht, CC_TMP0, C_irr, 28, 26, 28, 26, 0, "Ht_WAbEi (Ab,Ei)");
    dpd_contract244(&tIA, &Z, &Ht, 0, 0, 0, -1, 1);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&Ht);

    /** Ht(aB,eI) <--- -<aM|eF>*C[I][F]*t1[M][B] **/

    dpd_buf4_init(&Z, CC_TMP0, C_irr, 25, 25, 25, 25, 0, "Z(aM,eI)");
    dpd_buf4_init(&F, CC_FINTS, 0, 25, 29, 25, 29, 0, "F <aI|bC>");
    dpd_contract424(&F, &CME, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&F);

    dpd_buf4_init(&Ht, CC_TMP0, C_irr, 25, 29, 25, 29, 0, "Ht_WaBeI (eI,aB)");
    dpd_contract424(&Z, &tIA, &Ht, 1, 0, 0, -1, 1);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&Ht);

    /** Ht(aB,eI) <--- -C[I][F] * <mB|eF> * t1[m][a] **/

    dpd_buf4_init(&Z, CC_TMP0, C_irr, 27, 25, 27, 25, 0, "Z(mB,eI)");
    dpd_buf4_init(&F, CC_FINTS, 0, 27, 29, 27, 29, 0, "F <iA|bC>");
    dpd_contract424(&F, &CME, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&F);

    dpd_buf4_init(&Ht, CC_TMP0, C_irr, 29, 25, 29, 25, 0, "Ht_WaBeI (aB,eI)");
    dpd_contract244(&tia, &Z, &Ht, 0, 0, 0, -1, 1);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&Ht);


    /** term 6 **/

    /** Ht(A>B,EI) <---  0.5 * P(A/B) <MN||EF> * t1[M][A]*C[I][F]*t1[N][B] **/
    dpd_buf4_init(&Z, CC_TMP0, C_irr, 2, 21, 2, 21, 0, "Z (M>N,EI)");
    dpd_buf4_init(&D, CC_DINTS, 0, 2, 5, 2, 5, 0, "D <IJ||AB> (I>J,AB)");
    dpd_contract424(&D, &CME, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&Z);

    dpd_buf4_init(&Z, CC_TMP0, C_irr, 0, 21, 2, 21, 0, "Z (M>N,EI)");
    dpd_buf4_init(&Z1, CC_TMP0, C_irr, 21, 20, 21, 20, 0, "Z1 (EI,MB)");
    dpd_contract424(&Z, &tIA, &Z1, 1, 0, 0, 1, 0);
    dpd_buf4_close(&Z);

    dpd_buf4_init(&Ht, CC_TMP0, C_irr, 5, 21, 7, 21, 0, "Ht_WABEI (A>B,EI)");
    dpd_contract244(&tIA, &Z1, &Ht, 0, 2, 0, 1, 1);
    dpd_buf4_close(&Ht);
    dpd_buf4_close(&Z1);

    /** W(a>b,ei) <---  0.5 * P(a/b) <nm||ef> * t1[m][b]*C[i][f]*t1[n][a] **/
    dpd_buf4_init(&Z, CC_TMP0, C_irr, 12, 31, 12, 31, 0, "Z (m>n,ei)");
    dpd_buf4_init(&D, CC_DINTS, 0, 12, 15, 12, 15, 0, "D <ij||ab> (i>j,ab)");
    dpd_contract424(&D, &Cme, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&Z);

    dpd_buf4_init(&Z, CC_TMP0, C_irr, 10, 31, 12, 31, 0, "Z (m>n,ei)");
    dpd_buf4_init(&Z1, CC_TMP0, C_irr, 31, 30, 31, 30, 0, "Z1 (ei,mb)");
    dpd_contract424(&Z, &tia, &Z1, 1, 0, 0, 1, 0);
    dpd_buf4_close(&Z);

    dpd_buf4_init(&Ht, CC_TMP0, C_irr, 15, 31, 17, 31, 0, "Ht_Wabei (a>b,ei)");
    dpd_contract244(&tia, &Z1, &Ht, 0, 2, 0, 1, 1);
    dpd_buf4_close(&Ht);
    dpd_buf4_close(&Z1);


    /** W(Ab,Ei) <---  <Nm|Ef> * t1[m][b]*C[i][f]*t1[N][A] **/

    dpd_buf4_init(&Z, CC_TMP0, C_irr, 22, 26, 22, 26, 0, "Z (Nm,Ei)");
    dpd_buf4_init(&D, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    dpd_contract424(&D, &Cme, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&D);

    dpd_buf4_init(&Z1, CC_TMP0, C_irr, 26, 24, 26, 24, 0, "Z1 (Ei,Nb)");
    dpd_contract424(&Z, &tia, &Z1, 1, 0, 0, 1, 0);
    dpd_buf4_close(&Z);

    dpd_buf4_init(&Ht, CC_TMP0, C_irr, 28, 26, 28, 26, 0, "Ht_WAbEi (Ab,Ei)");
    dpd_contract244(&tIA, &Z1, &Ht, 0, 2, 0, 1, 1);
    dpd_buf4_close(&Ht);

    /** W(aB,eI) <---  <nM|eF> * t1[M][B]*C[I][F]*t1[n][a] **/

    dpd_buf4_init(&Z, CC_TMP0, C_irr, 23, 25, 23, 25, 0, "Z (nM,eI)");
    dpd_buf4_init(&D, CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");
    dpd_contract424(&D, &CME, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&D);

    dpd_buf4_init(&Z1, CC_TMP0, C_irr, 25, 27, 25, 27, 0, "Z1 (eI,nB)");
    dpd_contract424(&Z, &tIA, &Z1, 1, 0, 0, 1, 0);
    dpd_buf4_close(&Z);

    dpd_buf4_init(&Ht, CC_TMP0, C_irr, 29, 25, 29, 25, 0, "Ht_WaBeI (aB,eI)");
    dpd_contract244(&tia, &Z1, &Ht, 0, 2, 0, 1, 1);
    dpd_buf4_close(&Ht);


    dpd_file2_close(&CME);
    dpd_file2_close(&Cme);
    dpd_file2_close(&tIA);
    dpd_file2_close(&tia);


    /* combine W(Ab,Ei) and W(Ei,Ab)  */

    dpd_buf4_init(&Ht, CC_TMP0, C_irr, 28, 26, 28, 26, 0, "Ht_WAbEi (Ab,Ei)");
    dpd_buf4_sort_axpy(&Ht, CC_TMP0, rspq, 26, 28, "Ht_WAbEi (Ei,Ab)", 1); 
    dpd_buf4_close(&Ht);

    dpd_buf4_init(&Ht, CC_TMP0, C_irr, 29, 25, 29, 25, 0, "Ht_WaBeI (aB,eI)");
    dpd_buf4_sort_axpy(&Ht, CC_TMP0, rspq, 25, 29, "Ht_WaBeI (eI,aB)", 1); 
    dpd_buf4_close(&Ht);

    /* sort to Wabei (ei,ab) */
    dpd_buf4_init(&W, CC_TMP0, C_irr, 7, 21, 7, 21, 0, "Ht_WABEI (A>B,EI)");
    dpd_buf4_sort(&W, CC_TMP2, rspq, 21, 7, "Ht_WABEI (EI,A>B)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC_TMP0, C_irr, 17, 31, 17, 31, 0, "Ht_Wabei (a>b,ei)");
    dpd_buf4_sort(&W, CC_TMP2, rspq, 31, 17, "Ht_Wabei (ei,a>b)");
    dpd_buf4_close(&W);

    /* sort to Wabei (ie,ba) */
    dpd_buf4_init(&W, CC_TMP2, C_irr, 21, 7, 21, 7, 0, "Ht_WABEI (EI,A>B)");
    dpd_buf4_sort(&W, CC3_HC1ET1, qprs, 20, 7, "Ht_WABEI (IE,B>A)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC3_HC1ET1, C_irr, 20, 7, 20, 7, 0, "Ht_WABEI (IE,B>A)");
    dpd_buf4_scm(&W, -1.0);
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP2, C_irr, 31, 17, 31, 17, 0, "Ht_Wabei (ei,a>b)");
    dpd_buf4_sort(&W, CC3_HC1ET1, qprs, 30, 17, "Ht_Wabei (ie,b>a)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC3_HC1ET1, C_irr, 30, 17, 30, 17, 0, "Ht_Wabei (ie,b>a)");
    dpd_buf4_scm(&W, -1.0);
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP0, C_irr, 26, 28, 26, 28, 0, "Ht_WAbEi (Ei,Ab)");
    dpd_buf4_sort(&W, CC3_HC1ET1, qpsr, 27, 29, "Ht_WAbEi (iE,bA)");
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP0, C_irr, 25, 29, 25, 29, 0, "Ht_WaBeI (eI,aB)");
    dpd_buf4_sort(&W, CC3_HC1ET1, qpsr, 24, 28, "Ht_WaBeI (Ie,Ba)");
    dpd_buf4_close(&W);


    /************ TEST *************/
    /*
    dpd_buf4_init(&W, CC3_HC1ET1, 0, 20, 5, 20, 7, 0, "Ht_WABEI (IE,B>A)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    fprintf(outfile,"\t<WABEI(IE,B>A)|WABEI> = %15.10lf\n", dot);

    dpd_buf4_init(&W, CC3_HC1ET1, 0, 30, 15, 30, 17, 0, "Ht_Wabei (ie,b>a)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    fprintf(outfile,"\t<Wabei (ie,b>a)|Wabei> = %15.10lf\n", dot);

    dpd_buf4_init(&W, CC3_HC1ET1, 0, 27, 29, 27, 29, 0, "Ht_WAbEi (iE,bA)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    fprintf(outfile,"\t<WAbEi (iE,bA)|WAbEi> = %15.10lf\n", dot);

    dpd_buf4_init(&W, CC3_HC1ET1, 0, 24, 28, 24, 28, 0, "Ht_WaBeI (Ie,Ba)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    fprintf(outfile,"\t<WaBeI (Ie,Ba)|WaBeI> = %15.10lf\n", dot);
    */

  }

  return;
}

void HC1ET1_Wmbij_rhf(int i, int C_irr)
{
  double dot;
  dpdbuf4 C, D, E, F, Ht, W, W1, X, Z;
  dpdfile2 CME;
  char CME_lbl[32];
  sprintf(CME_lbl, "%s %d", "CME", i);

  dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);

  /***** Ht (Mb,Ij) <--- -WMnIj * Cnb *****/
  dpd_buf4_init(&Ht, CC3_HC1ET1, C_irr, 10, 0, 10, 0, 0, "Ht_WMbIj (Mb,Ij)");
  dpd_buf4_init(&W, CC3_HET1, 0, 0, 0, 0, 0, 0, "CC3 WMnIj (Mn,Ij)");
  dpd_contract424(&W, &CME, &Ht, 1, 0, 1, -1.0, 0.0);
  dpd_buf4_close(&W);
  dpd_buf4_close(&Ht);
  
  /***** Ht (Mb,Ij) <--- CIE * WMbEj *****/ 
  dpd_buf4_init(&Z, CC_TMP0, C_irr, 0, 10, 0, 10, 0, "Z (MI,jb)");
  dpd_buf4_init(&W, CC3_HET1, 0, 10, 10, 10, 10, 0, "CC3 WMbEj (ME,jb)");
  dpd_contract424(&W, &CME, &Z, 1, 1, 1, 1.0, 0.0);
  dpd_buf4_close(&W);
  dpd_buf4_sort_axpy(&Z, CC3_HC1ET1, psqr, 10, 0, "Ht_WMbIj (Mb,Ij)", 1);
  dpd_buf4_close(&Z);
    
  /***** Ht (Mb,Ij) <--- Cje * WMbeI *****/
  dpd_buf4_init(&Z, CC_TMP0, C_irr, 0, 10, 0, 10, 0, "Z (Mj,Ib)");
  dpd_buf4_init(&W, CC3_HET1, 0, 10, 10, 10, 10, 0, "CC3 WMbeJ (Me,Jb)");
  dpd_contract424(&W, &CME, &Z, 1, 1, 1, 1.0, 0.0);
  dpd_buf4_close(&W);
  dpd_buf4_sort_axpy(&Z, CC3_HC1ET1, psrq, 10, 0, "Ht_WMbIj (Mb,Ij)", -1);
  dpd_buf4_close(&Z);

  dpd_buf4_init(&Ht, CC3_HC1ET1, C_irr, 10, 0, 10, 0, 0, "Ht_WMbIj (Mb,Ij)");
  dpd_buf4_sort(&Ht, CC3_HC1ET1, rspq, 0, 10, "Ht_WMbIj (Ij,Mb)");
  dpd_buf4_close(&Ht);
  
  dpd_file2_close(&CME);

  return;
}

void HC1ET1_Wabei_rhf(int i, int C_irr)
{
  double dot;
  dpdfile2 CME, tIA;
  dpdbuf4 Ht, Z, Z1, Z2, Z3, B, C, D, E, F, W, X;
  char CME_lbl[32];
  sprintf(CME_lbl, "%s %d", "CME", i);

  dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);
  dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");

  /***** Ht_WAbEi <--- -CMA * WMbEi *****/
  dpd_buf4_init(&Ht, CC_TMP0, C_irr, 5, 11, 5, 11, 0, "Ht_WAbEi (Ab,Ei)");
  dpd_buf4_init(&W, CC3_HET1, 0, 10, 11, 10, 11, 0, "CC3 WMbEj (Mb,Ej)");
  dpd_contract244(&CME, &W, &Ht, 0, 0, 0, -1.0, 0.0);
  dpd_buf4_close(&W);
  dpd_buf4_close(&Ht);

  /***** Ht_WAbEi <--- WAmEi * Cmb *****/
  dpd_buf4_init(&Ht, CC_TMP0, C_irr, 11, 5, 11, 5, 0, "Ht_WAbEi (Ei,Ab)");
  dpd_buf4_init(&W, CC3_HET1, 0, 11, 11, 11, 11, 0, "CC3 WMbeJ (bM,eJ)");
  dpd_contract424(&W, &CME, &Ht, 1, 0, 0, 1.0, 0.0);
  dpd_buf4_close(&W);
  dpd_buf4_close(&Ht);

  /***** HAbEi <--- <Ab|Ef> * Cif *****/
  dpd_buf4_init(&Ht, CC_TMP0, C_irr, 5, 11, 5, 11, 0, "Ht_WAbEi (Ab,Ei)");
  dpd_buf4_init(&B, CC_BINTS, 0, 5, 5, 5, 5, 0, "B <ab|cd>");
  dpd_contract424(&B, &CME, &Ht, 3, 1, 0, 1.0, 1.0);
  dpd_buf4_close(&B);
  dpd_buf4_close(&Ht);

  /** H(Ab,Ei) <--- -<Am|Ef>*C[i][f]*t1[m][b] **/
  dpd_buf4_init(&Z, CC_TMP0, C_irr, 11, 11, 11, 11, 0, "Z(Am,Ei)");
  dpd_buf4_init(&F, CC_FINTS, 0, 11, 5, 11, 5, 0, "F <ai|bc>");
  dpd_contract424(&F, &CME, &Z, 3, 1, 0, 1, 0);
  dpd_buf4_close(&F);

  dpd_buf4_init(&Ht, CC_TMP0, C_irr, 11, 5, 11, 5, 0, "Ht_WAbEi (Ei,Ab)");
  dpd_contract424(&Z, &tIA, &Ht, 1, 0, 0, -1, 1);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&Ht);

  /** Ht(Ab,Ei) <--- -<Mb|Ef>*C[i][f]*t1[M][A] **/
  dpd_buf4_init(&Z, CC_TMP0, C_irr, 10, 11, 10, 11, 0, "Z(Mb,Ei)");
  dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
  dpd_contract424(&F, &CME, &Z, 3, 1, 0, 1, 0);
  dpd_buf4_close(&F);

  dpd_buf4_init(&Ht, CC_TMP0, C_irr, 5, 11, 5, 11, 0, "Ht_WAbEi (Ab,Ei)");
  dpd_contract244(&tIA, &Z, &Ht, 0, 0, 0, -1, 1);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&Ht);

  /** W(Ab,Ei) <---  <Nm|Ef> * t1[m][b]*C[i][f]*t1[N][A] **/
  dpd_buf4_init(&Z, CC_TMP0, C_irr, 0, 11, 0, 11, 0, "Z (Nm,Ei)");
  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_contract424(&D, &CME, &Z, 3, 1, 0, 1, 0);
  dpd_buf4_close(&D);

  dpd_buf4_init(&Z1, CC_TMP0, C_irr, 11, 10, 11, 10, 0, "Z1 (Ei,Nb)");
  dpd_contract424(&Z, &tIA, &Z1, 1, 0, 0, 1, 0);
  dpd_buf4_close(&Z);

  dpd_buf4_init(&Ht, CC_TMP0, C_irr, 5, 11, 5, 11, 0, "Ht_WAbEi (Ab,Ei)");
  dpd_contract244(&tIA, &Z1, &Ht, 0, 2, 0, 1, 1);
  dpd_buf4_close(&Ht);

  dpd_file2_close(&CME);
  dpd_file2_close(&tIA);

  /* combine W(Ab,Ei) and W(Ei,Ab)  */
  dpd_buf4_init(&Ht, CC_TMP0, C_irr, 5, 11, 5, 11, 0, "Ht_WAbEi (Ab,Ei)");
  dpd_buf4_sort_axpy(&Ht, CC_TMP0, rspq, 11, 5, "Ht_WAbEi (Ei,Ab)", 1);
  dpd_buf4_close(&Ht);

  /* sort to Wabei (ie,ba) */
  dpd_buf4_init(&W, CC_TMP0, C_irr, 11, 5, 11, 5, 0, "Ht_WAbEi (Ei,Ab)");
  dpd_buf4_sort(&W, CC3_HC1ET1, qpsr, 10, 5, "Ht_WAbEi (iE,bA)");
  dpd_buf4_close(&W);

  return;
}


}} // namespace psi::cceom
