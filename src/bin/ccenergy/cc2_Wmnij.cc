/*! \file
    \ingroup CCENERGY
    \brief Enter brief description of file here 
*/
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include "Params.h"
#define EXTERN
#include "MOInfo.h"
#include "globals.h"

namespace psi { namespace ccenergy {

/* cc2_Wmnij(): Compute the Wmnij matrix from CC2 theory, which is
** given in spin-orbitals as:
**
** Wmnij = <mn||ij> + P(ij) t_j^e <mn||ie> + t_i^e t_j^f <mn||ef>
**
** TDC, Feb 2004
*/

void purge_cc2_Wmnij(void);

void cc2_Wmnij_build(void)
{
  dpdbuf4 A, E, D, Z, W, Z1, X;
  dpdfile2 t1, tIA, tia;

  timer_on("A->Wmnij");
  if(params.ref == 0) { /** RHF **/
    /* Wmnij <- <mn||ij> */
    dpd_buf4_init(&A, CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
    dpd_buf4_copy(&A, CC2_HET1, "CC2 WMnIj");
    dpd_buf4_close(&A);
  }
  else if (params.ref == 1) { /** ROHF **/
    /** W(M>N,I>J) <--- <MN||IJ> **/
    /** W(m>n,i>j) <--- <mn||ij> **/
    dpd_buf4_init(&A, CC_AINTS, 0, 2, 2, 0, 0, 1, "A <ij|kl>");
    dpd_buf4_copy(&A, CC2_HET1, "CC2 WMNIJ (M>N,I>J)");
    dpd_buf4_copy(&A, CC2_HET1, "CC2 Wmnij (m>n,i>j)");
    dpd_buf4_close(&A);

    /** W(Mn,Ij) <--- <Mn|Ij> **/
    dpd_buf4_init(&A, CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
    dpd_buf4_copy(&A, CC2_HET1, "CC2 WMnIj");
    dpd_buf4_close(&A);
  }
  else if (params.ref == 2) { /** UHF **/
    /** W(M>N,I>J) <--- <MN||IJ> **/
    dpd_buf4_init(&A, CC_AINTS, 0, 2, 2, 0, 0, 1, "A <IJ|KL>");
    dpd_buf4_copy(&A, CC2_HET1, "CC2 WMNIJ (M>N,I>J)");
    dpd_buf4_close(&A);

    /** W(m>n,i>j) <--- <mn||ij> **/
    dpd_buf4_init(&A, CC_AINTS, 0, 12, 12, 10, 10, 1, "A <ij|kl>");
    dpd_buf4_copy(&A, CC2_HET1, "CC2 Wmnij (m>n,i>j)");
    dpd_buf4_close(&A);

    /** W(Mn,Ij) <--- <Mn|Ij> **/
    dpd_buf4_init(&A, CC_AINTS, 0, 22, 22, 22, 22, 0, "A <Ij|Kl>");
    dpd_buf4_copy(&A, CC2_HET1, "CC2 WMnIj");
    dpd_buf4_close(&A);
  }
  timer_off("A->Wmnij");

  timer_on("E->Wmnij");
  if(params.ref == 0) { /** RHF **/

    dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");

    /* Wmnij <- + P(ij) t(j,e) * <mn||ie> */
    dpd_buf4_init(&Z, CC_TMP0, 0, 0, 0, 0, 0, 0, "CC2 ZMnIj");
    dpd_buf4_init(&E, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    dpd_contract424(&E, &t1, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&E);

    dpd_buf4_init(&W, CC2_HET1, 0, 0, 0, 0, 0, 0, "CC2 WMnIj");
    dpd_buf4_axpy(&Z, &W, 1);
    dpd_buf4_close(&W);
    dpd_buf4_sort_axpy(&Z, CC2_HET1, qpsr, 0, 0, "CC2 WMnIj", 1);
    dpd_buf4_close(&Z);

    dpd_file2_close(&t1);
  }
  else if (params.ref == 1) { /** ROHF **/

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");

    /**** W(M>N,I>J) <-- ZMNIJ <-- P(I/J)( <MN||IE> * t1[J][E] ) ****/
    dpd_buf4_init(&Z, CC_TMP0, 0, 2, 0, 2, 0, 0, "Z (M>N,IJ)");
    dpd_buf4_init(&E, CC_EINTS, 0, 2, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
    dpd_contract424(&E, &tIA, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&E);

    dpd_buf4_sort(&Z, CC_TMP0, pqsr, 2, 0, "Z (M>N,JI)");
    dpd_buf4_init(&Z1, CC_TMP0, 0, 2, 0, 2, 0, 0, "Z (M>N,JI)");
    dpd_buf4_axpy(&Z1, &Z, -1);
    dpd_buf4_close(&Z1);

    dpd_buf4_init(&W, CC2_HET1, 0, 2, 0, 2, 2, 0, "CC2 WMNIJ (M>N,I>J)");
    dpd_buf4_axpy(&Z, &W, 1);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&W);

    /**** W(m>n,i>j) <-- Zmnij <-- P(i/j)( <mn||ie> * t1[j][e] ) ****/
    dpd_buf4_init(&Z, CC_TMP0, 0, 2, 0, 2, 0, 0, "Z (m>n,ij)");
    dpd_buf4_init(&E, CC_EINTS, 0, 2, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
    dpd_contract424(&E, &tia, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&E);

    dpd_buf4_sort(&Z, CC_TMP0, pqsr, 2, 0, "Z (m>n,ji)");
    dpd_buf4_init(&Z1, CC_TMP0, 0, 2, 0, 2, 0, 0, "Z (m>n,ji)");
    dpd_buf4_axpy(&Z1, &Z, -1);
    dpd_buf4_close(&Z1);

    dpd_buf4_init(&W, CC2_HET1, 0, 2, 0, 2, 2, 0, "CC2 Wmnij (m>n,i>j)");
    dpd_buf4_axpy(&Z, &W, 1);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&W);

    /**** W(Mn,Ij) <-- <Mn|Ie> * t1[j][e] ****/
    dpd_buf4_init(&W, CC2_HET1, 0, 0, 0, 0, 0, 0, "CC2 WMnIj");
    dpd_buf4_init(&E, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    dpd_contract424(&E, &tia, &W, 3, 1, 0, 1, 1);
    dpd_buf4_close(&E);
    dpd_buf4_close(&W);

    /**** W(Mn,Ij) <-- <Mn|Ej> * t1[I][E] ****/
    dpd_buf4_init(&W, CC2_HET1, 0, 0, 0, 0, 0, 0, "CC2 WMnIj");
    dpd_buf4_init(&E, CC_EINTS, 0, 0, 11, 0, 11, 0, "E <ij|ak>");
    dpd_contract244(&tIA, &E, &W, 1, 2, 1, 1, 1);
    dpd_buf4_close(&E);
    dpd_buf4_close(&W);

    dpd_file2_close(&tIA);
    dpd_file2_close(&tia);
  }
  else if (params.ref == 2) { /** UHF **/

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 2, 3, "tia");

    /**** W(M>N,I>J) <-- ZMNIJ <-- P(I/J)( <MN||IE> * t1[J][E] ) ****/
    dpd_buf4_init(&Z, CC_TMP0, 0, 2, 0, 2, 0, 0, "Z (M>N,IJ)");
    dpd_buf4_init(&E, CC_EINTS, 0, 2, 20, 2, 20, 0, "E <IJ||KA> (I>J,KA)");
    dpd_contract424(&E, &tIA, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&E);

    dpd_buf4_sort(&Z, CC_TMP0, pqsr, 2, 0, "Z (M>N,JI)");
    dpd_buf4_init(&Z1, CC_TMP0, 0, 2, 0, 2, 0, 0, "Z (M>N,JI)");
    dpd_buf4_axpy(&Z1, &Z, -1);
    dpd_buf4_close(&Z1);

    dpd_buf4_init(&W, CC2_HET1, 0, 2, 0, 2, 2, 0, "CC2 WMNIJ (M>N,I>J)");
    dpd_buf4_axpy(&Z, &W, 1);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&W);

    /**** W(m>n,i>j) <-- Zmnij <-- P(i/j)( <mn||ie> * t1[j][e] ) ****/
    dpd_buf4_init(&Z, CC_TMP0, 0, 12, 10, 12, 10, 0, "Z (m>n,ij)");
    dpd_buf4_init(&E, CC_EINTS, 0, 12, 30, 12, 30, 0, "E <ij||ka> (i>j,ka)");
    dpd_contract424(&E, &tia, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&E);

    dpd_buf4_sort(&Z, CC_TMP0, pqsr, 12, 10, "Z (m>n,ji)");
    dpd_buf4_init(&Z1, CC_TMP0, 0, 12, 10, 12, 10, 0, "Z (m>n,ji)");
    dpd_buf4_axpy(&Z1, &Z, -1);
    dpd_buf4_close(&Z1);

    dpd_buf4_init(&W, CC2_HET1, 0, 12, 10, 12, 12, 0, "CC2 Wmnij (m>n,i>j)");
    dpd_buf4_axpy(&Z, &W, 1);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&W);

    /**** W(Mn,Ij) <-- <Mn|Ie> * t1[j][e] ****/
    dpd_buf4_init(&W, CC2_HET1, 0, 22, 22, 22, 22, 0, "CC2 WMnIj");
    dpd_buf4_init(&E, CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
    dpd_contract424(&E, &tia, &W, 3, 1, 0, 1, 1);
    dpd_buf4_close(&E);
    dpd_buf4_close(&W);

    /**** W(Mn,Ij) <-- <Mn|Ej> * t1[I][E] ****/
    dpd_buf4_init(&W, CC2_HET1, 0, 22, 22, 22, 22, 0, "CC2 WMnIj");
    dpd_buf4_init(&E, CC_EINTS, 0, 22, 26, 22, 26, 0, "E <Ij|Ak>");
    dpd_contract244(&tIA, &E, &W, 1, 2, 1, 1, 1);
    dpd_buf4_close(&E);
    dpd_buf4_close(&W);

    dpd_file2_close(&tIA);
    dpd_file2_close(&tia);
  }
  timer_off("E->Wmnij");

  timer_on("D->Wmnij");
  if(params.ref == 0) { /** RHF **/

    dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");

    /* Wmnij<- +1/2 P(ij) t(i,e) t(j,f) * <mn||ef> */
    dpd_buf4_init(&Z, CC_TMP0, 0, 0, 10, 0, 10, 0, "CC2 ZMnIf (Mn,If)");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_contract244(&t1, &D, &Z, 1, 2, 1, 1, 0);
    dpd_buf4_close(&D);

    dpd_buf4_init(&Z1, CC_TMP0, 0, 0, 0, 0, 0, 0, "CC2 ZMnIj");
    dpd_contract424(&Z, &t1, &Z1, 3, 1, 0, 0.5, 0);
    dpd_buf4_close(&Z);
    dpd_buf4_init(&W, CC2_HET1, 0, 0, 0, 0, 0, 0, "CC2 WMnIj");
    dpd_buf4_axpy(&Z1, &W, 1);
    dpd_buf4_close(&W);
    dpd_buf4_sort_axpy(&Z1, CC2_HET1, qpsr, 0, 0, "CC2 WMnIj", 1);
    dpd_buf4_close(&Z1);

    dpd_file2_close(&t1);
  }
  else if (params.ref == 1) { /** ROHF **/

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");

    /**** W(M>N,I>J) <-- tIE tJF <MN||EF> ****/
    dpd_buf4_init(&Z, CC_TMP0, 0, 2, 11, 2, 11, 0, "Z (M>N,EJ)");
    dpd_buf4_init(&D, CC_DINTS, 0, 2, 5, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
    dpd_contract424(&D, &tIA, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&D);

    dpd_buf4_init(&W, CC2_HET1, 0, 2, 0, 2, 2, 0, "CC2 WMNIJ (M>N,I>J)");
    dpd_contract244(&tIA, &Z, &W, 1, 2, 1, 1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Z);

    /**** W(M>N,I>J) <-- tIE tJF <MN||EF> ****/
    dpd_buf4_init(&Z, CC_TMP0, 0, 2, 11, 2, 11, 0, "Z (m>n,ej)");
    dpd_buf4_init(&D, CC_DINTS, 0, 2, 5, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
    dpd_contract424(&D, &tia, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&D);

    dpd_buf4_init(&W, CC2_HET1, 0, 2, 0, 2, 2, 0, "CC2 Wmnij (m>n,i>j)");
    dpd_contract244(&tia, &Z, &W, 1, 2, 1, 1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Z);

    /**** W(Mn,Ij) <-- tIE tjf <Mn|Ef> ****/
    dpd_buf4_init(&Z, CC_TMP0, 0, 0, 11, 0, 11, 0, "Z (Mn,Ej)");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_contract424(&D, &tia, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&D);

    dpd_buf4_init(&W, CC2_HET1, 0, 0, 0, 0, 0, 0, "CC2 WMnIj");
    dpd_contract244(&tIA, &Z, &W, 1, 2, 1, 1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Z);

    dpd_file2_close(&tIA);
    dpd_file2_close(&tia);

    purge_cc2_Wmnij();
  }
  else if (params.ref == 2) { /** UHF **/

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 2, 3, "tia");

    /**** W(M>N,I>J) <-- tIE tJF <MN||EF> ****/
    dpd_buf4_init(&Z, CC_TMP0, 0, 2, 21, 2, 21, 0, "Z (M>N,EJ)");
    dpd_buf4_init(&D, CC_DINTS, 0, 2, 5, 2, 7, 0, "D <IJ||AB> (I>J,A>B)");
    dpd_contract424(&D, &tIA, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&D);

    dpd_buf4_init(&W, CC2_HET1, 0, 2, 0, 2, 2, 0, "CC2 WMNIJ (M>N,I>J)");
    dpd_contract244(&tIA, &Z, &W, 1, 2, 1, 1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Z);

    /**** W(M>N,I>J) <-- tIE tJF <MN||EF> ****/
    dpd_buf4_init(&Z, CC_TMP0, 0, 12, 31, 12, 31, 0, "Z (m>n,ej)");
    dpd_buf4_init(&D, CC_DINTS, 0, 12, 15, 12, 17, 0, "D <ij||ab> (i>j,a>b)");
    dpd_contract424(&D, &tia, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&D);

    dpd_buf4_init(&W, CC2_HET1, 0, 12, 10, 12, 12, 0, "CC2 Wmnij (m>n,i>j)");
    dpd_contract244(&tia, &Z, &W, 1, 2, 1, 1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Z);

    /**** W(Mn,Ij) <-- tIE tjf <Mn|Ef> ****/
    dpd_buf4_init(&Z, CC_TMP0, 0, 22, 26, 22, 26, 0, "Z (Mn,Ej)");
    dpd_buf4_init(&D, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    dpd_contract424(&D, &tia, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&D);

    dpd_buf4_init(&W, CC2_HET1, 0, 22, 22, 22, 22, 0, "CC2 WMnIj");
    dpd_contract244(&tIA, &Z, &W, 1, 2, 1, 1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Z);

    dpd_file2_close(&tIA);
    dpd_file2_close(&tia);
  }
  timer_off("D->Wmnij");

}


void purge_cc2_Wmnij(void) {
  dpdfile2 FAE, Fmi, FME, Fme;
  dpdfile4 W;
  int *occpi, *virtpi;
  int h, a, b, e, f, i, j, m, n, omit;
  int    A, B, E, F, I, J, M, N;
  int mn, ei, ma, ef, me, jb, mb, ij, ab;
  int asym, bsym, esym, fsym, isym, jsym, msym, nsym;
  int *occ_off, *vir_off;
  int *occ_sym, *vir_sym;
  int *openpi, nirreps;

  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi; virtpi = moinfo.virtpi;
  occ_off = moinfo.occ_off; vir_off = moinfo.vir_off;
  occ_sym = moinfo.occ_sym; vir_sym = moinfo.vir_sym;
  openpi = moinfo.openpi;

  /* Purge Wmnij matrix elements */
  dpd_file4_init(&W, CC2_HET1, 0, 2, 2,"CC2 Wmnij (m>n,i>j)");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
    for(mn=0; mn < W.params->rowtot[h]; mn++) {
      m = W.params->roworb[h][mn][0];
      n = W.params->roworb[h][mn][1];
      msym = W.params->psym[m];
      nsym = W.params->qsym[n];
      M = m - occ_off[msym];
      N = n - occ_off[nsym];
      for(ij=0; ij < W.params->coltot[h]; ij++) {
        i = W.params->colorb[h][ij][0];
        j = W.params->colorb[h][ij][1];
        isym = W.params->rsym[i];
        jsym = W.params->ssym[j];
        I = i - occ_off[isym];
        J = j - occ_off[jsym];
        if ((I >= (occpi[isym] - openpi[isym])) ||
            (J >= (occpi[jsym] - openpi[jsym])) ||
            (M >= (occpi[msym] - openpi[msym])) ||
            (N >= (occpi[nsym] - openpi[nsym])) )
          W.matrix[h][mn][ij] = 0.0;
      }
    }
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);

  dpd_file4_init(&W, CC2_HET1, 0, 0, 0,"CC2 WMnIj");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
    for(mn=0; mn < W.params->rowtot[h]; mn++) {
      n = W.params->roworb[h][mn][1];
      nsym = W.params->qsym[n];
      N = n - occ_off[nsym];
      for(ij=0; ij < W.params->coltot[h]; ij++) {
        j = W.params->colorb[h][ij][1];
        jsym = W.params->ssym[j];
        J = j - occ_off[jsym];
        if ((J >= (occpi[jsym] - openpi[jsym])) ||
            (N >= (occpi[nsym] - openpi[nsym])) )
          W.matrix[h][mn][ij] = 0.0;
      }
    }
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);
}
}} // namespace psi::ccenergy
