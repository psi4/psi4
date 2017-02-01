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
    \ingroup CCENERGY
    \brief Enter brief description of file here
*/
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#include "ccwave.h"

namespace psi { namespace ccenergy {

/* cc3_Wmnij(): Compute the Wmnij components of the
** T1-similarity-transformed Hamiltonian matrix, which is given in
** spin-orbitals as:
**
** Wmnij = <mn||ij> + P(ij) t_j^e <mn||ie> + t_i^e t_j^f <mn||ef>
**
** TDC, Feb 2004
*/

void purge_Wmnij(void);

void CCEnergyWavefunction::cc3_Wmnij(void)
{
  dpdbuf4 A, E, D, Z, W, Z1, X;
  dpdfile2 t1, tIA, tia;

  if(params_.ref == 0) { /** RHF **/

    global_dpd_->buf4_init(&A, PSIF_CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
    global_dpd_->buf4_copy(&A, PSIF_CC3_HET1, "CC3 WMnIj (Mn,Ij)");
    global_dpd_->buf4_close(&A);

    global_dpd_->file2_init(&t1, PSIF_CC_OEI, 0, 0, 1, "tIA");

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 0, 0, 0, 0, 0, "CC3 ZMnIj (Mn,Ij)");
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    global_dpd_->contract424(&E, &t1, &Z, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&E);
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 0, 0, 0, 0, 0, "CC3 WMnIj (Mn,Ij)");
    global_dpd_->buf4_axpy(&Z, &W, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_sort_axpy(&Z, PSIF_CC3_HET1, qpsr, 0, 0, "CC3 WMnIj (Mn,Ij)", 1);
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 0, 10, 0, 10, 0, "CC3 ZMnIf (Mn,If)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->contract244(&t1, &D, &Z, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 0, 0, 0, 0, 0, "CC3 WMnIj (Mn,Ij)");
    global_dpd_->contract424(&Z, &t1, &W, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Z);

    global_dpd_->file2_close(&t1);
  }

  else if (params_.ref == 1) {
    /** W(M>N,I>J) <--- <MN||IJ> **/
    /** W(m>n,i>j) <--- <mn||ij> **/
    global_dpd_->buf4_init(&A, PSIF_CC_AINTS, 0, 2, 2, 0, 0, 1, "A <ij|kl>");
    global_dpd_->buf4_copy(&A, PSIF_CC3_HET1, "CC3 WMNIJ (M>N,I>J)");
    global_dpd_->buf4_copy(&A, PSIF_CC3_HET1, "CC3 Wmnij (m>n,i>j)");
    global_dpd_->buf4_close(&A);

    /** W(Mn,Ij) <--- <Mn|Ij> **/
    global_dpd_->buf4_init(&A, PSIF_CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
    global_dpd_->buf4_copy(&A, PSIF_CC3_HET1, "CC3 WMnIj (Mn,Ij)");
    global_dpd_->buf4_close(&A);

    /** term 2 **/

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");

    /**** W(M>N,I>J) <-- ZMNIJ <-- P(I/J)( <MN||IE> * t1[J][E] ) ****/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 2, 0, 2, 0, 0, "Z (M>N,IJ)");
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 2, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
    global_dpd_->contract424(&E, &tIA, &Z, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&E);

    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, pqsr, 2, 0, "Z (M>N,JI)");
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 2, 0, 2, 0, 0, "Z (M>N,JI)");
    global_dpd_->buf4_axpy(&Z1, &Z, -1);
    global_dpd_->buf4_close(&Z1);

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 2, 0, 2, 2, 0, "CC3 WMNIJ (M>N,I>J)");
    global_dpd_->buf4_axpy(&Z, &W, 1);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&W);

    /**** W(m>n,i>j) <-- Zmnij <-- P(i/j)( <mn||ie> * t1[j][e] ) ****/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 2, 0, 2, 0, 0, "Z (m>n,ij)");
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 2, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
    global_dpd_->contract424(&E, &tia, &Z, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&E);

    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, pqsr, 2, 0, "Z (m>n,ji)");
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 2, 0, 2, 0, 0, "Z (m>n,ji)");
    global_dpd_->buf4_axpy(&Z1, &Z, -1);
    global_dpd_->buf4_close(&Z1);

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 2, 0, 2, 2, 0, "CC3 Wmnij (m>n,i>j)");
    global_dpd_->buf4_axpy(&Z, &W, 1);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&W);

    /**** W(Mn,Ij) <-- <Mn|Ie> * t1[j][e] ****/
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 0, 0, 0, 0, 0, "CC3 WMnIj (Mn,Ij)");
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    global_dpd_->contract424(&E, &tia, &W, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&E);
    global_dpd_->buf4_close(&W);

    /**** W(Mn,Ij) <-- <Mn|Ej> * t1[I][E] ****/
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 0, 0, 0, 0, 0, "CC3 WMnIj (Mn,Ij)");
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 0, 11, 0, 11, 0, "E <ij|ak>");
    global_dpd_->contract244(&tIA, &E, &W, 1, 2, 1, 1.0, 1);
    global_dpd_->buf4_close(&E);
    global_dpd_->buf4_close(&W);

    /** term 3 **/

    /**** W(M>N,I>J) <-- tIE tJF <MN||EF> ****/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 2, 11, 2, 11, 0, "Z (M>N,EJ)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 2, 5, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
    global_dpd_->contract424(&D, &tIA, &Z, 3, 1, 0, 1.0, 0);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 2, 0, 2, 2, 0, "CC3 WMNIJ (M>N,I>J)");
    global_dpd_->contract244(&tIA, &Z, &W, 1, 2, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Z);

    /**** W(M>N,I>J) <-- tIE tJF <MN||EF> ****/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 2, 11, 2, 11, 0, "Z (m>n,ej)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 2, 5, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
    global_dpd_->contract424(&D, &tia, &Z, 3, 1, 0, 1.0, 0);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 2, 0, 2, 2, 0, "CC3 Wmnij (m>n,i>j)");
    global_dpd_->contract244(&tia, &Z, &W, 1, 2, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Z);

    /**** W(Mn,Ij) <-- tIE tjf <Mn|Ef> ****/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 0, 11, 0, 11, 0, "Z (Mn,Ej)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->contract424(&D, &tia, &Z, 3, 1, 0, 1.0, 0);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 0, 0, 0, 0, 0, "CC3 WMnIj (Mn,Ij)");
    global_dpd_->contract244(&tIA, &Z, &W, 1, 2, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Z);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);

    purge_Wmnij();
  }

  else if (params_.ref == 2) {

    /** W(M>N,I>J) <--- <MN||IJ> **/
    global_dpd_->buf4_init(&A, PSIF_CC_AINTS, 0, 2, 2, 0, 0, 1, "A <IJ|KL>");
    global_dpd_->buf4_copy(&A, PSIF_CC3_HET1, "CC3 WMNIJ (M>N,I>J)");
    global_dpd_->buf4_close(&A);

    /** W(m>n,i>j) <--- <mn||ij> **/
    global_dpd_->buf4_init(&A, PSIF_CC_AINTS, 0, 12, 12, 10, 10, 1, "A <ij|kl>");
    global_dpd_->buf4_copy(&A, PSIF_CC3_HET1, "CC3 Wmnij (m>n,i>j)");
    global_dpd_->buf4_close(&A);

    /** W(Mn,Ij) <--- <Mn|Ij> **/
    global_dpd_->buf4_init(&A, PSIF_CC_AINTS, 0, 22, 22, 22, 22, 0, "A <Ij|Kl>");
    global_dpd_->buf4_copy(&A, PSIF_CC3_HET1, "CC3 WMnIj (Mn,Ij)");
    global_dpd_->buf4_close(&A);

    /** term 2 **/

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");

    /**** W(M>N,I>J) <-- ZMNIJ <-- P(I/J)( <MN||IE> * t1[J][E] ) ****/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 2, 0, 2, 0, 0, "Z (M>N,IJ)");
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 2, 20, 2, 20, 0, "E <IJ||KA> (I>J,KA)");
    global_dpd_->contract424(&E, &tIA, &Z, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&E);

    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, pqsr, 2, 0, "Z (M>N,JI)");
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 2, 0, 2, 0, 0, "Z (M>N,JI)");
    global_dpd_->buf4_axpy(&Z1, &Z, -1);
    global_dpd_->buf4_close(&Z1);

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 2, 0, 2, 2, 0, "CC3 WMNIJ (M>N,I>J)");
    global_dpd_->buf4_axpy(&Z, &W, 1);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&W);

    /**** W(m>n,i>j) <-- Zmnij <-- P(i/j)( <mn||ie> * t1[j][e] ) ****/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 12, 10, 12, 10, 0, "Z (m>n,ij)");
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 12, 30, 12, 30, 0, "E <ij||ka> (i>j,ka)");
    global_dpd_->contract424(&E, &tia, &Z, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&E);

    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, pqsr, 12, 10, "Z (m>n,ji)");
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 12, 10, 12, 10, 0, "Z (m>n,ji)");
    global_dpd_->buf4_axpy(&Z1, &Z, -1);
    global_dpd_->buf4_close(&Z1);

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 12, 10, 12, 12, 0, "CC3 Wmnij (m>n,i>j)");
    global_dpd_->buf4_axpy(&Z, &W, 1);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&W);

    /**** W(Mn,Ij) <-- <Mn|Ie> * t1[j][e] ****/
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 22, 22, 22, 22, 0, "CC3 WMnIj (Mn,Ij)");
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
    global_dpd_->contract424(&E, &tia, &W, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&E);
    global_dpd_->buf4_close(&W);

    /**** W(Mn,Ij) <-- <Mn|Ej> * t1[I][E] ****/
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 22, 22, 22, 22, 0, "CC3 WMnIj (Mn,Ij)");
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 22, 26, 22, 26, 0, "E <Ij|Ak>");
    global_dpd_->contract244(&tIA, &E, &W, 1, 2, 1, 1.0, 1);
    global_dpd_->buf4_close(&E);
    global_dpd_->buf4_close(&W);

    /** term 3 **/

    /**** W(M>N,I>J) <-- tIE tJF <MN||EF> ****/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 2, 21, 2, 21, 0, "Z (M>N,EJ)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 2, 5, 2, 7, 0, "D <IJ||AB> (I>J,A>B)");
    global_dpd_->contract424(&D, &tIA, &Z, 3, 1, 0, 1.0, 0);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 2, 0, 2, 2, 0, "CC3 WMNIJ (M>N,I>J)");
    global_dpd_->contract244(&tIA, &Z, &W, 1, 2, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Z);

    /**** W(M>N,I>J) <-- tIE tJF <MN||EF> ****/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 12, 31, 12, 31, 0, "Z (m>n,ej)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 12, 15, 12, 17, 0, "D <ij||ab> (i>j,a>b)");
    global_dpd_->contract424(&D, &tia, &Z, 3, 1, 0, 1.0, 0);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 12, 10, 12, 12, 0, "CC3 Wmnij (m>n,i>j)");
    global_dpd_->contract244(&tia, &Z, &W, 1, 2, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Z);

    /**** W(Mn,Ij) <-- tIE tjf <Mn|Ef> ****/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 22, 26, 22, 26, 0, "Z (Mn,Ej)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    global_dpd_->contract424(&D, &tia, &Z, 3, 1, 0, 1.0, 0);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 22, 22, 22, 22, 0, "CC3 WMnIj (Mn,Ij)");
    global_dpd_->contract244(&tIA, &Z, &W, 1, 2, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Z);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);
  }
}


void CCEnergyWavefunction::purge_Wmnij(void) {
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

  nirreps = moinfo_.nirreps;
  occpi = moinfo_.occpi; virtpi = moinfo_.virtpi;
  occ_off = moinfo_.occ_off; vir_off = moinfo_.vir_off;
  occ_sym = moinfo_.occ_sym; vir_sym = moinfo_.vir_sym;
  openpi = moinfo_.openpi;

  /* Purge Wmnij matrix elements */
  global_dpd_->file4_init(&W, PSIF_CC3_HET1, 0, 2, 2,"CC3 Wmnij (m>n,i>j)");
  for(h=0; h < nirreps; h++) {
    global_dpd_->file4_mat_irrep_init(&W, h);
    global_dpd_->file4_mat_irrep_rd(&W, h);
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
    global_dpd_->file4_mat_irrep_wrt(&W, h);
    global_dpd_->file4_mat_irrep_close(&W, h);
  }
  global_dpd_->file4_close(&W);

  global_dpd_->file4_init(&W, PSIF_CC3_HET1, 0, 0, 0,"CC3 WMnIj (Mn,Ij)");
  for(h=0; h < nirreps; h++) {
    global_dpd_->file4_mat_irrep_init(&W, h);
    global_dpd_->file4_mat_irrep_rd(&W, h);
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
    global_dpd_->file4_mat_irrep_wrt(&W, h);
    global_dpd_->file4_mat_irrep_close(&W, h);
  }
  global_dpd_->file4_close(&W);
}
}} // namespace psi::ccenergy
