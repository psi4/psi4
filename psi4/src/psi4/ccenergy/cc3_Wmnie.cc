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
#include "Params.h"
#include "MOInfo.h"
#include "ccwave.h"

namespace psi { namespace ccenergy {

/* cc3_Wmnie(): Compute the Wmnie matrix from CC3 theory, which is
** given in spin-orbitals as:
**
** Wmnie = <mn||ie> + t_i^f <mn||fe>
**
** TDC, Feb 2004
*/

void purge_Wmnie(void);

void CCEnergyWavefunction::cc3_Wmnie(void)
{
  dpdbuf4 E, D, W, Z;
  dpdfile2 tIA, tia;

  if(params_.ref == 0) { /** RHF **/

    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    global_dpd_->buf4_copy(&E, PSIF_CC3_HET1, "CC3 WMnIe (Mn,Ie)");
    global_dpd_->buf4_close(&E);

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 0, 10, 0, 10, 0, "CC3 WMnIe (Mn,Ie)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract244(&tIA, &D, &W, 1, 2, 1, 1, 1);
    global_dpd_->file2_close(&tIA);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&W);
  }

  else if (params_.ref == 1) { /* ROHF */

    /** W(M>N,IE) <--- <MN||IE> **/
    /** W(m>n,ie) <--- <mn||ie> **/
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 2, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
    global_dpd_->buf4_sort(&E, PSIF_CC3_HET1, pqsr, 2, 11, "CC3 WMNIE (M>N,EI)");
    global_dpd_->buf4_sort(&E, PSIF_CC3_HET1, pqsr, 2, 11, "CC3 Wmnie (m>n,ei)");
    global_dpd_->buf4_close(&E);

    /** W(Mn,Ie) <--- <Mn|Ie> **/
    /** W(mN,iE) <--- <mN|iE> **/
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    global_dpd_->buf4_sort(&E, PSIF_CC3_HET1, pqsr, 0, 11, "CC3 WMnIe (Mn,eI)");
    global_dpd_->buf4_sort(&E, PSIF_CC3_HET1, pqsr, 0, 11, "CC3 WmNiE (mN,Ei)");
    global_dpd_->buf4_close(&E);

    /**** Term 2 ****/

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");

    /* <M>N||EF> T(I,F) --> W(M>N,EI) */
    /* <m>n||ef> T(i,f) --> W(m>n,ei) */
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 2, 11, 2, 11, 0, "CC3 WMNIE (M>N,EI)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 2, 5, 2, 5, 0, "D <ij||ab> (i>j,ab)");
    global_dpd_->contract424(&D, &tIA, &W, 3, 1, 0, -1, 1.0);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 2, 11, 2, 11, 0, "CC3 Wmnie (m>n,ei)");
    global_dpd_->contract424(&D, &tia, &W, 3, 1, 0, -1, 1.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&W);

    /* Z(nM,eI) = <nM|eF> T(I,F) */
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, 0, 0, 11, 0, 11, 0, "Z(nM,eI)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->contract424(&D, &tIA, &Z, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&D);
    /* Z(nM,eI) --> W(Mn,eI) */
    global_dpd_->buf4_sort_axpy(&Z, PSIF_CC3_HET1, qprs, 0, 11, "CC3 WMnIe (Mn,eI)", 1);
    global_dpd_->buf4_close(&Z);

    /* Z(Nm,Ei) = <Nm|Ef> T(i,f) */
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, 0, 0, 11, 0, 11, 0, "Z(Nm,Ei)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->contract424(&D, &tia, &Z, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&D);
    /* Z(Nm,Ei) --> W(mN,Ei) */
    global_dpd_->buf4_sort_axpy(&Z, PSIF_CC3_HET1, qprs, 0, 11, "CC3 WmNiE (mN,Ei)", 1);
    global_dpd_->buf4_close(&Z);

    /* purge (mn,ei)'s before sorting */
    purge_Wmnie();

    /* also put "normal" sorted versions in CC3_HET1 */
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 2, 11, 2, 11, 0, "CC3 WMNIE (M>N,EI)");
    global_dpd_->buf4_sort(&W, PSIF_CC3_HET1, pqsr, 2, 10, "CC3 WMNIE (M>N,IE)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 2, 11, 2, 11, 0, "CC3 Wmnie (m>n,ei)");
    global_dpd_->buf4_sort(&W, PSIF_CC3_HET1, pqsr, 2, 10, "CC3 Wmnie (m>n,ie)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 0, 11, 0, 11, 0, "CC3 WMnIe (Mn,eI)");
    global_dpd_->buf4_sort(&W, PSIF_CC3_HET1, pqsr, 0, 10, "CC3 WMnIe (Mn,Ie)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 0, 11, 0, 11, 0, "CC3 WmNiE (mN,Ei)");
    global_dpd_->buf4_sort(&W, PSIF_CC3_HET1, pqsr, 0, 10, "CC3 WmNiE (mN,iE)");
    global_dpd_->buf4_close(&W);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);

  }

  else if (params_.ref == 2) {

    /** W(M>N,IE) <--- <MN||IE> **/
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 2, 20, 2, 20, 0, "E <IJ||KA> (I>J,KA)");
    global_dpd_->buf4_sort(&E, PSIF_CC3_HET1, pqsr, 2, 21, "CC3 WMNIE (M>N,EI)");
    global_dpd_->buf4_close(&E);

    /** W(m>n,ie) <--- <mn||ie> **/
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 12, 30, 12, 30, 0, "E <ij||ka> (i>j,ka)");
    global_dpd_->buf4_sort(&E, PSIF_CC3_HET1, pqsr, 12, 31, "CC3 Wmnie (m>n,ei)");
    global_dpd_->buf4_close(&E);

    /** W(Mn,Ie) <--- <Mn|Ie> **/
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
    global_dpd_->buf4_sort(&E, PSIF_CC3_HET1, pqsr, 22, 25, "CC3 WMnIe (Mn,eI)");
    global_dpd_->buf4_close(&E);

    /** W(mN,iE) <--- <mN|iE> **/
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 23, 27, 23, 27, 0, "E <iJ|kA>");
    global_dpd_->buf4_sort(&E, PSIF_CC3_HET1, pqsr, 23, 26, "CC3 WmNiE (mN,Ei)");
    global_dpd_->buf4_close(&E);

    /**** Term 2 ****/

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");

    /* <M>N||EF> T(I,F) --> W(M>N,EI) */
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 2, 21, 2, 21, 0, "CC3 WMNIE (M>N,EI)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 2, 5, 2, 5, 0, "D <IJ||AB> (I>J,AB)");
    global_dpd_->contract424(&D, &tIA, &W, 3, 1, 0, -1, 1.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&W);

    /* <m>n||ef> T(i,f) --> W(m>n,ei) */
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 12, 31, 12, 31, 0, "CC3 Wmnie (m>n,ei)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 12, 15, 12, 15, 0, "D <ij||ab> (i>j,ab)");
    global_dpd_->contract424(&D, &tia, &W, 3, 1, 0, -1, 1.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&W);

    /* Z(nM,eI) = <nM|eF> T(I,F) */
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, 0, 23, 25, 23, 25, 0, "Z(nM,eI)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");
    global_dpd_->contract424(&D, &tIA, &Z, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&D);
    /* Z(nM,eI) --> W(Mn,eI) */
    global_dpd_->buf4_sort_axpy(&Z, PSIF_CC3_HET1, qprs, 22, 25, "CC3 WMnIe (Mn,eI)", 1);
    global_dpd_->buf4_close(&Z);

    /* Z(Nm,Ei) = <Nm|Ef> T(i,f) */
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, 0, 22, 26, 22, 26, 0, "Z(Nm,Ei)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    global_dpd_->contract424(&D, &tia, &Z, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&D);
    /* Z(Nm,Ei) --> W(mN,Ei) */
    global_dpd_->buf4_sort_axpy(&Z, PSIF_CC3_HET1, qprs, 23, 26, "CC3 WmNiE (mN,Ei)", 1);
    global_dpd_->buf4_close(&Z);

    /* also put "normal" sorted versions in CC3_HET1 */
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 2, 21, 2, 21, 0, "CC3 WMNIE (M>N,EI)");
    global_dpd_->buf4_sort(&W, PSIF_CC3_HET1, pqsr, 2, 20, "CC3 WMNIE (M>N,IE)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 12, 31, 12, 31, 0, "CC3 Wmnie (m>n,ei)");
    global_dpd_->buf4_sort(&W, PSIF_CC3_HET1, pqsr, 12, 30, "CC3 Wmnie (m>n,ie)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 22, 25, 22, 25, 0, "CC3 WMnIe (Mn,eI)");
    global_dpd_->buf4_sort(&W, PSIF_CC3_HET1, pqsr, 22, 24, "CC3 WMnIe (Mn,Ie)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 23, 26, 23, 26, 0, "CC3 WmNiE (mN,Ei)");
    global_dpd_->buf4_sort(&W, PSIF_CC3_HET1, pqsr, 23, 27, "CC3 WmNiE (mN,iE)");
    global_dpd_->buf4_close(&W);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);

  }
}


/* Purge Wmnie matrix elements */
void CCEnergyWavefunction::purge_Wmnie(void) {
  dpdfile4 W;
  int *occpi, *virtpi;
  int h, a, b, e, f, i, j, m, n;
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

  global_dpd_->file4_init(&W, PSIF_CC3_HET1, 0, 0, 11,"CC3 WMnIe (Mn,eI)");
  for(h=0; h < nirreps; h++) {
    global_dpd_->file4_mat_irrep_init(&W, h);
    global_dpd_->file4_mat_irrep_rd(&W, h);
    for(mn=0; mn<W.params->rowtot[h]; mn++) {
      n = W.params->roworb[h][mn][1];
      nsym = W.params->qsym[n];
      N = n - occ_off[nsym];
      for(ei=0; ei<W.params->coltot[h]; ei++) {
        if (N >= (occpi[nsym] - openpi[nsym]))
          W.matrix[h][mn][ei] = 0.0;
      }
    }
    global_dpd_->file4_mat_irrep_wrt(&W, h);
    global_dpd_->file4_mat_irrep_close(&W, h);
  }

  global_dpd_->file4_init(&W, PSIF_CC3_HET1, 0, 2, 11, "CC3 WMNIE (M>N,EI)");
  for(h=0; h < W.params->nirreps; h++) {
    global_dpd_->file4_mat_irrep_init(&W, h);
    global_dpd_->file4_mat_irrep_rd(&W, h);
    for(mn=0; mn<W.params->rowtot[h]; mn++) {
      for(ei=0; ei<W.params->coltot[h]; ei++) {
        e = W.params->colorb[h][ei][0];
        esym = W.params->rsym[e];
        E = e - vir_off[esym];
        if (E >= (virtpi[esym] - openpi[esym]))
          W.matrix[h][mn][ei] = 0.0;
      }
    }
    global_dpd_->file4_mat_irrep_wrt(&W, h);
    global_dpd_->file4_mat_irrep_close(&W, h);
  }
  global_dpd_->file4_close(&W);
  global_dpd_->file4_init(&W, PSIF_CC3_HET1, 0, 2, 11,"CC3 Wmnie (m>n,ei)");
  for(h=0; h < nirreps; h++) {
    global_dpd_->file4_mat_irrep_init(&W, h);
    global_dpd_->file4_mat_irrep_rd(&W, h);
    for(mn=0; mn<W.params->rowtot[h]; mn++) {
      m = W.params->roworb[h][mn][0];
      n = W.params->roworb[h][mn][1];
      msym = W.params->psym[m];
      nsym = W.params->qsym[n];
      M = m - occ_off[msym];
      N = n - occ_off[nsym];
      for(ei=0; ei<W.params->coltot[h]; ei++) {
        i = W.params->colorb[h][ei][1];
        isym = W.params->ssym[i];
        I = i - occ_off[isym];
        if ((M >= (occpi[msym] - openpi[msym])) ||
          (N >= (occpi[nsym] - openpi[nsym])) ||
          (I >= (occpi[isym] - openpi[isym])) )
          W.matrix[h][mn][ei] = 0.0;
      }
    }
    global_dpd_->file4_mat_irrep_wrt(&W, h);
    global_dpd_->file4_mat_irrep_close(&W, h);
  }
  global_dpd_->file4_close(&W);

  global_dpd_->file4_init(&W, PSIF_CC3_HET1, 0, 0, 11,"CC3 WmNiE (mN,Ei)");
  for(h=0; h < nirreps; h++) {
    global_dpd_->file4_mat_irrep_init(&W, h);
    global_dpd_->file4_mat_irrep_rd(&W, h);
    for(mn=0; mn<W.params->rowtot[h]; mn++) {
      m = W.params->roworb[h][mn][0];
      msym = W.params->psym[m];
      M = m - occ_off[msym];
      for(ei=0; ei<W.params->coltot[h]; ei++) {
        e = W.params->colorb[h][ei][0];
        i = W.params->colorb[h][ei][1];
        esym = W.params->rsym[e];
        isym = W.params->ssym[i];
        E = e - vir_off[esym];
        I = i - occ_off[isym];
        if ((M >= (occpi[msym] - openpi[msym])) ||
            (E >= (virtpi[esym] - openpi[esym])) ||
            (I >= (occpi[isym] - openpi[isym])) )
          W.matrix[h][mn][ei] = 0.0;
      }
    }
    global_dpd_->file4_mat_irrep_wrt(&W, h);
    global_dpd_->file4_mat_irrep_close(&W, h);
  }
  global_dpd_->file4_close(&W);
  return;
}

}} // namespace psi::ccenergy
