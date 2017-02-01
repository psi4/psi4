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

/* cc3_Wamef(): Compute the Wamef matrix from CC3 theory, which is
** given in spin-orbitals as:
**
** Wamef = <am||ef> - t_n^a <nm||ef>
**
** TDC, Feb 2004
*/

void purge_Wamef(void);

void CCEnergyWavefunction::cc3_Wamef(void)
{
  dpdbuf4 F, D, W;
  dpdfile2 t1,tia,tIA;

  if(params_.ref == 0) { /** RHF **/

    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    global_dpd_->buf4_sort(&F, PSIF_CC3_HET1, qpsr, 11, 5, "CC3 WAmEf (Am,Ef)");
    global_dpd_->buf4_close(&F);

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 11, 5, 11, 5, 0, "CC3 WAmEf (Am,Ef)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->file2_init(&t1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract244(&t1, &D, &W, 0, 0, 0, -1, 1);
    global_dpd_->file2_close(&t1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_sort(&W, PSIF_CC3_HET1, qprs, 10, 5, "CC3 WAmEf (mA,Ef)");
    global_dpd_->buf4_close(&W);
  }

  else if (params_.ref == 1) { /** ROHF **/

    /** W(AM,E>F) <--- <AM||EF> **/
    /** W(am,e>f) <--- <am||ef> **/
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 11, 7, 11, 5, 1, "F <ai|bc>");
    global_dpd_->buf4_copy(&F, PSIF_CC3_HET1, "CC3 WAMEF (AM,E>F)");
    global_dpd_->buf4_copy(&F, PSIF_CC3_HET1, "CC3 Wamef (am,e>f)");
    global_dpd_->buf4_close(&F);

    /** W(Am,Ef) <--- <Am|Ef> **/
    /** W(aM,eF) <--- <aM|eF> **/
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 11, 5, 11, 5, 0, "F <ai|bc>");
    global_dpd_->buf4_copy(&F, PSIF_CC3_HET1, "CC3 WAmEf (Am,Ef)");
    global_dpd_->buf4_copy(&F, PSIF_CC3_HET1, "CC3 WaMeF (aM,eF)");
    global_dpd_->buf4_close(&F);

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");

    /* t(N,A) <NM||EF> --> W(AM,E>F) */
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 11, 7, 11, 7, 0, "CC3 WAMEF (AM,E>F)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 7, 0, 7, 0, "D <ij||ab> (ij,a>b)");
    global_dpd_->contract244(&tIA, &D, &W, 0, 0, 0, -1, 1.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_sort(&W, PSIF_CC3_HET1, qprs, 10, 7, "CC3 WAMEF (MA,F>E)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 10, 7, 10, 7, 0, "CC3 WAMEF (MA,F>E)");
    global_dpd_->buf4_scm(&W, -1.0);
    global_dpd_->buf4_close(&W);

    /* t(n,a) <nm||ef> --> W(am,e>f) */
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 11, 7, 11, 7, 0, "CC3 Wamef (am,e>f)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 7, 0, 7, 0, "D <ij||ab> (ij,a>b)");
    global_dpd_->contract244(&tia, &D, &W, 0, 0, 0, -1, 1.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_sort(&W, PSIF_CC3_HET1, qprs, 10, 7, "CC3 Wamef (ma,f>e)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 10, 7, 10, 7, 0, "CC3 Wamef (ma,f>e)");
    global_dpd_->buf4_scm(&W, -1.0);
    global_dpd_->buf4_close(&W);

    /* t(N,A) <Nm|Ef> --> W(Am,Ef) */
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 11, 5, 11, 5, 0, "CC3 WAmEf (Am,Ef)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->contract244(&tIA, &D, &W, 0, 0, 0, -1, 1.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_sort(&W, PSIF_CC3_HET1, qpsr, 10, 5, "CC3 WAmEf (mA,fE)");
    global_dpd_->buf4_close(&W);

    /* t(n,a) <nM|eF> --> W(aM,eF) */
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 11, 5, 11, 5, 0, "CC3 WaMeF (aM,eF)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->contract244(&tia, &D, &W, 0, 0, 0, -1, 1.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_sort(&W, PSIF_CC3_HET1, qpsr, 10, 5, "CC3 WaMeF (Ma,Fe)");
    global_dpd_->buf4_close(&W);

    global_dpd_->file2_close(&tia);
    global_dpd_->file2_close(&tIA);

    purge_Wamef();
  }

  else if (params_.ref == 2) {

    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");

    /** W(AM,E>F) <--- <AM||EF> **/
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 21, 7, 21, 5, 1, "F <AI|BC>");
    global_dpd_->buf4_copy(&F, PSIF_CC3_HET1, "CC3 WAMEF (AM,E>F)");
    global_dpd_->buf4_close(&F);

    /** W(am,e>f) <--- <am||ef> **/
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 31, 17, 31, 15, 1, "F <ai|bc>");
    global_dpd_->buf4_copy(&F, PSIF_CC3_HET1, "CC3 Wamef (am,e>f)");
    global_dpd_->buf4_close(&F);

    /** W(Am,Ef) <--- <Am|Ef> **/
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 26, 28, 26, 28, 0, "F <Ai|Bc>");
    global_dpd_->buf4_copy(&F, PSIF_CC3_HET1, "CC3 WAmEf (Am,Ef)");
    global_dpd_->buf4_close(&F);

    /** W(aM,eF) <--- <aM|eF> **/
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 25, 29, 25, 29, 0, "F <aI|bC>");
    global_dpd_->buf4_copy(&F, PSIF_CC3_HET1, "CC3 WaMeF (aM,eF)");
    global_dpd_->buf4_close(&F);

    /** W(AM,E>F) <--- tNA * <NM||EF> **/
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 21, 7, 21, 7, 0, "CC3 WAMEF (AM,E>F)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 7, 0, 7, 0, "D <IJ||AB> (IJ,A>B)");
    global_dpd_->contract244(&tIA, &D, &W, 0, 0, 0, -1, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_sort(&W, PSIF_CC3_HET1, qprs, 20, 7, "CC3 WAMEF (MA,F>E)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 20, 7, 20, 7, 0, "CC3 WAMEF (MA,F>E)");
    global_dpd_->buf4_scm(&W, -1.0);
    global_dpd_->buf4_close(&W);

    /** W(am,e>f) <--- tna * <nm||ef> **/
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 31, 17, 31, 17, 0, "CC3 Wamef (am,e>f)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 17, 10, 17, 0, "D <ij||ab> (ij,a>b)");
    global_dpd_->contract244(&tia, &D, &W, 0, 0, 0, -1, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_sort(&W, PSIF_CC3_HET1, qprs, 30, 17, "CC3 Wamef (ma,f>e)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 30, 17, 30, 17, 0, "CC3 Wamef (ma,f>e)");
    global_dpd_->buf4_scm(&W, -1.0);
    global_dpd_->buf4_close(&W);

    /** W(Am,Ef) <--- tNA * <Nm|Ef> **/
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 26, 28, 26, 28, 0, "CC3 WAmEf (Am,Ef)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    global_dpd_->contract244(&tIA, &D, &W, 0, 0, 0, -1, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_sort(&W, PSIF_CC3_HET1, qpsr, 27, 29, "CC3 WAmEf (mA,fE)");
    global_dpd_->buf4_close(&W);

    /** W(aM,eF) <--- tna * <nM|eF> **/
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 25, 29, 25, 29, 0, "CC3 WaMeF (aM,eF)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");
    global_dpd_->contract244(&tia, &D, &W, 0, 0, 0, -1, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_sort(&W, PSIF_CC3_HET1, qpsr, 24, 28, "CC3 WaMeF (Ma,Fe)");
    global_dpd_->buf4_close(&W);

    global_dpd_->file2_close(&tia);
    global_dpd_->file2_close(&tIA);
  }
}

void CCEnergyWavefunction::purge_Wamef(void) {
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

  /* Purge Wamef matrix elements */
  global_dpd_->file4_init(&W, PSIF_CC3_HET1, 0, 11, 7,"CC3 WAMEF (AM,E>F)");
  for(h=0; h < nirreps; h++) {
    global_dpd_->file4_mat_irrep_init(&W, h);
    global_dpd_->file4_mat_irrep_rd(&W, h);
    for(ma=0; ma < W.params->rowtot[h]; ma++) {
      a = W.params->roworb[h][ma][0];
      asym = W.params->psym[a];
      A = a - vir_off[asym];
      for(ef=0; ef< W.params->coltot[h]; ef++) {
        e = W.params->colorb[h][ef][0];
        f = W.params->colorb[h][ef][1];
        esym = W.params->rsym[e];
        fsym = W.params->ssym[f];
        E = e - vir_off[esym];
        F = f - vir_off[fsym];
        if ((A >= (virtpi[asym] - openpi[asym])) ||
            (E >= (virtpi[esym] - openpi[esym])) ||
            (F >= (virtpi[fsym] - openpi[fsym])) )
          W.matrix[h][ma][ef] = 0.0;
      }
    }
    global_dpd_->file4_mat_irrep_wrt(&W, h);
    global_dpd_->file4_mat_irrep_close(&W, h);
  }
  global_dpd_->file4_close(&W);

  global_dpd_->file4_init(&W, PSIF_CC3_HET1, 0, 11, 7,"CC3 Wamef (am,e>f)");
  for(h=0; h < nirreps; h++) {
    global_dpd_->file4_mat_irrep_init(&W, h);
    global_dpd_->file4_mat_irrep_rd(&W, h);
    for(ma=0; ma < W.params->rowtot[h]; ma++) {
      m = W.params->roworb[h][ma][1];
      msym = W.params->qsym[m];
      M = m - occ_off[msym];
      for(ef=0; ef< W.params->coltot[h]; ef++) {
        if (M >=  (occpi[msym] - openpi[msym]))
          W.matrix[h][ma][ef] = 0.0;
      }
    }
    global_dpd_->file4_mat_irrep_wrt(&W, h);
    global_dpd_->file4_mat_irrep_close(&W, h);
  }
  global_dpd_->file4_close(&W);

  global_dpd_->file4_init(&W, PSIF_CC3_HET1, 0, 11, 5,"CC3 WAmEf (Am,Ef)");
  for(h=0; h < nirreps; h++) {
    global_dpd_->file4_mat_irrep_init(&W, h);
    global_dpd_->file4_mat_irrep_rd(&W, h);
    for(ma=0; ma < W.params->rowtot[h]; ma++) {
      a = W.params->roworb[h][ma][0];
      m = W.params->roworb[h][ma][1];
      asym = W.params->psym[a];
      msym = W.params->qsym[m];
      M = m - occ_off[msym];
      A = a - vir_off[asym];
      for(ef=0; ef< W.params->coltot[h]; ef++) {
        e = W.params->colorb[h][ef][0];
        esym = W.params->rsym[e];
        E = e - vir_off[esym];
        if ((A >= (virtpi[asym] - openpi[asym])) ||
            (M >=  (occpi[msym] - openpi[msym])) ||
            (E >= (virtpi[esym] - openpi[esym])) )
          W.matrix[h][ma][ef] = 0.0;
      }
    }
    global_dpd_->file4_mat_irrep_wrt(&W, h);
    global_dpd_->file4_mat_irrep_close(&W, h);
  }
  global_dpd_->file4_close(&W);

  global_dpd_->file4_init(&W, PSIF_CC3_HET1, 0, 11, 5,"CC3 WaMeF (aM,eF)");
  for(h=0; h < nirreps; h++) {
    global_dpd_->file4_mat_irrep_init(&W, h);
    global_dpd_->file4_mat_irrep_rd(&W, h);
    for(ma=0; ma < W.params->rowtot[h]; ma++) {
      for(ef=0; ef< W.params->coltot[h]; ef++) {
        f = W.params->colorb[h][ef][1];
        fsym = W.params->ssym[f];
        F = f - vir_off[fsym];
        if (F >= (virtpi[fsym] - openpi[fsym]))
          W.matrix[h][ma][ef] = 0.0;
      }
    }
    global_dpd_->file4_mat_irrep_wrt(&W, h);
    global_dpd_->file4_mat_irrep_close(&W, h);
  }
  global_dpd_->file4_close(&W);

  return;
}
}} // namespace psi::ccenergy
