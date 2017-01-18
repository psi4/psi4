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
    \ingroup CCEOM
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cceom {

/* This function computes the extra contributions to sigma_1 and sigma_2
  for EOM_CC3 computations that are not normally present in a EOM_CCSD
  calculation */

/* The additional terms are:
 * Term 1:
 * <S| H    <T| (Uhat C2)c   |0> |T> / (w-wt) -> sigma_1
 * <D| Hhat <T| (Uhat C2)c   |0> |T> / (w-wt) -> sigma_2
 * Term 2:
 * <S| H    <T| (Utilde T2)c |0> |T> / (w-wt) -> sigma_1
 * <D| Hhat <T| (Utilde T2)c |0> |T> / (w-wt) -> sigma_2
 * Term 3:
 * <D| H'   <T| (Uhat T2)c   |0> |T> / (-wt) -> sigma_2
 *
 *  See Eqn. (83) in JCP, 103, 7429, 1995
 *  All three terms can be evaluated by the same functions in
 *  cc3_sigma_RHF given different matrix elements.
 *
 * - RHF case modifed from ccenergy cc3 code, RAK 2006
 * */

/* extern void cc3_sigma_RHF(dpdbuf4 *CIjAb, dpdbuf4 *WAbEi, dpdbuf4 *WMbIj,
    int do_singles, dpdbuf4 *Dints, dpdfile2 *SIA,
    int do_doubles, dpdfile2 *FME, dpdbuf4 *WAmEf, dpdbuf4 *WMnIe,
    dpdbuf4 *SIjAb, double energy); */

void sigmaCC3_RHF_obsolete(int i, int C_irr, double omega)
{
  int ii,j,a,b,A,B,Ga,Gb,Gij=0,ab,ab0,ab1;
  dpdfile2 SIA, FME;
  dpdbuf4 CMnEf, WAbEi, WMbIj, Dints, WmAEf, WMnIe, SIjAb;
  dpdbuf4 tIjAb;
  char lbl[32];

  sprintf(lbl, "%s %d", "SIA", i);
  global_dpd_->file2_init(&SIA, PSIF_EOM_SIA, C_irr, 0, 1, lbl);
  sprintf(lbl, "%s %d", "SIjAb", i);
  global_dpd_->buf4_init(&SIjAb, PSIF_EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, lbl);

  /*** alpha-alpha-beta term 1 ***/
  /* quantities to compute X3 */
  sprintf(lbl, "%s %d", "CMnEf", i);
  global_dpd_->buf4_init(&CMnEf, PSIF_EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, lbl);
  global_dpd_->buf4_init(&WAbEi, PSIF_CC3_HET1, 0, 10, 5, 10, 5, 0, "CC3 WAbEi (iE,bA)");
  global_dpd_->buf4_init(&WMbIj, PSIF_CC3_HET1, 0, 0, 10, 0, 10, 0, "CC3 WMbIj (Ij,Mb)");
  /* quantities to compute sigma */
  global_dpd_->buf4_init(&Dints, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
  global_dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");
  global_dpd_->buf4_init(&WmAEf, PSIF_CC3_HET1, 0, 10, 5, 10, 5, 0, "CC3 WAmEf (mA,Ef)");
  global_dpd_->buf4_init(&WMnIe, PSIF_CC3_HET1, 0, 0, 10, 0, 10, 0, "CC3 WMnIe (Mn,Ie)");

       /* * <S| H    <T| (Uhat C2)c   |0> |T> / (w-wt) -> sigma_1
          * <D| Hhat <T| (Uhat C2)c   |0> |T> / (w-wt) -> sigma_2 */

  global_dpd_->cc3_sigma_RHF(&CMnEf, &WAbEi, &WMbIj, 1,  &Dints, &SIA,
    1, &FME, &WmAEf, &WMnIe, &SIjAb, moinfo.occpi, moinfo.occ_off,
    moinfo.virtpi, moinfo.vir_off, omega, "outfile", params.newtrips);

  global_dpd_->buf4_close(&CMnEf);
  global_dpd_->buf4_close(&WAbEi);
  global_dpd_->buf4_close(&WMbIj);
  global_dpd_->buf4_close(&Dints);
  global_dpd_->file2_close(&FME);
  global_dpd_->buf4_close(&WmAEf);
  global_dpd_->buf4_close(&WMnIe);

#ifdef EOM_DEBUG
  dpd_file2_close(&SIA);
  dpd_buf4_close(&SIjAb);
  check_sum("<Psi|Hhat<T|(Uhat C2)c|0>|T>/(w-wt)", i, C_irr);
  sprintf(lbl, "%s %d", "SIA", i);
  dpd_file2_init(&SIA, EOM_SIA, C_irr, 0, 1, lbl);
  sprintf(lbl, "%s %d", "SIjAb", i);
  dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, lbl);
#endif

  /* do alpha-alpha-beta term 2 */
  global_dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  global_dpd_->buf4_init(&WAbEi, PSIF_CC3_HC1ET1, C_irr, 10, 5, 10, 5, 0, "Ht_WAbEi (iE,bA)");
  global_dpd_->buf4_init(&WMbIj, PSIF_CC3_HC1ET1, C_irr, 0, 10, 0, 10, 0, "Ht_WMbIj (Ij,Mb)");

  global_dpd_->buf4_init(&Dints, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
  global_dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");
  global_dpd_->buf4_init(&WmAEf, PSIF_CC3_HET1, 0, 10, 5, 10, 5, 0, "CC3 WAmEf (mA,Ef)");
  global_dpd_->buf4_init(&WMnIe, PSIF_CC3_HET1, 0, 0, 10, 0, 10, 0, "CC3 WMnIe (Mn,Ie)");

         /* * <S| H    <T| (Utilde T2)c |0> |T> / (w-wt) -> sigma_1
            * <D| Hhat <T| (Utilde T2)c |0> |T> / (w-wt) -> sigma_2 */

  global_dpd_->cc3_sigma_RHF(&tIjAb, &WAbEi, &WMbIj, 1,  &Dints, &SIA,
     1, &FME, &WmAEf, &WMnIe, &SIjAb, moinfo.occpi, moinfo.occ_off,
     moinfo.virtpi, moinfo.vir_off, omega, "outfile", params.newtrips);

  global_dpd_->buf4_close(&tIjAb);
  global_dpd_->buf4_close(&WAbEi);
  global_dpd_->buf4_close(&WMbIj);
  global_dpd_->buf4_close(&Dints);
  global_dpd_->file2_close(&FME);
  global_dpd_->buf4_close(&WmAEf);
  global_dpd_->buf4_close(&WMnIe);

#ifdef EOM_DEBUG
  dpd_file2_close(&SIA);
  dpd_buf4_close(&SIjAb);
  check_sum("<Psi|Hhat<T|(Utilde T2)c|0>|T>/(w-wt)", i, C_irr);
  sprintf(lbl, "%s %d", "SIA", i);
  dpd_file2_init(&SIA, EOM_SIA, C_irr, 0, 1, lbl);
  sprintf(lbl, "%s %d", "SIjAb", i);
  dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, lbl);
#endif

  /* alpha-alpha-beta term 3 */
  global_dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  global_dpd_->buf4_init(&WAbEi, PSIF_CC3_HET1, 0, 10, 5, 10, 5, 0, "CC3 WAbEi (iE,bA)");
  global_dpd_->buf4_init(&WMbIj, PSIF_CC3_HET1, 0, 0, 10, 0, 10, 0, "CC3 WMbIj (Ij,Mb)");

  global_dpd_->file2_init(&FME, PSIF_CC3_HC1, C_irr, 0, 1, "HC1 FME");
  global_dpd_->buf4_init(&WmAEf, PSIF_CC3_HC1, C_irr, 10, 5, 10, 5, 0, "HC1 WAmEf (mA,Ef)");
  global_dpd_->buf4_init(&WMnIe, PSIF_CC3_HC1, C_irr, 0, 10, 0, 10, 0, "HC1 WMnIe (Mn,Ie)");

         /* <D| H'   <T| (Uhat T2)c   |0> |T> / (-wt) -> sigma_2 */

  global_dpd_->cc3_sigma_RHF(&tIjAb, &WAbEi, &WMbIj, 0, NULL, NULL,
     1, &FME, &WmAEf, &WMnIe, &SIjAb, moinfo.occpi, moinfo.occ_off,
     moinfo.virtpi, moinfo.vir_off, 0.0, "outfile", params.newtrips);

  global_dpd_->buf4_close(&tIjAb);
  global_dpd_->buf4_close(&WAbEi);
  global_dpd_->buf4_close(&WMbIj);
  global_dpd_->file2_close(&FME);
  global_dpd_->buf4_close(&WmAEf);
  global_dpd_->buf4_close(&WMnIe);

#ifdef EOM_DEBUG
  dpd_file2_close(&SIA);
  dpd_buf4_close(&SIjAb);
  check_sum("<Psi|H'<T|(Uhat T2)c|0>|T>/(w-wt)", i, C_irr);
  sprintf(lbl, "%s %d", "SIA", i);
  dpd_file2_init(&SIA, EOM_SIA, C_irr, 0, 1, lbl);
  sprintf(lbl, "%s %d", "SIjAb", i);
  dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, lbl);
#endif

  global_dpd_->file2_close(&SIA);
  global_dpd_->buf4_close(&SIjAb);
  return;
}

}} // namespace psi::cceom
