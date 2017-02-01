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
#include <cstdio>
#include <cstdlib>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "Params.h"
#include "ccwave.h"

namespace psi { namespace ccenergy {


double CCEnergyWavefunction::mp2_energy(void)
{
	/* Note that if we reach this point and ref=1 (ROHF), then we aren't using
	 * semicanonical orbitals and so we can't compute a non-iterative MBPT(2)
	 * energy */
  if(params_.ref == 0) return(rhf_mp2_energy());
  else if(params_.ref == 2) return(uhf_mp2_energy());
  else return 0.0;

}

double CCEnergyWavefunction::rhf_mp2_energy(void)
{
  double T2_energy, T1_energy;
  dpdfile2 F, T1, D1;
  dpdbuf4 T2, D;
  dpdbuf4 S;
  double os_energy, ss_energy, scs_energy;

  /* Initialize MP2 T1 Amps (zero for true HF references) */
  global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 1, "fIA");
  global_dpd_->file2_copy(&F, PSIF_CC_OEI, "tIA (MP2)");
  global_dpd_->file2_close(&F);
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA (MP2)");
  global_dpd_->file2_init(&D1, PSIF_CC_OEI, 0, 0, 1, "dIA");
  global_dpd_->file2_dirprd(&D1, &T1);
  global_dpd_->file2_close(&D1);
  global_dpd_->file2_close(&T1);

  global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 1, "fIA");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA (MP2)");
  T1_energy = 2.0 * global_dpd_->file2_dot(&F, &T1);
  global_dpd_->file2_close(&F);
  global_dpd_->file2_close(&T1);

  /* Initialize MP2 T2 Amps */
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  global_dpd_->buf4_copy(&D, PSIF_CC_TAMPS, "tIjAb (MP2)");
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb (MP2)");
  global_dpd_->buf4_init(&D, PSIF_CC_DENOM, 0, 0, 5, 0, 5, 0, "dIjAb");
  global_dpd_->buf4_dirprd(&D, &T2);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_close(&T2);

  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb (MP2)");
  T2_energy = global_dpd_->buf4_dot(&D, &T2);

  global_dpd_->buf4_init(&S, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  os_energy = global_dpd_->buf4_dot(&S, &T2);
  global_dpd_->buf4_close(&S);
  ss_energy = (T2_energy - os_energy);

  moinfo_.emp2_ss = ss_energy;
  moinfo_.emp2_os = os_energy;

  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&D);

  return (T2_energy+T1_energy);
}

double CCEnergyWavefunction::uhf_mp2_energy(void)
{
  double E2AA, E2BB, E2AB, T1A, T1B;
  dpdbuf4 T2, D;
  dpdfile2 T1, F, D1;

  /* Initialize MP2 T1 Amps (zero for true HF references) */
  global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 1, "fIA");
  global_dpd_->file2_copy(&F, PSIF_CC_OEI, "tIA (MP2)");
  global_dpd_->file2_close(&F);
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA (MP2)");
  global_dpd_->file2_init(&D1, PSIF_CC_OEI, 0, 0, 1, "dIA");
  global_dpd_->file2_dirprd(&D1, &T1);
  global_dpd_->file2_close(&D1);
  global_dpd_->file2_close(&T1);

  global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 2, 3, "fia");
  global_dpd_->file2_copy(&F, PSIF_CC_OEI, "tia (MP2)");
  global_dpd_->file2_close(&F);
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia (MP2)");
  global_dpd_->file2_init(&D1, PSIF_CC_OEI, 0, 2, 3, "dia");
  global_dpd_->file2_dirprd(&D1, &T1);
  global_dpd_->file2_close(&D1);
  global_dpd_->file2_close(&T1);

  global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 1, "fIA");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA (MP2)");
  T1A = global_dpd_->file2_dot(&F, &T1);
  global_dpd_->file2_close(&F);
  global_dpd_->file2_close(&T1);

  global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 2, 3, "fia");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia (MP2)");
  T1B = global_dpd_->file2_dot(&F, &T1);
  global_dpd_->file2_close(&F);
  global_dpd_->file2_close(&T1);

  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 2, 7, 2, 7, 0, "D <IJ||AB> (I>J,A>B)");
  global_dpd_->buf4_copy(&D, PSIF_CC_TAMPS, "tIJAB (MP2)");
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tIJAB (MP2)");
  global_dpd_->buf4_init(&D, PSIF_CC_DENOM, 0, 2, 7, 2, 7, 0, "dIJAB");
  global_dpd_->buf4_dirprd(&D, &T2);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_close(&T2);

  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 12, 17, 12, 17, 0, "D <ij||ab> (i>j,a>b)");
  global_dpd_->buf4_copy(&D, PSIF_CC_TAMPS, "tijab (MP2)");
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 12, 17, 12, 17, 0, "tijab (MP2)");
  global_dpd_->buf4_init(&D, PSIF_CC_DENOM, 0, 12, 17, 12, 17, 0, "dijab");
  global_dpd_->buf4_dirprd(&D, &T2);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_close(&T2);

  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
  global_dpd_->buf4_copy(&D, PSIF_CC_TAMPS, "tIjAb (MP2)");
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb (MP2)");
  global_dpd_->buf4_init(&D, PSIF_CC_DENOM, 0, 22, 28, 22, 28, 0, "dIjAb");
  global_dpd_->buf4_dirprd(&D, &T2);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_close(&T2);

  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tIJAB (MP2)");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 2, 7, 2, 7, 0, "D <IJ||AB> (I>J,A>B)");
  E2AA = global_dpd_->buf4_dot(&D, &T2);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_close(&T2);

  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 12, 17, 12, 17, 0, "tijab (MP2)");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 12, 17, 12, 17, 0, "D <ij||ab> (i>j,a>b)");
  E2BB = global_dpd_->buf4_dot(&D, &T2);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_close(&T2);

  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb (MP2)");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
  E2AB = global_dpd_->buf4_dot(&D, &T2);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_close(&T2);

  // We define EMP2_SS as the same-spin pair energy, and EMP2_OS as the
  // opposite-spin pair energy (singles not included)
  moinfo_.emp2_ss = E2AA + E2BB;
  moinfo_.emp2_os = E2AB;

  return(T1A + T1B + E2AA + E2BB + E2AB);
}
}} // namespace psi::ccenergy
