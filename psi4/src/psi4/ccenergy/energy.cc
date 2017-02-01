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

double CCEnergyWavefunction::energy(void)
{
  double e = 0.0;
  if(params_.ref == 0) e = rhf_energy();
  else if(params_.ref == 1) e = rohf_energy();
  else if(params_.ref == 2) e = uhf_energy();

  return e;
}

double CCEnergyWavefunction::rhf_energy(void)
{
  double tauIjAb_energy, tIA_energy;
  dpdfile2 fIA, tIA;
  dpdbuf4 tauIjAb, D, E;
  dpdbuf4 S;
  double os_energy, ss_energy, Energy;

  global_dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
  global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
  /*   dpd_file2_print(&tIA, outfile); */
  tIA_energy = 2.0 * global_dpd_->file2_dot(&fIA, &tIA);
  global_dpd_->file2_close(&fIA);
  global_dpd_->file2_close(&tIA);

  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
  global_dpd_->buf4_init(&tauIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
  tauIjAb_energy = global_dpd_->buf4_dot(&D, &tauIjAb);

  global_dpd_->buf4_init(&S, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  os_energy = global_dpd_->buf4_dot(&S, &tauIjAb);
  ss_energy = (tauIjAb_energy - os_energy);

  moinfo_.ecc_ss = ss_energy;
  moinfo_.ecc_os = os_energy;

  global_dpd_->buf4_close(&S);
  global_dpd_->buf4_close(&tauIjAb);
  global_dpd_->buf4_close(&D);

  /*
    outfile->Printf( "Two AB Energy = %20.14f\n", tauIjAb_energy);
  */

  return (tauIjAb_energy+tIA_energy);
}

double CCEnergyWavefunction::rohf_energy(void)
{
  double tIA_energy, tia_energy, tauIJAB_energy, tauijab_energy, tauIjAb_energy;
  dpdfile2 tIA, tia, fIA, fia;
  dpdbuf4 tauIJAB, tauijab, tauIjAb, D;

  global_dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
  global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
/*  dpd_file2_print(&tIA, outfile);  */
  tIA_energy = global_dpd_->file2_dot(&fIA, &tIA);
  global_dpd_->file2_close(&fIA);
  global_dpd_->file2_close(&tIA);

  global_dpd_->file2_init(&fia, PSIF_CC_OEI, 0, 0, 1, "fia");
  global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");
/*  dpd_file2_print(&tia, outfile); */
  tia_energy = global_dpd_->file2_dot(&fia, &tia);
  global_dpd_->file2_close(&fia);
  global_dpd_->file2_close(&tia);

  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
  global_dpd_->buf4_init(&tauIJAB, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
/*  dpd_buf4_print(&tauIJAB, outfile);  */
  tauIJAB_energy = global_dpd_->buf4_dot(&D, &tauIJAB);
  global_dpd_->buf4_close(&tauIJAB);
  global_dpd_->buf4_init(&tauijab, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauijab");
/*  dpd_buf4_print(&tauijab, outfile); */
  tauijab_energy = global_dpd_->buf4_dot(&D, &tauijab);
  global_dpd_->buf4_close(&tauijab);
  global_dpd_->buf4_close(&D);

  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  global_dpd_->buf4_init(&tauIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
/*  dpd_buf4_print(&tauIjAb, outfile);  */
  tauIjAb_energy = global_dpd_->buf4_dot(&D, &tauIjAb);
  global_dpd_->buf4_close(&tauIjAb);
  global_dpd_->buf4_close(&D);

  /*
  outfile->Printf( "One A Energy = %20.14f\n", tIA_energy);
  outfile->Printf( "One B Energy = %20.14f\n", tia_energy);
  outfile->Printf( "Two AA Energy = %20.14f\n", tauIJAB_energy);
  outfile->Printf( "Two BB Energy = %20.14f\n", tauijab_energy);
  outfile->Printf( "Two AB Energy = %20.14f\n", tauIjAb_energy);
  */

  // Store the same-spin and opposite-spin pair energies
  // (not including singles here)
  moinfo_.ecc_ss = tauIJAB_energy + tauijab_energy;
  moinfo_.ecc_os = tauIjAb_energy;

  return (tIA_energy + tia_energy +
      tauIJAB_energy + tauijab_energy + tauIjAb_energy);
}

double CCEnergyWavefunction::uhf_energy(void)
{
  double E2AA, E2BB, E2AB, T1A, T1B;
  dpdbuf4 T2, D;
  dpdfile2 T1, F;

  global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 1, "fIA");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
/*  dpd_file2_print(&tIA, outfile);  */
  T1A = global_dpd_->file2_dot(&F, &T1);
  global_dpd_->file2_close(&F);
  global_dpd_->file2_close(&T1);

  global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 2, 3, "fia");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
/*  dpd_file2_print(&tIA, outfile);  */
  T1B = global_dpd_->file2_dot(&F, &T1);
  global_dpd_->file2_close(&F);
  global_dpd_->file2_close(&T1);

  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 2, 7, 2, 7, 0, "D <IJ||AB> (I>J,A>B)");
  E2AA = global_dpd_->buf4_dot(&D, &T2);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_close(&T2);

  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 12, 17, 12, 17, 0, "tauijab");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 12, 17, 12, 17, 0, "D <ij||ab> (i>j,a>b)");
  E2BB = global_dpd_->buf4_dot(&D, &T2);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_close(&T2);

  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tauIjAb");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
  E2AB = global_dpd_->buf4_dot(&D, &T2);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_close(&T2);

  /*
  outfile->Printf( "One A Energy = %20.14f\n", T1A);
  outfile->Printf( "One B Energy = %20.14f\n", T1B);
  outfile->Printf( "Two AA Energy = %20.14f\n", E2AA);
  outfile->Printf( "Two BB Energy = %20.14f\n", E2BB);
  outfile->Printf( "Two AB Energy = %20.14f\n", E2AB);
  */

  /*
  outfile->Printf("\n    Opposite-spin energy  = %20.15f\n",E2AB);
  outfile->Printf("    Same-spin energy  = %20.15f\n",E2AA+E2BB);
  */

  // Store the same-spin and opposite-spin pair energies
  // (not including singles here)
  moinfo_.ecc_ss = E2AA + E2BB;
  moinfo_.ecc_os = E2AB;

  return(T1A + T1B + E2AA + E2BB + E2AB);
}
}} // namespace psi::ccenergy
