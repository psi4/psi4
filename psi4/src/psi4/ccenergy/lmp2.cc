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
#include <cmath>
#include "psi4/libdpd/dpd.h"
#include "psi4/libpsio/psio.h"
#include "psi4/psifiles.h"
#include "Params.h"
#include "Local.h"
#include "MOInfo.h"
#include "ccwave.h"

namespace psi { namespace ccenergy {

/* lmp2(): Computes the local-MP2 energy and the local-MP2 weak-pair energy.
**
** Given a set of non-canonical occupied orbitals and canonical virtual
** orbitals, this routine computes the local-MP2 energy as described originally
** by Saebo and Pulay, J. Chem. Phys. 86, 914 (1987).
**
** The primary purpose of this code is actually to compute the
** so-called "weak pair" energy for use in a subsequent local-CCSD
** calculation.  Given a set of t_ij^ab doubles-excitation cluster
** amplitudes, a pair ij is considered weak if the domains of atoms
** belonging to occupied orbitals i and j contain no atoms in common.
** These pairs are originally identified in local_init(), which must
** be run prior to this routine.
**
** In this code, we first compute the LMP2 energy treating all pairs
** ij equivalently.  Due to the fact that the occupied orbitals are
** non-canonical, the LMP2 amplitude equations must be solved
** iteratively (i.e., the amplitudes are coupled).  The weak-pair
** energy is defined as the MP2 energy contribution associated with only
** the weak pairs.  This energy is used later to correct the LCCSD
** energy in which all such pairs are completely neglected.
**
** TDC, June 2002
*/

void CCEnergyWavefunction::lmp2(void)
{
  int i, j, k, ij, ab, iter, conv, row, col, nocc, nvir, natom, weak;
  double energy, rms, weak_pair_energy;
  dpdbuf4 T2, newT2, D;
  dpdfile2 fij, fab;
  psio_address next;

  nocc = local_.nocc;
  nvir = local_.nvir;
  natom = local_.natom;

  local_.domain = (int **) malloc(local_.nocc*sizeof(int *));
  next = PSIO_ZERO;
  for(i=0; i<nocc; i++) {
    local_.domain[i] = (int *) malloc(local_.natom*sizeof(int));
    psio_read(PSIF_CC_INFO, "Local Domains", (char *) local_.domain[i],
              natom*sizeof(int), next, &next);
  }

  /* First, turn on all weak pairs for the LMP2 */
  for(ij=0; ij < nocc*nocc; ij++) local_.weak_pairs[ij] = 0;

  /* Clean out diagonal element of occ-occ Fock matrix */
  global_dpd_->file2_init(&fij, PSIF_CC_OEI, 0, 0, 0, "fIJ");
  global_dpd_->file2_copy(&fij, PSIF_CC_OEI, "fIJ (non-diagonal)");
  global_dpd_->file2_close(&fij);

  global_dpd_->file2_init(&fij, PSIF_CC_OEI, 0, 0, 0, "fIJ (non-diagonal)");
  global_dpd_->file2_mat_init(&fij);
  global_dpd_->file2_mat_rd(&fij);
  for(i=0; i < nocc; i++) fij.matrix[0][i][i] = 0.0;
  global_dpd_->file2_mat_wrt(&fij);
  global_dpd_->file2_close(&fij);

  /* Build initial LMP2 amplitudes */
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  global_dpd_->buf4_copy(&D, PSIF_CC_TAMPS, "LMP2 tIjAb");
  global_dpd_->buf4_close(&D);

  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "LMP2 tIjAb");
  if(params_.local) {
    local_filter_T2(&T2);
  }
  else {
    global_dpd_->buf4_init(&D, PSIF_CC_DENOM, 0, 0, 5, 0, 5, 0, "dIjAb");
    global_dpd_->buf4_dirprd(&D, &T2);
    global_dpd_->buf4_close(&D);
  }
  global_dpd_->buf4_close(&T2);

  /* Compute the LMP2 energy */
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "LMP2 tIjAb");
  energy = global_dpd_->buf4_dot(&D, &T2);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&D);

  outfile->Printf( "\n    Computing LMP2 amplitudes:\n");
  outfile->Printf(   "    --------------------------\n");
  outfile->Printf( "    iter = %d    LMP2 Energy = %20.14f\n", 0, energy);

  conv = 0;
  int lmp2_maxiter=1000;
  for(iter=1; iter < lmp2_maxiter; iter++) {

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->buf4_copy(&D, PSIF_CC_TAMPS, "New LMP2 tIjAb Increment");
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&newT2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New LMP2 tIjAb Increment");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "LMP2 tIjAb");

    global_dpd_->file2_init(&fij, PSIF_CC_OEI, 0, 0, 0, "fIJ");
    global_dpd_->contract424(&T2, &fij, &newT2, 1, 0, 1, -1, 1);
    global_dpd_->contract244(&fij, &T2, &newT2, 0, 0, 0, -1, 1);
    global_dpd_->file2_close(&fij);

    global_dpd_->file2_init(&fab, PSIF_CC_OEI, 0, 1, 1, "fAB");
    global_dpd_->contract244(&fab, &T2, &newT2, 1, 2, 1, 1, 1);
    global_dpd_->contract424(&T2, &fab, &newT2, 3, 1, 0, 1, 1);
    global_dpd_->file2_close(&fab);

    global_dpd_->buf4_copy(&T2, PSIF_CC_TAMPS, "New LMP2 tIjAb");
    global_dpd_->buf4_close(&T2);

    if(params_.local) {
      local_filter_T2(&newT2);
    }
    else {
      global_dpd_->buf4_init(&D, PSIF_CC_DENOM, 0, 0, 5, 0, 5, 0, "dIjAb");
      global_dpd_->buf4_dirprd(&D, &newT2);
      global_dpd_->buf4_close(&D);
    }
    global_dpd_->buf4_close(&newT2);

    global_dpd_->buf4_init(&newT2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New LMP2 tIjAb");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New LMP2 tIjAb Increment");
    global_dpd_->buf4_axpy(&T2, &newT2, 1);
    global_dpd_->buf4_close(&T2);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    energy = global_dpd_->buf4_dot(&D, &newT2);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_close(&newT2);

    /* Check for convergence */
    global_dpd_->buf4_init(&newT2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New LMP2 tIjAb");
    global_dpd_->buf4_mat_irrep_init(&newT2, 0);
    global_dpd_->buf4_mat_irrep_rd(&newT2, 0);

    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "LMP2 tIjAb");
    global_dpd_->buf4_mat_irrep_init(&T2, 0);
    global_dpd_->buf4_mat_irrep_rd(&T2, 0);

    rms = 0.0;
    for(row=0; row < T2.params->rowtot[0]; row++)
      for(col=0; col < T2.params->coltot[0]; col++)
        rms += (newT2.matrix[0][row][col] - T2.matrix[0][row][col]) *
          (newT2.matrix[0][row][col] - T2.matrix[0][row][col]);

    global_dpd_->buf4_mat_irrep_close(&T2, 0);
    global_dpd_->buf4_mat_irrep_close(&newT2, 0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&newT2);

    rms = sqrt(rms);

    outfile->Printf( "    iter = %d    LMP2 Energy = %20.14f   RMS = %4.3e\n", iter, energy, rms);

    if(rms < params_.convergence) {
      conv = 1;
      outfile->Printf( "\n    LMP2 Iterations converged.\n");
      break;
    }
    else {
      global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New LMP2 tIjAb");
      global_dpd_->buf4_copy(&T2, PSIF_CC_TAMPS, "LMP2 tIjAb");
      global_dpd_->buf4_close(&T2);
    }
  }

  if(!conv) {
    outfile->Printf( "\n    LMP2 Iterative procedure failed.\n");
    throw ConvergenceError<int>("LMP2 interative procedure failed.", lmp2_maxiter, params_.convergence, rms, __FILE__, __LINE__);
  }

  /* Turn off weak pairs again for the LCCSD */

  for(i=0,ij=0; i < nocc; i++)
    for(j=0; j < nocc; j++,ij++) {
      weak = 1;
      for(k=0; k < natom; k++)
        if(local_.domain[i][k] && local_.domain[j][k]) weak = 0;

      if(weak) local_.weak_pairs[ij] = 1;
      else local_.weak_pairs[ij] = 0;
    }

  /* Compute the MP2 weak-pair energy */
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
  global_dpd_->buf4_mat_irrep_init(&D, 0);
  global_dpd_->buf4_mat_irrep_rd(&D, 0);
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New LMP2 tIjAb");
  global_dpd_->buf4_mat_irrep_init(&T2, 0);
  global_dpd_->buf4_mat_irrep_rd(&T2, 0);

  weak_pair_energy = 0.0;
  for(ij=0; ij < nocc*nocc; ij++)
    if(local_.weak_pairs[ij])
      for(ab=0; ab < nvir*nvir; ab++)
        weak_pair_energy += D.matrix[0][ij][ab] * T2.matrix[0][ij][ab];

  global_dpd_->buf4_mat_irrep_close(&T2, 0);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_mat_irrep_close(&D, 0);
  global_dpd_->buf4_close(&D);

  outfile->Printf( "\n    LMP2 Weak Pair Energy   = %20.14f\n", weak_pair_energy);
  outfile->Printf( "    LMP2 Correlation Energy = %20.14f\n", energy);
  outfile->Printf( "    LMP2 Total Energy       = %20.14f\n\n", energy+moinfo_.eref);


  local_.weak_pair_energy = weak_pair_energy;

  for(i=0; i<nocc; i++)
    free(local_.domain[i]);
  free(local_.domain);
}
}} // namespace psi::ccenergy
