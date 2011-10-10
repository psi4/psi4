/*! \file
    \ingroup CCENERGY
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <libdpd/dpd.h>
#include <libpsio/psio.h>
#include <psifiles.h>
#include "Params.h"
#include "Local.h"
#include "MOInfo.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccenergy {

void local_filter_T2(dpdbuf4 *T2);

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

void lmp2(void)
{
  int i, j, k, ij, ab, iter, conv, row, col, nocc, nvir, natom, weak;
  double energy, rms, weak_pair_energy;
  dpdbuf4 T2, newT2, D;
  dpdfile2 fij, fab;
  psio_address next;

  nocc = local.nocc;
  nvir = local.nvir;
  natom = local.natom;

  local.domain = (int **) malloc(local.nocc*sizeof(int *));
  next = PSIO_ZERO;
  for(i=0; i<nocc; i++) {
    local.domain[i] = (int *) malloc(local.natom*sizeof(int));
    psio_read(CC_INFO, "Local Domains", (char *) local.domain[i],
              natom*sizeof(int), next, &next);
  }

  /* First, turn on all weak pairs for the LMP2 */
  for(ij=0; ij < nocc*nocc; ij++) local.weak_pairs[ij] = 0;

  /* Clean out diagonal element of occ-occ Fock matrix */
  dpd_file2_init(&fij, CC_OEI, 0, 0, 0, "fIJ");
  dpd_file2_copy(&fij, CC_OEI, "fIJ (non-diagonal)");
  dpd_file2_close(&fij);

  dpd_file2_init(&fij, CC_OEI, 0, 0, 0, "fIJ (non-diagonal)");
  dpd_file2_mat_init(&fij);
  dpd_file2_mat_rd(&fij);
  for(i=0; i < nocc; i++) fij.matrix[0][i][i] = 0.0;
  dpd_file2_mat_wrt(&fij);
  dpd_file2_close(&fij);

  /* Build initial LMP2 amplitudes */
  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_buf4_copy(&D, CC_TAMPS, "LMP2 tIjAb");
  dpd_buf4_close(&D);

  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "LMP2 tIjAb");
  if(params.local) {
    local_filter_T2(&T2);
  }
  else {
    dpd_buf4_init(&D, CC_DENOM, 0, 0, 5, 0, 5, 0, "dIjAb");
    dpd_buf4_dirprd(&D, &T2);
    dpd_buf4_close(&D);
  }
  dpd_buf4_close(&T2);

  /* Compute the LMP2 energy */
  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "LMP2 tIjAb");
  energy = dpd_buf4_dot(&D, &T2);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&D);

  fprintf(outfile, "\n\tComputing LMP2 amplitudes:\n");
  fprintf(outfile,   "\t--------------------------\n");
  fprintf(outfile, "\titer = %d    LMP2 Energy = %20.14f\n", 0, energy);

  conv = 0;
  int lmp2_maxiter=1000;
  for(iter=1; iter < lmp2_maxiter; iter++) {

    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_buf4_copy(&D, CC_TAMPS, "New LMP2 tIjAb Increment");
    dpd_buf4_close(&D);

    dpd_buf4_init(&newT2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New LMP2 tIjAb Increment");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "LMP2 tIjAb");

    dpd_file2_init(&fij, CC_OEI, 0, 0, 0, "fIJ");
    dpd_contract424(&T2, &fij, &newT2, 1, 0, 1, -1, 1);
    dpd_contract244(&fij, &T2, &newT2, 0, 0, 0, -1, 1);
    dpd_file2_close(&fij);

    dpd_file2_init(&fab, CC_OEI, 0, 1, 1, "fAB");
    dpd_contract244(&fab, &T2, &newT2, 1, 2, 1, 1, 1);
    dpd_contract424(&T2, &fab, &newT2, 3, 1, 0, 1, 1);
    dpd_file2_close(&fab);

    dpd_buf4_copy(&T2, CC_TAMPS, "New LMP2 tIjAb");
    dpd_buf4_close(&T2);

    if(params.local) {
      local_filter_T2(&newT2);
    }
    else {
      dpd_buf4_init(&D, CC_DENOM, 0, 0, 5, 0, 5, 0, "dIjAb");
      dpd_buf4_dirprd(&D, &newT2);
      dpd_buf4_close(&D);
    }
    dpd_buf4_close(&newT2);

    dpd_buf4_init(&newT2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New LMP2 tIjAb");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New LMP2 tIjAb Increment");
    dpd_buf4_axpy(&T2, &newT2, 1);
    dpd_buf4_close(&T2);

    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    energy = dpd_buf4_dot(&D, &newT2);
    dpd_buf4_close(&D);

    dpd_buf4_close(&newT2);

    /* Check for convergence */
    dpd_buf4_init(&newT2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New LMP2 tIjAb");
    dpd_buf4_mat_irrep_init(&newT2, 0);
    dpd_buf4_mat_irrep_rd(&newT2, 0);

    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "LMP2 tIjAb");
    dpd_buf4_mat_irrep_init(&T2, 0);
    dpd_buf4_mat_irrep_rd(&T2, 0);

    rms = 0.0;
    for(row=0; row < T2.params->rowtot[0]; row++)
      for(col=0; col < T2.params->coltot[0]; col++)
        rms += (newT2.matrix[0][row][col] - T2.matrix[0][row][col]) *
          (newT2.matrix[0][row][col] - T2.matrix[0][row][col]);

    dpd_buf4_mat_irrep_close(&T2, 0);
    dpd_buf4_mat_irrep_close(&newT2, 0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&newT2);

    rms = sqrt(rms);

    fprintf(outfile, "\titer = %d    LMP2 Energy = %20.14f   RMS = %4.3e\n", iter, energy, rms);

    if(rms < params.convergence) {
      conv = 1;
      fprintf(outfile, "\n\tLMP2 Iterations converged.\n");
      break;
    }
    else {
      dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New LMP2 tIjAb");
      dpd_buf4_copy(&T2, CC_TAMPS, "LMP2 tIjAb");
      dpd_buf4_close(&T2);
    }
  }

  if(!conv) {
    fprintf(outfile, "\n\tLMP2 Iterative procedure failed.\n");
    throw ConvergenceError<int>("LMP2 interative procedure failed.", lmp2_maxiter, params.convergence, rms, __FILE__, __LINE__);
  }

  /* Turn off weak pairs again for the LCCSD */

  for(i=0,ij=0; i < nocc; i++)
    for(j=0; j < nocc; j++,ij++) {
      weak = 1;
      for(k=0; k < natom; k++)
        if(local.domain[i][k] && local.domain[j][k]) weak = 0;

      if(weak) local.weak_pairs[ij] = 1;
      else local.weak_pairs[ij] = 0;
    }

  /* Compute the MP2 weak-pair energy */
  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
  dpd_buf4_mat_irrep_init(&D, 0);
  dpd_buf4_mat_irrep_rd(&D, 0);
  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New LMP2 tIjAb");
  dpd_buf4_mat_irrep_init(&T2, 0);
  dpd_buf4_mat_irrep_rd(&T2, 0);

  weak_pair_energy = 0.0;
  for(ij=0; ij < nocc*nocc; ij++)
    if(local.weak_pairs[ij])
      for(ab=0; ab < nvir*nvir; ab++)
        weak_pair_energy += D.matrix[0][ij][ab] * T2.matrix[0][ij][ab];

  dpd_buf4_mat_irrep_close(&T2, 0);
  dpd_buf4_close(&T2);
  dpd_buf4_mat_irrep_close(&D, 0);
  dpd_buf4_close(&D);

  fprintf(outfile, "\n\tLMP2 Weak Pair Energy   = %20.14f\n", weak_pair_energy);
  fprintf(outfile, "\tLMP2 Correlation Energy = %20.14f\n", energy);
  fprintf(outfile, "\tLMP2 Total Energy       = %20.14f\n\n", energy+moinfo.eref);
  fflush(outfile);

  local.weak_pair_energy = weak_pair_energy;

  for(i=0; i<nocc; i++)
    free(local.domain[i]);
  free(local.domain);
}
}} // namespace psi::ccenergy
