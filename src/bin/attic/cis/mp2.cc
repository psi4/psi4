/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

/*! \file
    \ingroup CIS
    \brief Enter brief description of file here 
*/
#include <cmath>
#include <cstdlib>
#include <libdpd/dpd.h>
#include <psifiles.h>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cis {

void local_filter_T1(dpdfile2 *T1);
void local_filter_T2(dpdbuf4 *T2);

void mp2(void)
{
  int iter, h, nirreps, row, col;
  double energy, conv, rms, value;
  dpdfile2 F;
  dpdbuf4 D, T2, newT2, Z;

  nirreps = moinfo.nirreps;

  if(params.ref == 0) { /** RHF **/

    /* build initial guess amplitudes */
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->buf4_copy(&D, PSIF_CC_MISC, "MP2 tIjAb");
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&T2, PSIF_CC_MISC, 0, 0, 5, 0, 5, 0, "MP2 tIjAb");
    if(params.local) local_filter_T2(&T2);
    else {
      global_dpd_->buf4_init(&D, PSIF_CC_DENOM, 0, 0, 5, 0, 5, 0, "dIjAb");
      global_dpd_->buf4_dirprd(&D, &T2);
      global_dpd_->buf4_close(&D);
    }

    global_dpd_->buf4_copy(&T2, PSIF_CC_MISC, "New MP2 tIjAb");
    global_dpd_->buf4_close(&T2);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    global_dpd_->buf4_init(&T2, PSIF_CC_MISC, 0, 0, 5, 0, 5, 0, "MP2 tIjAb");
    energy = global_dpd_->buf4_dot(&D, &T2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&D);

    if(params.local) {
      fprintf(outfile, "\n\tSolving for LMP2 wave function:\n");
      fprintf(outfile,   "\t-------------------------------\n");
      fprintf(outfile, "\titer = %d  LMP2 Energy = %20.14f\n", 0, energy);
    }
    else {
      fprintf(outfile, "\n\tSolving for MP2 wave function:\n");
      fprintf(outfile,   "\t-------------------------------\n");
      fprintf(outfile, "\titer = %d  MP2 Energy = %20.14f\n", 0, energy);
    }

    conv = 0;
    for(iter=1; iter < params.maxiter; iter++) {

      global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
      global_dpd_->buf4_copy(&D, PSIF_CC_MISC, "New MP2 tIjAb Increment");
      global_dpd_->buf4_close(&D);

      global_dpd_->buf4_init(&newT2, PSIF_CC_MISC, 0, 0, 5, 0, 5, 0, "New MP2 tIjAb Increment");
      global_dpd_->buf4_init(&T2, PSIF_CC_MISC, 0, 0, 5, 0, 5, 0, "MP2 tIjAb");

      global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 0, "fIJ");
      global_dpd_->contract424(&T2, &F, &newT2, 1, 0, 1, -1, 1);
      global_dpd_->contract244(&F, &T2, &newT2, 0, 0, 0, -1, 1);
      global_dpd_->file2_close(&F);

      global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 1, 1, "fAB");
      global_dpd_->contract244(&F, &T2, &newT2, 1, 2, 1, 1, 1);
      global_dpd_->contract424(&T2, &F, &newT2, 3, 1, 0, 1, 1);
      global_dpd_->file2_close(&F);

      global_dpd_->buf4_close(&T2);

      if(params.local) {
	local_filter_T2(&newT2);
      }
      else {
	global_dpd_->buf4_init(&D, PSIF_CC_DENOM, 0, 0, 5, 0, 5, 0, "dIjAb");
	global_dpd_->buf4_dirprd(&D, &newT2);
	global_dpd_->buf4_close(&D);
      }

      global_dpd_->buf4_close(&newT2);

      global_dpd_->buf4_init(&newT2, PSIF_CC_MISC, 0, 0, 5, 0, 5, 0, "New MP2 tIjAb");
      global_dpd_->buf4_init(&T2, PSIF_CC_MISC, 0, 0, 5, 0, 5, 0, "New MP2 tIjAb Increment");
      global_dpd_->buf4_axpy(&T2, &newT2, 1);
      global_dpd_->buf4_close(&T2);

      global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
      energy = global_dpd_->buf4_dot(&D, &newT2);
      global_dpd_->buf4_close(&D);
      global_dpd_->buf4_close(&newT2);

      global_dpd_->buf4_init(&newT2, PSIF_CC_MISC, 0, 0, 5, 0, 5, 0, "New MP2 tIjAb");
      global_dpd_->buf4_init(&T2, PSIF_CC_MISC, 0, 0, 5, 0, 5, 0, "MP2 tIjAb");
      rms = 0.0;
      for(h=0; h < nirreps; h++) {
	global_dpd_->buf4_mat_irrep_init(&newT2, h);
	global_dpd_->buf4_mat_irrep_rd(&newT2, h);
	global_dpd_->buf4_mat_irrep_init(&T2, h);
	global_dpd_->buf4_mat_irrep_rd(&T2, h);

	for(row=0; row < T2.params->rowtot[h]; row++)
	  for(col=0; col < T2.params->coltot[h]; col++) {
	    value = newT2.matrix[h][row][col] - T2.matrix[h][row][col];
	    rms += value * value;
	  }

	global_dpd_->buf4_mat_irrep_close(&T2, h);
	global_dpd_->buf4_mat_irrep_close(&newT2, h);
      }
      global_dpd_->buf4_close(&T2);
      global_dpd_->buf4_close(&newT2);
      rms = sqrt(rms);

      if(params.local) {
	fprintf(outfile, "\titer = %d   LMP2 Energy = %20.14f   RMS = %4.3e\n", iter, energy, rms);
      }
      else {
	fprintf(outfile, "\titer = %d   MP2 Energy = %20.14f   RMS = %4.3e\n", iter, energy, rms);
      }

      if(rms < params.convergence) {
	conv = 1;
	fprintf(outfile, "\n\tMP2 iterations converged.\n\n");
	break;
      }
      else {
	global_dpd_->buf4_init(&T2, PSIF_CC_MISC, 0, 0, 5, 0, 5, 0, "New MP2 tIjAb");
	global_dpd_->buf4_copy(&T2, PSIF_CC_MISC, "MP2 tIjAb");
	global_dpd_->buf4_close(&T2);
      }
    }

    if(!conv) {
      fprintf(outfile, "\n\tMP2 iterative procedure failed.\n");
      throw PsiException("cis MP2 iteration error", __FILE__, __LINE__);
    }

    /* spin adapt the final amplitudes */
    global_dpd_->buf4_init(&T2, PSIF_CC_MISC, 0, 0, 5, 0, 5, 0, "MP2 tIjAb");
    global_dpd_->buf4_sort(&T2, PSIF_CC_TMP0, pqsr, 0, 5, "MP2 tIjbA");
    global_dpd_->buf4_copy(&T2, PSIF_CC_MISC, "MP2 2 tIjAb - tIjbA");
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_MISC, 0, 0, 5, 0, 5, 0, "MP2 2 tIjAb - tIjbA");
    global_dpd_->buf4_scm(&T2, 2);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "MP2 tIjbA");
    global_dpd_->buf4_axpy(&Z, &T2, -1);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&T2);
  }
  else if(params.ref == 2) { /** UHF **/
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 2, 7, 2, 7, 0, "D <IJ||AB> (I>J,A>B)");
    global_dpd_->buf4_copy(&D, PSIF_CC_MISC, "MP2 tIJAB");
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&T2, PSIF_CC_MISC, 0, 2, 7, 2, 7, 0, "MP2 tIJAB");
    global_dpd_->buf4_init(&D, PSIF_CC_DENOM, 0, 1, 6, 1, 6, 0, "dIJAB");
    global_dpd_->buf4_dirprd(&D, &T2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&T2);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 12, 17, 12, 17, 0, "D <ij||ab> (i>j,a>b)");
    global_dpd_->buf4_copy(&D, PSIF_CC_MISC, "MP2 tijab");
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&T2, PSIF_CC_MISC, 0, 12, 17, 12, 17, 0, "MP2 tijab");
    global_dpd_->buf4_init(&D, PSIF_CC_DENOM, 0, 11, 16, 11, 16, 0, "dijab");
    global_dpd_->buf4_dirprd(&D, &T2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&T2);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    global_dpd_->buf4_copy(&D, PSIF_CC_MISC, "MP2 tIjAb");
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&T2, PSIF_CC_MISC, 0, 22, 28, 22, 28, 0, "MP2 tIjAb");
    global_dpd_->buf4_init(&D, PSIF_CC_DENOM, 0, 22, 28, 22, 28, 0, "dIjAb");
    global_dpd_->buf4_dirprd(&D, &T2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&T2);
  }
}


}} // namespace psi::cis
