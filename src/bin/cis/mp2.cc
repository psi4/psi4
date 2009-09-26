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
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_buf4_copy(&D, CC_MISC, "MP2 tIjAb");
    dpd_buf4_close(&D);

    dpd_buf4_init(&T2, CC_MISC, 0, 0, 5, 0, 5, 0, "MP2 tIjAb");
    if(params.local) local_filter_T2(&T2);
    else {
      dpd_buf4_init(&D, CC_DENOM, 0, 0, 5, 0, 5, 0, "dIjAb");
      dpd_buf4_dirprd(&D, &T2);
      dpd_buf4_close(&D);
    }

    dpd_buf4_copy(&T2, CC_MISC, "New MP2 tIjAb");
    dpd_buf4_close(&T2);

    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    dpd_buf4_init(&T2, CC_MISC, 0, 0, 5, 0, 5, 0, "MP2 tIjAb");
    energy = dpd_buf4_dot(&D, &T2);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&D);

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

      dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
      dpd_buf4_copy(&D, CC_MISC, "New MP2 tIjAb Increment");
      dpd_buf4_close(&D);

      dpd_buf4_init(&newT2, CC_MISC, 0, 0, 5, 0, 5, 0, "New MP2 tIjAb Increment");
      dpd_buf4_init(&T2, CC_MISC, 0, 0, 5, 0, 5, 0, "MP2 tIjAb");

      dpd_file2_init(&F, CC_OEI, 0, 0, 0, "fIJ");
      dpd_contract424(&T2, &F, &newT2, 1, 0, 1, -1, 1);
      dpd_contract244(&F, &T2, &newT2, 0, 0, 0, -1, 1);
      dpd_file2_close(&F);

      dpd_file2_init(&F, CC_OEI, 0, 1, 1, "fAB");
      dpd_contract244(&F, &T2, &newT2, 1, 2, 1, 1, 1);
      dpd_contract424(&T2, &F, &newT2, 3, 1, 0, 1, 1);
      dpd_file2_close(&F);

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

      dpd_buf4_init(&newT2, CC_MISC, 0, 0, 5, 0, 5, 0, "New MP2 tIjAb");
      dpd_buf4_init(&T2, CC_MISC, 0, 0, 5, 0, 5, 0, "New MP2 tIjAb Increment");
      dpd_buf4_axpy(&T2, &newT2, 1);
      dpd_buf4_close(&T2);

      dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
      energy = dpd_buf4_dot(&D, &newT2);
      dpd_buf4_close(&D);
      dpd_buf4_close(&newT2);

      dpd_buf4_init(&newT2, CC_MISC, 0, 0, 5, 0, 5, 0, "New MP2 tIjAb");
      dpd_buf4_init(&T2, CC_MISC, 0, 0, 5, 0, 5, 0, "MP2 tIjAb");
      rms = 0.0;
      for(h=0; h < nirreps; h++) {
	dpd_buf4_mat_irrep_init(&newT2, h);
	dpd_buf4_mat_irrep_rd(&newT2, h);
	dpd_buf4_mat_irrep_init(&T2, h);
	dpd_buf4_mat_irrep_rd(&T2, h);

	for(row=0; row < T2.params->rowtot[h]; row++)
	  for(col=0; col < T2.params->coltot[h]; col++) {
	    value = newT2.matrix[h][row][col] - T2.matrix[h][row][col];
	    rms += value * value;
	  }

	dpd_buf4_mat_irrep_close(&T2, h);
	dpd_buf4_mat_irrep_close(&newT2, h);
      }
      dpd_buf4_close(&T2);
      dpd_buf4_close(&newT2);
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
	dpd_buf4_init(&T2, CC_MISC, 0, 0, 5, 0, 5, 0, "New MP2 tIjAb");
	dpd_buf4_copy(&T2, CC_MISC, "MP2 tIjAb");
	dpd_buf4_close(&T2);
      }
    }

    if(!conv) {
      fprintf(outfile, "\n\tMP2 iterative procedure failed.\n");
      throw PsiException("cis MP2 iteration error", __FILE__, __LINE__);
    }

    /* spin adapt the final amplitudes */
    dpd_buf4_init(&T2, CC_MISC, 0, 0, 5, 0, 5, 0, "MP2 tIjAb");
    dpd_buf4_sort(&T2, CC_TMP0, pqsr, 0, 5, "MP2 tIjbA");
    dpd_buf4_copy(&T2, CC_MISC, "MP2 2 tIjAb - tIjbA");
    dpd_buf4_close(&T2);
    dpd_buf4_init(&T2, CC_MISC, 0, 0, 5, 0, 5, 0, "MP2 2 tIjAb - tIjbA");
    dpd_buf4_scm(&T2, 2);
    dpd_buf4_init(&Z, CC_TMP0, 0, 0, 5, 0, 5, 0, "MP2 tIjbA");
    dpd_buf4_axpy(&Z, &T2, -1);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&T2);
  }
  else if(params.ref == 2) { /** UHF **/
    dpd_buf4_init(&D, CC_DINTS, 0, 2, 7, 2, 7, 0, "D <IJ||AB> (I>J,A>B)");
    dpd_buf4_copy(&D, CC_MISC, "MP2 tIJAB");
    dpd_buf4_close(&D);
    dpd_buf4_init(&T2, CC_MISC, 0, 2, 7, 2, 7, 0, "MP2 tIJAB");
    dpd_buf4_init(&D, CC_DENOM, 0, 1, 6, 1, 6, 0, "dIJAB");
    dpd_buf4_dirprd(&D, &T2);
    dpd_buf4_close(&D);
    dpd_buf4_close(&T2);

    dpd_buf4_init(&D, CC_DINTS, 0, 12, 17, 12, 17, 0, "D <ij||ab> (i>j,a>b)");
    dpd_buf4_copy(&D, CC_MISC, "MP2 tijab");
    dpd_buf4_close(&D);
    dpd_buf4_init(&T2, CC_MISC, 0, 12, 17, 12, 17, 0, "MP2 tijab");
    dpd_buf4_init(&D, CC_DENOM, 0, 11, 16, 11, 16, 0, "dijab");
    dpd_buf4_dirprd(&D, &T2);
    dpd_buf4_close(&D);
    dpd_buf4_close(&T2);

    dpd_buf4_init(&D, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    dpd_buf4_copy(&D, CC_MISC, "MP2 tIjAb");
    dpd_buf4_close(&D);
    dpd_buf4_init(&T2, CC_MISC, 0, 22, 28, 22, 28, 0, "MP2 tIjAb");
    dpd_buf4_init(&D, CC_DENOM, 0, 22, 28, 22, 28, 0, "dIjAb");
    dpd_buf4_dirprd(&D, &T2);
    dpd_buf4_close(&D);
    dpd_buf4_close(&T2);
  }
}


}} // namespace psi::cis
