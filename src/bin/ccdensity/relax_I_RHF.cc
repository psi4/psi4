/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here 
*/
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

/* relax_I_RHF(): Add the ROHF orbital-response contributions from
** the one-electron density matrix to the I(I,J) and I(I,A) blocks of
** the Lagrangian.  These terms arise from the first-order CPHF
** equations.  I *think* the following code is general enough to deal
** with both RHF and ROHF cases. */

void relax_I_RHF(void)
{
  dpdfile2 I, D, f;
  dpdbuf4 E;
  int h, nirreps, i, j, e, *occpi, *virtpi, *openpi;

  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi;
  virtpi = moinfo.virtpi;
  openpi = moinfo.openpi;

  /* I(I,A) = I'(I,A) - sum_M f(I,M) D(orb)(A,M) */
  dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I'IA");
  dpd_file2_copy(&I, CC_OEI, "I(I,A)");
  dpd_file2_close(&I);
  dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I(I,A)");
  dpd_file2_init(&D, CC_OEI, 0, 1, 0, "D(orb)(A,I)");
  dpd_file2_init(&f, CC_OEI, 0, 0, 0, "fIJ");
  dpd_contract222(&f, &D, &I, 0, 0, -1.0, 1.0);
  dpd_file2_close(&f);
  dpd_file2_close(&D);
  dpd_file2_close(&I);

  /* RHF Case: I(i,j) = I'(i,j) - D(orb)(e,c) [4 <ei|mj> - <ei|jm> - <ej|im>] */
  dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'IJ");
  dpd_file2_copy(&I, CC_OEI, "I(I,J)");
  dpd_file2_close(&I);

  dpd_buf4_init(&E, CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
  dpd_buf4_scmcopy(&E, CC_EINTS, "4 <ei|mj> - <ei|jm> - <ej|im>", 4);
  dpd_buf4_sort_axpy(&E, CC_EINTS, pqsr, 11, 0, "4 <ei|mj> - <ei|jm> - <ej|im>", -1);
  dpd_buf4_sort_axpy(&E, CC_EINTS, psqr, 11, 0, "4 <ei|mj> - <ei|jm> - <ej|im>", -1);
  dpd_buf4_close(&E);

  dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I(I,J)");
  dpd_file2_init(&D, CC_OEI, 0, 1, 0, "D(orb)(A,I)");
  dpd_buf4_init(&E, CC_EINTS, 0, 11, 0, 11, 0, 0, "4 <ei|mj> - <ei|jm> - <ej|im>");
  dpd_dot13(&D, &E, &I, 0, 0, -1.0, 1.0);
  dpd_buf4_close(&E);
  dpd_file2_close(&D);

  dpd_file2_close(&I);
}

}} // namespace psi::ccdensity
