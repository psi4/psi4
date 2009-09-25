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

/* relax_I_ROHF(): Add the ROHF orbital-response contributions from
** the one-electron density matrix to the I(I,J) and I(I,A) blocks of
** the Lagrangian.  These terms arise from the first-order CPHF
** equations.  I *think* the following code is general enough to deal
** with both RHF and ROHF cases. */

void relax_I_ROHF(void)
{
  dpdfile2 I, D, f;
  dpdbuf4 E;
  int h, nirreps, i, j, e, *occpi, *virtpi, *openpi;

  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi;
  virtpi = moinfo.virtpi;
  openpi = moinfo.openpi;

  /*** occupied-virtual relaxation terms */

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

  /* I(i,a) = I'(i,a) - sum_m f(i,m) D(orb)(a,m) */
  dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I'ia");
  dpd_file2_copy(&I, CC_OEI, "I(i,a)");
  dpd_file2_close(&I);
  dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I(i,a)");
  dpd_file2_init(&D, CC_OEI, 0, 1, 0, "D(orb)(a,i)");
  dpd_file2_init(&f, CC_OEI, 0, 0, 0, "fij");
  dpd_contract222(&f, &D, &I, 0, 0, -1.0, 1.0);
  dpd_file2_close(&f);
  dpd_file2_close(&D);
  dpd_file2_close(&I);

  /*** occupied-occupied relaxtion terms */

  /* I(I,J) <-- I'(I,J) - sum_E,M D(orb)(E,M) [<EI||MJ> + <EJ||MI>]
                      - 2 sum_e,m D(orb)(e,m) <eI|mJ> */
  dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'IJ");
  dpd_file2_copy(&I, CC_OEI, "I(I,J)");
  dpd_file2_close(&I);
  dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I(I,J)");
  dpd_file2_init(&D, CC_OEI, 0, 1, 0, "D(orb)(A,I)");
  dpd_buf4_init(&E, CC_EINTS, 0, 11, 0, 11, 0, 1, "E <ai|jk>");
  dpd_dot13(&D, &E, &I, 0, 0, -1.0, 1.0);
  dpd_dot13(&D, &E, &I, 0, 1, -1.0, 1.0);
  dpd_buf4_close(&E);
  dpd_file2_close(&D);
  dpd_file2_init(&D, CC_OEI, 0, 1, 0, "D(orb)(a,i)");
  dpd_buf4_init(&E, CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
  dpd_dot13(&D, &E, &I, 0, 0, -2.0, 1.0);
  dpd_buf4_close(&E);
  dpd_file2_close(&D);

  /* I(I,J) <-- - 2 sum_E  f(E,I) D(orb)(E,J) (J,socc)

   Note that this same term is not needed in the I(i,j) block since J
   is required to be a singly occupied orbital */
  dpd_file2_mat_init(&I);
  dpd_file2_mat_rd(&I);

  dpd_file2_init(&f, CC_OEI, 0, 0, 1, "fIA");
  dpd_file2_mat_init(&f);
  dpd_file2_mat_rd(&f);
  dpd_file2_init(&D, CC_OEI, 0, 1, 0, "D(orb)(A,I)");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);

  for(h=0; h < nirreps; h++) {
      for(i=0; i < occpi[h]; i++)
	  for(j=(occpi[h] - openpi[h]); j < occpi[h]; j++)
	      for(e=0; e < virtpi[h]; e++)
		  I.matrix[h][i][j] -= 2 * f.matrix[h][i][e] * D.matrix[h][e][j];
    }

  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);
  dpd_file2_mat_close(&f);
  dpd_file2_close(&f);

  dpd_file2_mat_wrt(&I);
  dpd_file2_mat_close(&I);
  dpd_file2_close(&I);

  /* I(i,j) <-- I'(i,j) - sum_e,m D(orb)(e,m) [<ei||mj> + <ej||mi>]
                      - 2 sum_E,M D(orb)(E,M) <Ei|Mj> */
  dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'ij");
  dpd_file2_copy(&I, CC_OEI, "I(i,j)");
  dpd_file2_close(&I);
  dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I(i,j)");
  dpd_file2_init(&D, CC_OEI, 0, 1, 0, "D(orb)(a,i)");
  dpd_buf4_init(&E, CC_EINTS, 0, 11, 0, 11, 0, 1, "E <ai|jk>");
  dpd_dot13(&D, &E, &I, 0, 0, -1.0, 1.0);
  dpd_dot13(&D, &E, &I, 0, 1, -1.0, 1.0);
  dpd_buf4_close(&E);
  dpd_file2_close(&D);
  dpd_file2_init(&D, CC_OEI, 0, 1, 0, "D(orb)(A,I)");
  dpd_buf4_init(&E, CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
  dpd_dot13(&D, &E, &I, 0, 0, -2.0, 1.0);
  dpd_buf4_close(&E);
  dpd_file2_close(&D);
  dpd_file2_close(&I);

  /* Clean the I(i,j) block yet again */
  dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I(i,j)");
  dpd_file2_mat_init(&I);
  dpd_file2_mat_rd(&I);

  for(h=0; h < nirreps; h++)
      for(i=0; i < occpi[h]; i++)
	  for(j=0; j < occpi[h]; j++)
	      if((i >= (occpi[h] - openpi[h])) ||
		 (j >= (occpi[h] - openpi[h])) )
		  I.matrix[h][i][j] = 0.0;

  dpd_file2_mat_wrt(&I);
  dpd_file2_mat_close(&I);
  dpd_file2_close(&I);
}

}} // namespace psi::ccdensity
