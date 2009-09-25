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

/* relax_I_UHF(): Add the UHF orbital-response contributions from the
** one-electron density matrix to the I(I,J) and I(I,A) blocks of the
** Lagrangian.  These terms arise from the first-order CPHF equations.
*/

void relax_I_UHF(void)
{
  dpdfile2 I, D, f;
  dpdbuf4 E;
  int h, nirreps, i, a;
  int *aoccpi, *avirtpi;
  int *boccpi, *bvirtpi;

  nirreps = moinfo.nirreps;
  aoccpi = moinfo.aoccpi;  avirtpi = moinfo.avirtpi;
  boccpi = moinfo.boccpi;  bvirtpi = moinfo.bvirtpi;

  /*** occupied-virtual relaxation terms */

  /* I(I,A) = I'(I,A) - f(I,I) D(orb)(A,I) */
  dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I'IA");
  dpd_file2_copy(&I, CC_OEI, "I(I,A)");
  dpd_file2_close(&I);

  dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I(I,A)");
  dpd_file2_mat_init(&I);
  dpd_file2_mat_rd(&I);

  dpd_file2_init(&D, CC_OEI, 0, 1, 0, "D(orb)(A,I)");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);

  dpd_file2_init(&f, CC_OEI, 0, 0, 0, "fIJ");
  dpd_file2_mat_init(&f);
  dpd_file2_mat_rd(&f);

  for(h=0; h < nirreps; h++)
    for(i=0; i < aoccpi[h]; i++)
      for(a=0; a < avirtpi[h]; a++)
	I.matrix[h][i][a] -= D.matrix[h][a][i] * f.matrix[h][i][i];

  dpd_file2_mat_close(&f);
  dpd_file2_close(&f);
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  dpd_file2_mat_wrt(&I);
  dpd_file2_mat_close(&I);
  dpd_file2_close(&I);

  /* I(i,a) = I'(i,a) - f(i,i) D(orb)(a,i) */
  dpd_file2_init(&I, CC_OEI, 0, 2, 3, "I'ia");
  dpd_file2_copy(&I, CC_OEI, "I(i,a)");
  dpd_file2_close(&I);

  dpd_file2_init(&I, CC_OEI, 0, 2, 3, "I(i,a)");
  dpd_file2_mat_init(&I);
  dpd_file2_mat_rd(&I);

  dpd_file2_init(&D, CC_OEI, 0, 3, 2, "D(orb)(a,i)");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);

  dpd_file2_init(&f, CC_OEI, 0, 2, 2, "fij");
  dpd_file2_mat_init(&f);
  dpd_file2_mat_rd(&f);

  for(h=0; h < nirreps; h++)
    for(i=0; i < boccpi[h]; i++)
      for(a=0; a < bvirtpi[h]; a++)
	I.matrix[h][i][a] -= D.matrix[h][a][i] * f.matrix[h][i][i];

  dpd_file2_mat_close(&f);
  dpd_file2_close(&f);
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  dpd_file2_mat_wrt(&I);
  dpd_file2_mat_close(&I);
  dpd_file2_close(&I);

  /*** occupied-occupied relaxtion terms */

  /* I(I,J) <-- I'(I,J) - sum_A,K D(orb)(A,K) [<AI||KJ> + <AJ||KI>] - 2 sum_a,k D(orb)(a,k) <aI|kJ> */
  dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'IJ");
  dpd_file2_copy(&I, CC_OEI, "I(I,J)");
  dpd_file2_close(&I);
  dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I(I,J)");
  dpd_file2_init(&D, CC_OEI, 0, 1, 0, "D(orb)(A,I)");
  dpd_buf4_init(&E, CC_EINTS, 0, 21, 0, 21, 0, 1, "E <AI|JK>");
  dpd_dot13(&D, &E, &I, 0, 0, -1.0, 1.0);
  dpd_dot13(&D, &E, &I, 0, 1, -1.0, 1.0);
  dpd_buf4_close(&E);
  dpd_file2_close(&D);
  dpd_file2_init(&D, CC_OEI, 0, 3, 2, "D(orb)(a,i)");
  dpd_buf4_init(&E, CC_EINTS, 0, 24, 22, 24, 22, 0, "E <Ia|Jk>");
  dpd_dot24(&D, &E, &I, 0, 0, -2.0, 1.0);
  dpd_buf4_close(&E);
  dpd_file2_close(&D);
  dpd_file2_close(&I);

  /* I(i,j) <-- I'(i,j) - sum_a,k D(orb)(a,k) [<ai||kj> + <aj||ki>] - 2 sum_A,K D(orb)(A,K) <Ai|Kj> */
  dpd_file2_init(&I, CC_OEI, 0, 2, 2, "I'ij");
  dpd_file2_copy(&I, CC_OEI, "I(i,j)");
  dpd_file2_close(&I);
  dpd_file2_init(&I, CC_OEI, 0, 2, 2, "I(i,j)");
  dpd_file2_init(&D, CC_OEI, 0, 3, 2, "D(orb)(a,i)");
  dpd_buf4_init(&E, CC_EINTS, 0, 31, 10, 31, 10, 1, "E <ai|jk>");
  dpd_dot13(&D, &E, &I, 0, 0, -1.0, 1.0);
  dpd_dot13(&D, &E, &I, 0, 1, -1.0, 1.0);
  dpd_buf4_close(&E);
  dpd_file2_close(&D);
  dpd_file2_init(&D, CC_OEI, 0, 1, 0, "D(orb)(A,I)");
  dpd_buf4_init(&E, CC_EINTS, 0, 26, 22, 26, 22, 0, "E <Ai|Jk>");
  dpd_dot13(&D, &E, &I, 0, 0, -2.0, 1.0);
  dpd_buf4_close(&E);
  dpd_file2_close(&D);
  dpd_file2_close(&I);

}

}} // namespace psi::ccdensity
