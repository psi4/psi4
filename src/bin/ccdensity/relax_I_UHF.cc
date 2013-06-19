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
  dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 1, "I'IA");
  dpd_->file2_copy(&I, PSIF_CC_OEI, "I(I,A)");
  dpd_->file2_close(&I);

  dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 1, "I(I,A)");
  dpd_->file2_mat_init(&I);
  dpd_->file2_mat_rd(&I);

  dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 0, "D(orb)(A,I)");
  dpd_->file2_mat_init(&D);
  dpd_->file2_mat_rd(&D);

  dpd_->file2_init(&f, PSIF_CC_OEI, 0, 0, 0, "fIJ");
  dpd_->file2_mat_init(&f);
  dpd_->file2_mat_rd(&f);

  for(h=0; h < nirreps; h++)
    for(i=0; i < aoccpi[h]; i++)
      for(a=0; a < avirtpi[h]; a++)
	I.matrix[h][i][a] -= D.matrix[h][a][i] * f.matrix[h][i][i];

  dpd_->file2_mat_close(&f);
  dpd_->file2_close(&f);
  dpd_->file2_mat_close(&D);
  dpd_->file2_close(&D);

  dpd_->file2_mat_wrt(&I);
  dpd_->file2_mat_close(&I);
  dpd_->file2_close(&I);

  /* I(i,a) = I'(i,a) - f(i,i) D(orb)(a,i) */
  dpd_->file2_init(&I, PSIF_CC_OEI, 0, 2, 3, "I'ia");
  dpd_->file2_copy(&I, PSIF_CC_OEI, "I(i,a)");
  dpd_->file2_close(&I);

  dpd_->file2_init(&I, PSIF_CC_OEI, 0, 2, 3, "I(i,a)");
  dpd_->file2_mat_init(&I);
  dpd_->file2_mat_rd(&I);

  dpd_->file2_init(&D, PSIF_CC_OEI, 0, 3, 2, "D(orb)(a,i)");
  dpd_->file2_mat_init(&D);
  dpd_->file2_mat_rd(&D);

  dpd_->file2_init(&f, PSIF_CC_OEI, 0, 2, 2, "fij");
  dpd_->file2_mat_init(&f);
  dpd_->file2_mat_rd(&f);

  for(h=0; h < nirreps; h++)
    for(i=0; i < boccpi[h]; i++)
      for(a=0; a < bvirtpi[h]; a++)
	I.matrix[h][i][a] -= D.matrix[h][a][i] * f.matrix[h][i][i];

  dpd_->file2_mat_close(&f);
  dpd_->file2_close(&f);
  dpd_->file2_mat_close(&D);
  dpd_->file2_close(&D);

  dpd_->file2_mat_wrt(&I);
  dpd_->file2_mat_close(&I);
  dpd_->file2_close(&I);

  /*** occupied-occupied relaxtion terms */

  /* I(I,J) <-- I'(I,J) - sum_A,K D(orb)(A,K) [<AI||KJ> + <AJ||KI>] - 2 sum_a,k D(orb)(a,k) <aI|kJ> */
  dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'IJ");
  dpd_->file2_copy(&I, PSIF_CC_OEI, "I(I,J)");
  dpd_->file2_close(&I);
  dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I(I,J)");
  dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 0, "D(orb)(A,I)");
  dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 21, 0, 21, 0, 1, "E <AI|JK>");
  dpd_->dot13(&D, &E, &I, 0, 0, -1.0, 1.0);
  dpd_->dot13(&D, &E, &I, 0, 1, -1.0, 1.0);
  dpd_->buf4_close(&E);
  dpd_->file2_close(&D);
  dpd_->file2_init(&D, PSIF_CC_OEI, 0, 3, 2, "D(orb)(a,i)");
  dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 24, 22, 24, 22, 0, "E <Ia|Jk>");
  dpd_->dot24(&D, &E, &I, 0, 0, -2.0, 1.0);
  dpd_->buf4_close(&E);
  dpd_->file2_close(&D);
  dpd_->file2_close(&I);

  /* I(i,j) <-- I'(i,j) - sum_a,k D(orb)(a,k) [<ai||kj> + <aj||ki>] - 2 sum_A,K D(orb)(A,K) <Ai|Kj> */
  dpd_->file2_init(&I, PSIF_CC_OEI, 0, 2, 2, "I'ij");
  dpd_->file2_copy(&I, PSIF_CC_OEI, "I(i,j)");
  dpd_->file2_close(&I);
  dpd_->file2_init(&I, PSIF_CC_OEI, 0, 2, 2, "I(i,j)");
  dpd_->file2_init(&D, PSIF_CC_OEI, 0, 3, 2, "D(orb)(a,i)");
  dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 31, 10, 31, 10, 1, "E <ai|jk>");
  dpd_->dot13(&D, &E, &I, 0, 0, -1.0, 1.0);
  dpd_->dot13(&D, &E, &I, 0, 1, -1.0, 1.0);
  dpd_->buf4_close(&E);
  dpd_->file2_close(&D);
  dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 0, "D(orb)(A,I)");
  dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 26, 22, 26, 22, 0, "E <Ai|Jk>");
  dpd_->dot13(&D, &E, &I, 0, 0, -2.0, 1.0);
  dpd_->buf4_close(&E);
  dpd_->file2_close(&D);
  dpd_->file2_close(&I);

}

}} // namespace psi::ccdensity
