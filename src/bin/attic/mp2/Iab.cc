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

/*! 
** \file
** \ingroup MP2
** \brief Enter brief description of file here 
*/

#include <cmath>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

namespace psi{ namespace mp2{

void rhf_sf_Iab(void);
void uhf_sf_Iab(void);

void Iab(void)
{
  if(params.ref == 0) rhf_sf_Iab();
  else if(params.ref == 2) uhf_sf_Iab();
}

void rhf_sf_Iab(void)
{
  dpdfile2 F, D, I;
  dpdbuf4 G, Dints;

  /* I'AB <-- sum_I fAI (DBI + DIB) + sum_C fAC (DBC + DCB) */
  global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 1, "I'AB");

  global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 1, 1, "fAB");
  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, "DAB");
  global_dpd_->contract222(&F, &D, &I, 0, 0, 1.0, 0.0);
  global_dpd_->contract222(&F, &D, &I, 0, 1, 1.0, 1.0);
  global_dpd_->file2_close(&D);
  global_dpd_->file2_close(&F);

  global_dpd_->file2_close(&I);

  /* I'ab <-- sum_i fai (Dbi + Dib) + sum_c fac (Dbc + Dcb) */
  global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 1, "I'ab");

  global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 1, 1, "fab");
  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, "Dab");
  global_dpd_->contract222(&F, &D, &I, 0, 0, 1.0, 0.0);
  global_dpd_->contract222(&F, &D, &I, 0, 1, 1.0, 1.0);
  global_dpd_->file2_close(&D);
  global_dpd_->file2_close(&F);

  global_dpd_->file2_close(&I);

  /* I'AB <-- sum_CIJ <IJ||CA> G(IJ,CB) + 2 sum_Ijc <Ij|Ac> G(Ij,Bc) */
  global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 1, "I'AB");

  global_dpd_->buf4_init(&Dints, PSIF_CC_DINTS, 0, 2, 5, 2, 5, 0, "D <ij||ab> (i>j,ab)");
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 2, 5, 2, 7, 0, "GIJAB");
  global_dpd_->contract442(&Dints, &G, &I, 3, 3, 2.0, 1.0);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&Dints);

  global_dpd_->buf4_init(&Dints, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
  global_dpd_->contract442(&Dints, &G, &I, 2, 2, 2.0, 1.0);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&Dints);

  global_dpd_->file2_close(&I);

  /* I'ab <-- sum_cij <ij||ca> G(ij,cb) + 2 sum_IjC <Ij|Ca> G(Ij,Cb) */
  global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 1, "I'ab");

  global_dpd_->buf4_init(&Dints, PSIF_CC_DINTS, 0, 2, 5, 2, 5, 0, "D <ij||ab> (i>j,ab)");
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 2, 5, 2, 7, 0, "Gijab");
  global_dpd_->contract442(&Dints, &G, &I, 3, 3, 2.0, 1.0);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&Dints);

  global_dpd_->buf4_init(&Dints, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
  global_dpd_->contract442(&Dints, &G, &I, 3, 3, 2.0, 1.0);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&Dints);

  global_dpd_->file2_close(&I);
}

void uhf_sf_Iab(void)
{

}


}} /* End namespaces */
