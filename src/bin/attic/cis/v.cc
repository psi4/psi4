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
#include <cstdio>
#include <cstring>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cis {

void v_build(int irrep, int root, enum Spin spin)
{
  char lbl[32];
  dpdfile2 B, V, B_A, B_B, V_A, V_B, F;
  dpdbuf4 T2;

  if(params.ref == 0) { /** RHF **/

    if(spin == singlet)
      sprintf(lbl, "BIA(%d)[%d] singlet", root, irrep);
    else
      sprintf(lbl, "BIA(%d)[%d] triplet", root, irrep);
    dpd_file2_init(&B, PSIF_CC_OEI, irrep, 0, 1, lbl);

    sprintf(lbl, "VIA[%d]", irrep);
    dpd_file2_init(&V, PSIF_CC_MISC, irrep, 0, 1, lbl);

    dpd_file2_init(&F, PSIF_CC_MISC, 0, 1, 1, "FAB");
    dpd_contract222(&B, &F, &V, 0, 0, 1, 0);
    dpd_file2_close(&F);

    dpd_file2_init(&F, PSIF_CC_MISC, 0, 0, 0, "FIJ");
    dpd_contract222(&F, &B, &V, 0, 1, 1, 1);
    dpd_file2_close(&F);

    sprintf(lbl, "FKC(%d)[%d]", root, irrep);
    dpd_file2_init(&F, PSIF_CC_MISC, irrep, 0, 1, lbl);

    if(spin == singlet) {
      dpd_buf4_init(&T2, PSIF_CC_MISC, 0, 0, 5, 0, 5, 0, "MP2 2 tIjAb - tIjbA");
      dpd_dot24(&F, &T2, &V, 0, 0, 1, 1);
    }
    else {
      dpd_buf4_init(&T2, PSIF_CC_MISC, 0, 0, 5, 0, 5, 0, "MP2 tIjAb");
      dpd_dot23(&F, &T2, &V, 0, 0, -1, 1);
    }

    dpd_buf4_close(&T2);
    dpd_file2_close(&F);

    dpd_file2_close(&V);
    dpd_file2_close(&B);

  }
  else if(params.ref == 2) { /** UHF **/

    sprintf(lbl, "BIA(%d)[%d]", root, irrep);
    dpd_file2_init(&B_A, PSIF_CC_OEI, irrep, 0, 1, lbl);

    sprintf(lbl, "Bia(%d)[%d]", root, irrep);
    dpd_file2_init(&B_B, PSIF_CC_OEI, irrep, 2, 3, lbl);

    sprintf(lbl, "VIA[%d]", irrep);
    dpd_file2_init(&V_A, PSIF_CC_MISC, irrep, 0, 1, lbl);
    sprintf(lbl, "Via[%d]", irrep);
    dpd_file2_init(&V_B, PSIF_CC_MISC, irrep, 2, 3, lbl);

    dpd_file2_init(&F, PSIF_CC_MISC, 0, 1, 1, "FAB");
    dpd_contract222(&B_A, &F, &V_A, 0, 0, 1, 0);
    dpd_file2_close(&F);

    dpd_file2_init(&F, PSIF_CC_MISC, 0, 0, 0, "FIJ");
    dpd_contract222(&F, &B_A, &V_A, 0, 1, 1, 1);
    dpd_file2_close(&F);

    sprintf(lbl, "FKC(%d)[%d]", root, irrep);
    dpd_file2_init(&F, PSIF_CC_MISC, irrep, 0, 1, lbl);
    dpd_buf4_init(&T2, PSIF_CC_MISC, 0, 0, 5, 2, 7, 0, "MP2 tIJAB");
    dpd_dot24(&F, &T2, &V_A, 0, 0, 1, 1);
    dpd_buf4_close(&T2);
    dpd_file2_close(&F);

    sprintf(lbl, "Fkc(%d)[%d]", root, irrep);
    dpd_file2_init(&F, PSIF_CC_MISC, irrep, 2, 3, lbl);
    dpd_buf4_init(&T2, PSIF_CC_MISC, 0, 22, 28, 22, 28, 0, "MP2 tIjAb");
    dpd_dot24(&F, &T2, &V_A, 0, 0, 1, 1);
    dpd_buf4_close(&T2);
    dpd_file2_close(&F);

    dpd_file2_init(&F, PSIF_CC_MISC, 0, 3, 3, "Fab");
    dpd_contract222(&B_B, &F, &V_B, 0, 0, 1, 0);
    dpd_file2_close(&F);

    dpd_file2_init(&F, PSIF_CC_MISC, 0, 2, 2, "Fij");
    dpd_contract222(&F, &B_B, &V_B, 0, 1, 1, 1);
    dpd_file2_close(&F);

    sprintf(lbl, "Fkc(%d)[%d]", root, irrep);
    dpd_file2_init(&F, PSIF_CC_MISC, irrep, 2, 3, lbl);
    dpd_buf4_init(&T2, PSIF_CC_MISC, 0, 10, 15, 12, 17, 0, "MP2 tijab");
    dpd_dot24(&F, &T2, &V_B, 0, 0, 1, 1);
    dpd_buf4_close(&T2);
    dpd_file2_close(&F);

    sprintf(lbl, "FKC(%d)[%d]", root, irrep);
    dpd_file2_init(&F, PSIF_CC_MISC, irrep, 0, 1, lbl);
    dpd_buf4_init(&T2, PSIF_CC_MISC, 0, 22, 28, 22, 28, 0, "MP2 tIjAb");
    dpd_dot13(&F, &T2, &V_B, 0, 0, 1, 1);
    dpd_buf4_close(&T2);
    dpd_file2_close(&F);

    dpd_file2_close(&V_A);
    dpd_file2_close(&V_B);

    dpd_file2_close(&B_A);
    dpd_file2_close(&B_B);

  }

}

}} // namespace psi::cis
