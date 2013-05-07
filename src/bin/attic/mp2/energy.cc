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
    \ingroup MP2
    \brief Enter brief description of file here
*/
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"
#include "moinfo.h"

namespace psi{ namespace mp2{

double rhf_energy(void);
double uhf_energy(void);

double energy(void)
{
  double e = 0.0;
  if(params.ref == 0) e = rhf_energy();
  else if(params.ref == 2) e = uhf_energy();
  return e;
}

double rhf_energy(void)
{
  double E = 0;
  dpdbuf4 tIjAb;
  dpdbuf4 D;
  dpdbuf4 S;
  double os_energy, ss_energy;

  dpd_buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
  E = dpd_buf4_dot(&D, &tIjAb);
  dpd_buf4_close(&D);
  dpd_buf4_init(&S, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  os_energy = dpd_buf4_dot(&S, &tIjAb);
  dpd_buf4_close(&tIjAb);
  dpd_buf4_close(&S);
  ss_energy = (E - os_energy);

  mo.emp2_os = os_energy;
  mo.emp2_ss = ss_energy;

  if (params.scs == 1) {
    os_energy = params.scs_scale_os * os_energy;
    ss_energy = params.scs_scale_ss * ss_energy;
  }
  else {
    os_energy = (6.0/5.6) * os_energy;
    ss_energy = (1.0/3.0) * ss_energy;
  }

  mo.escsmp2_os = os_energy;
  mo.escsmp2_ss = ss_energy;

  return(E);
}

double uhf_energy(void)
{
  double E1A = 0; double E1B = 0;
  double E2AA = 0; double E2BB = 0; double E2AB = 0;
  dpdfile2 T1, F;
  dpdbuf4 tIJAB, tijab, tIjAb,  D;

  if(params.semicanonical) {
    dpd_file2_init(&F, PSIF_CC_OEI, 0, 0, 1, "fIA");
    dpd_file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    E1A = dpd_file2_dot(&F, &T1);
    dpd_file2_close(&F);
    dpd_file2_close(&T1);

    dpd_file2_init(&F, PSIF_CC_OEI, 0, 2, 3, "fia");
    dpd_file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
    E1B = dpd_file2_dot(&F, &T1);
    dpd_file2_close(&F);
    dpd_file2_close(&T1);
  }

  dpd_buf4_init(&tIJAB, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tIJAB");
  dpd_buf4_init(&D, PSIF_CC_DINTS, 0, 2, 7, 2, 7, 0, "D <IJ||AB> (I>J,A>B)");
  E2AA = dpd_buf4_dot(&D, &tIJAB);
  dpd_buf4_close(&D);
  dpd_buf4_close(&tIJAB);

  dpd_buf4_init(&tijab, PSIF_CC_TAMPS, 0, 12, 17, 12, 17, 0, "tijab");
  dpd_buf4_init(&D, PSIF_CC_DINTS, 0, 12, 17, 12, 17, 0, "D <ij||ab> (i>j,a>b)");
  E2BB = dpd_buf4_dot(&D, &tijab);
  dpd_buf4_close(&D);
  dpd_buf4_close(&tijab);

  dpd_buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
  dpd_buf4_init(&D, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
  E2AB = dpd_buf4_dot(&D, &tIjAb);
  dpd_buf4_close(&D);
  dpd_buf4_close(&tIjAb);

  if(params.semicanonical)
    return(E1A+E1B+E2AA+E2BB+E2AB);
  else
    return(E2AA+E2BB+E2AB);
}

}} /* End namespaces */
