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
    \ingroup CCLAMBDA
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cclambda {

/* spinad_amps(): For RHF references, build the T2 AA and BB amplitudes from 
** the existing T2 AB amplitudes and copy the existing T1 A amplitudes 
** into B.
**
** T2(IJ,AB) = T2(ij,ab) = T2(Ij,Ab) - T2(Ij,Ba)
**
** T1(I,A) = T1(i,a)
**
*/

void spinad_amps(void)
{
  dpdfile2 T1;
  dpdbuf4 T2;

  if(params.ref == 0) { /** RHF **/

    dpd_file2_init(&T1, PSIF_CC_LAMBDA, 0, 0, 1, "LIA");
    dpd_file2_copy(&T1, PSIF_CC_LAMBDA, "Lia");
    dpd_file2_close(&T1);

    dpd_buf4_init(&T2, PSIF_CC_LAMBDA, 0, 2, 7, 0, 5, 1, "LIjAb");
    dpd_buf4_copy(&T2, PSIF_CC_LAMBDA, "LIJAB");
    dpd_buf4_copy(&T2, PSIF_CC_LAMBDA, "Lijab");
    dpd_buf4_close(&T2);

  }
}

}} // namespace psi::cclambda
