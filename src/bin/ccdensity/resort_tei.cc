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
#include <cstdio>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

void distribute(void);

int file_build(dpdfile4 *File, int inputfile, double tolerance,
	       int perm_pr, int perm_qs, int perm_prqs, int keep);

void resort_tei(void)
{
  double tolerance;
  dpdfile4 A, B, C, D, E, F;

  tolerance = params.tolerance;

  distribute();

  dpd_->file4_init_nocache(&A, PSIF_CC_AINTS_NEW, 0, 0, 0, "A <ij|kl>");
  file_build(&A, 90, tolerance, 1, 1, 1, 0);
  dpd_->file4_close(&A);

  dpd_->file4_init_nocache(&B, PSIF_CC_BINTS_NEW, 0, 5, 5, "B <ab|cd>");
  file_build(&B, 91, tolerance, 1, 1, 1, 0);
  dpd_->file4_close(&B);

  dpd_->file4_init_nocache(&C, PSIF_CC_CINTS_NEW, 0, 10, 10, "C <ia|jb>");
  file_build(&C, 92, tolerance, 1, 1, 0, 0);
  dpd_->file4_close(&C);

  dpd_->file4_init_nocache(&D, PSIF_CC_DINTS_NEW, 0, 0, 5, "D <ij|ab>");
  file_build(&D, 93, tolerance, 0, 0, 1, 0);
  dpd_->file4_close(&D);

  dpd_->file4_init_nocache(&E, PSIF_CC_EINTS_NEW, 0, 11, 0, "E <ai|jk>");
  file_build(&E, 94, tolerance, 0, 1, 0, 0);
  dpd_->file4_close(&E);

  dpd_->file4_init_nocache(&F, PSIF_CC_FINTS_NEW, 0, 10, 5, "F <ia|bc>");
  file_build(&F, 95, tolerance, 0, 1, 0, 0);
  dpd_->file4_close(&F);

  fflush(outfile);

}

}} // namespace psi::ccdensity
