/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
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
 * @END LICENSE
 */

/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {
#include "psi4/physconst.h"

#define _hartree2nm 0.02194746313710

void ex_td_print(std::vector<struct XTD_Params> xtd_list)
{
  int i;

  outfile->Printf("\n\t                   Excited State -> Excited State Transitions\n");
  outfile->Printf("\n\t                        Excitation Energy          OS       RS        RS     Einstein A\n");
  outfile->Printf("\tTransition   (eV)    (cm^-1)    (nm)     (au)              (l,au)   (v,au)     (s^-1)\n");
  for(i=0; i<xtd_list.size(); i++) {
    //outfile->Printf("\t  %d%s->%d%s %7.3lf %9.1lf %7.1lf %10.6lf %8.4lf %8.4lf %8.4lf  %12.1lf\n",
    outfile->Printf("\t  %d%s->%d%s %7.3lf %9.1lf %7.1lf %10.6lf %8.4lf %8.4lf %8.4lf  %7.6E\n",
            xtd_list[i].root1+1,moinfo.labels[xtd_list[i].irrep1],
            xtd_list[i].root2+1,moinfo.labels[xtd_list[i].irrep2],
            xtd_list[i].cceom_energy*pc_hartree2ev,
            xtd_list[i].cceom_energy*pc_hartree2wavenumbers,
            1/(xtd_list[i].cceom_energy*_hartree2nm),
            xtd_list[i].cceom_energy,xtd_list[i].OS,
            xtd_list[i].RS_length,xtd_list[i].RS_velocity,
            xtd_list[i].einstein_a);
  }
  outfile->Printf("\n");
}

}} // namespace psi::ccdensity
