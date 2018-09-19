/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2018 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
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

void td_print(void)
{
  int i;

  outfile->Printf("\n\t                   Ground State -> Excited State Transitions\n");
  outfile->Printf("\n\t                   Excitation Energy          OS       RS        RS     Einstein A\n");
  outfile->Printf("\tState   (eV)    (cm^-1)    (nm)     (au)              (l,au)   (v,au)     (s^-1)\n");
  for(i=0; i<params.nstates; i++) {
    //outfile->Printf("\t %d%3s %7.3lf %9.1lf %7.1lf %10.6lf %8.4lf %8.4lf %8.4lf  %12.1lf\n",
    outfile->Printf("\t %d%3s %7.3lf %9.1lf %7.1lf %10.6lf %8.4lf %8.4lf %8.4lf  %7.6E\n",
            td_params[i].root+1,moinfo.labels[td_params[i].irrep].c_str(),
            td_params[i].cceom_energy*pc_hartree2ev,
            td_params[i].cceom_energy*pc_hartree2wavenumbers,
            1/(td_params[i].cceom_energy*_hartree2nm),
            td_params[i].cceom_energy, td_params[i].OS,
            td_params[i].RS_length,td_params[i].RS_velocity,
            td_params[i].einstein_a);
  }
  outfile->Printf("\n");
}

}} // namespace psi::ccdensity
