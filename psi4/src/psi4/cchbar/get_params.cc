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
    \ingroup CCHBAR
    \brief Enter brief description of file here
*/

#include <cstdio>
#include <cstdlib>
#include <string>
#include <cmath>
#include "psi4/liboptions/liboptions.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/psi4-dec.h"
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cchbar {

void get_params(Options &options)
{
  params.memory = Process::environment.get_memory();

  /* compute the Tamplitude equation matrix elements (usually 0) */
//  params.Tamplitude = 0;
//  errcod = ip_boolean("TAMPLITUDE", &(params.Tamplitude),0);
  params.Tamplitude = options.get_bool("T_AMPS");

//  params.cachelev = 2;
//  errcod = ip_data("CACHELEVEL", "%d", &(params.cachelev),0);
  params.cachelev = options.get_int("CACHELEVEL");

//  params.print = 0;
//  errcod = ip_data("PRINT", "%d", &(params.print),0);
  params.print = options.get_int("PRINT");

//  errcod = ip_string("WFN", &(params.wfn), 0);
  params.wfn = options.get_str("WFN");

//  params.dertype = 0;
  std::string junk = options.get_str("DERTYPE");
  if(junk == "NONE") params.dertype = 0;
  else if(junk == "FIRST") params.dertype = 1;
  else if(junk == "RESPONSE") params.dertype = 3; /* linear response */
  else {
//    printf("Invalid value of input keyword DERTYPE: %s\n", junk);
//    return PSI_RETURN_FAILURE;
      throw PsiException("CCHBAR: Invalid value of input keyword DERTYPE",__FILE__,__LINE__);
  }

  /* Should we use the minimal-disk algorithm for Wabei?  It's VERY slow! */
//  params.wabei_lowdisk = 0;
//  errcod = ip_boolean("WABEI_LOWDISK", &params.wabei_lowdisk, 0);
  params.wabei_lowdisk = options.get_bool("WABEI_LOWDISK");
}

}} // namespace psi::cchbar
