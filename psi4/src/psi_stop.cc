/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
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

/*!
** \file
** \brief Close input and output, stop input parser
** \ingroup
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "psi4/psifiles.h"
#include "psi4/psi4-dec.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libplugin/plugin.h"
#include "psi4/libparallel/ParallelPrinter.h"
namespace psi {

/*!
** psi_stop(): This function closes input and output files and
** deinitializes Input Parsing library.
**
** Arguments: none
**
** Returns: one of standard PSI error codes
** \ingroup CIOMR
*/

int psi_stop(FILE* infile, std::string OutFileRMR, char* psi_file_prefix)
{
  free(psi_file_prefix);

  // Success Flag, so a user can tell via grep that the outfile worked (or at least didn't segfault)
  // With a little Psi4 flavor to it.
  outfile->Printf( "\n*** Psi4 exiting successfully. Buy a developer a beer!\n");


  //if (outfile)
  //    fclose(outfile);
  if (infile)
      fclose(infile);

  infile = NULL;
  outfile = boost::shared_ptr<OutFile>();
  psi_file_prefix = NULL;

  //psi::yetiEnv.free();

  return(PSI_RETURN_SUCCESS);
}

}
