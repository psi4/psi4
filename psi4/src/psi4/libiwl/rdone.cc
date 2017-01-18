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

/*!
  \file
  \ingroup IWL
*/
#include <cstdio>
#include <cmath>
#include "psi4/libpsio/psio.h"
#include "psi4/libciomr/libciomr.h"
#include "iwl.h"
#include "iwl.hpp"
#include "psi4/libparallel/ParallelPrinter.h"
namespace psi {

void IWL::read_one(PSIO *psio, int itap, const char *label, double *ints,
    int ntri, int erase, int printflg, std::string out)
{
   std::shared_ptr<psi::PsiOutStream> printer=(out=="outfile"?outfile:
         std::shared_ptr<OutFile>(new OutFile(out)));
    int nmo;

    psio->open(itap, PSIO_OPEN_OLD);
    psio->read_entry(itap, label, (char *) ints, ntri*sizeof(double));
    psio->close(itap, !erase);

    if (printflg) {
        nmo = (int) (sqrt((double) (1 + 8 * ntri)) - 1)/2;
        print_array(ints, nmo, out);
    }
}

/*!
** IWL_RDONE()
**
** This function reads the one-electron integrals in the MO basis
**  from disk and stores them in core.  Substantially revised on
**  29 April 1998 to filter out frozen orbitals if requested.
**  This change requires a very different argument list from the
**  previous version of this code.
**
** David Sherrill, January 1994
** Revised by David Sherrill, April 1998
** Revised by TDC, June 2001
**
**   \param itap       = tape to read ints from
**   \param label      = the PSIO label
**   \param ints       = buffer (already allocated) to store the integrals
**   \param ntri       = number of packed integrals
**   \param erase      = erase itap (1=yes, 0=no)
**   \param printflg   = printing flag.  Set to 1 to print ints;
**                       otherwise, set to 0
**   \param out    = file pointer for output of ints or error messages
** \ingroup IWL
*/
int iwl_rdone(int itap, const char *label, double *ints, int ntri, int erase,
              int printflg,std::string out)
{
   std::shared_ptr<psi::PsiOutStream> printer=(out=="outfile"?outfile:
         std::shared_ptr<OutFile>(new OutFile(out)));
  int nmo;

  psio_open(itap, PSIO_OPEN_OLD);
  psio_read_entry(itap, label, (char *) ints, ntri*sizeof(double));
  psio_close(itap, !erase);

  if (printflg) {
    nmo = (int) (sqrt((double) (1 + 8 * ntri)) - 1)/2;
    print_array(ints, nmo, out);
  }

  return(1);
}

}
