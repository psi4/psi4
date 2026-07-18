/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2026 The Psi4 Developers.
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

/*!
  \file
  \ingroup IWL
*/
#include <cmath>
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libciomr/libciomr.h"
#include "iwl.h"
#include "iwl.hpp"

namespace psi {

namespace {

void read_one_impl(PSIO *psio, int itap, const char *label, double *ints, int ntri, int erase, int printflg,
                   const std::string &out) {
    psio->open(itap, PSIO_OPEN_OLD);
    psio->read_entry(itap, label, reinterpret_cast<char *>(ints), ntri * sizeof(double));
    psio->close(itap, erase ? 0 : 1);

    if (printflg) {
        const int nmo = static_cast<int>(std::sqrt(static_cast<double>(1 + 8 * ntri)) - 1) / 2;
        print_array(ints, nmo, out);
    }
}

}  // namespace

void IWL::read_one(PSIO *psio, int itap, const char *label, double *ints, int ntri, int erase, int printflg,
                   std::string out) {
    read_one_impl(psio, itap, label, ints, ntri, erase, printflg, out);
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
int iwl_rdone(int itap, const char *label, double *ints, int ntri, int erase, int printflg, std::string out) {
    read_one_impl(_default_psio_lib_.get(), itap, label, ints, ntri, erase, printflg, out);
    return 1;
}
}
