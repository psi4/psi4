/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2023 The Psi4 Developers.
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
#include <cstdio>
#include "psi4/libpsio/psio.h"
#include "iwl.h"
#include "iwl.hpp"

namespace psi {

/*!
** IWL_WRTONE()
**
** This function writes one-electron integrals.
**
**   itap       = tape to read ints from
**   label      = the PSIO label
**   ntri       = the size of the array (lower triangle)
**   onel_ints  = array to hold the one-electron integrals.
**
** David Sherrill, March 1995
** Revised by TDC, June 2001
** \ingroup IWL
*/
void IWL::write_one(PSIO *psio, int itap, const char *label, int ntri, double *onel_ints) {
    psio->open(itap, PSIO_OPEN_OLD);
    psio->write_entry(itap, label, (char *)onel_ints, ntri * sizeof(double));
    psio->close(itap, 1);
}
}
