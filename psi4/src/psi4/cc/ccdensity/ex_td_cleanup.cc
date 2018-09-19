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

void ex_td_cleanup(void)
{
  /*
   *  Clean out these files between computing x_ltd and x_rtd,
   *  because LHS and RHS eigenvectors change (get swapped).
   */

  psio_close(PSIF_CC_TMP,0);
  psio_close(PSIF_EOM_TMP,0);
  psio_close(PSIF_EOM_TMP0,0);
  psio_close(PSIF_EOM_TMP1,0);
  psio_close(PSIF_CC_GLG,0);
  psio_close(PSIF_CC_GL,0);
  psio_close(PSIF_CC_GR,0);

  psio_open(PSIF_CC_TMP,PSIO_OPEN_NEW);
  psio_open(PSIF_EOM_TMP,PSIO_OPEN_NEW);
  psio_open(PSIF_EOM_TMP0,PSIO_OPEN_NEW);
  psio_open(PSIF_EOM_TMP1,PSIO_OPEN_NEW);
  psio_open(PSIF_CC_GLG,PSIO_OPEN_NEW);
  psio_open(PSIF_CC_GL,PSIO_OPEN_NEW);
  psio_open(PSIF_CC_GR,PSIO_OPEN_NEW);

  return;
}

}} // namespace psi::ccdensity
