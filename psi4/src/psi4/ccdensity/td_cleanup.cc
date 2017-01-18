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

void td_cleanup(void)
{
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

  if((params.ref==0) || (params.ref==1)) {
    free_block(moinfo.ltd);
    free_block(moinfo.rtd);
  }
  else if(params.ref==2) {
    free_block(moinfo.ltd_a);
    free_block(moinfo.ltd_b);
    free_block(moinfo.rtd_a);
    free_block(moinfo.rtd_b);
  }

  return;
}

}} // namespace psi::ccdensity
