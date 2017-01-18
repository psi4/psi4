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
 \ingroup PSIO
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/psi4-dec.h"

namespace psi {

void PSIO::get_volpath(unsigned int unit, unsigned int volume, char **path) {
  std::string kval;
  char volumeX[20];
  sprintf(volumeX, "VOLUME%u", volume+1);

  kval = filecfg_kwd("PSI", volumeX, unit);
  if (!kval.empty()) {
    *path = strdup(kval.c_str());
    return;
  }
  kval = filecfg_kwd("PSI", volumeX, -1);
  if (!kval.empty()) {
    *path = strdup(kval.c_str());
    return;
  }
  kval = filecfg_kwd("DEFAULT", volumeX, unit);
  if (!kval.empty()) {
    *path = strdup(kval.c_str());
    return;
  }
  kval = filecfg_kwd("DEFAULT", volumeX, -1);
  if (!kval.empty()) {
    *path = strdup(kval.c_str());
    return;
  }

  // assume default has been provided
  abort();
}

  int psio_get_volpath_default(unsigned int volume, char **path) {
    std::string kval;
    char volumeX[20];
    sprintf(volumeX, "VOLUME%u", volume+1);

    kval = _default_psio_lib_->filecfg_kwd("PSI", volumeX, -1);
    if (!kval.empty()) {
      *path = strdup(kval.c_str());
      return (1);
    }
    kval = _default_psio_lib_->filecfg_kwd("DEFAULT", volumeX, -1);
    if (!kval.empty()) {
      *path = strdup(kval.c_str());
      return (1);
    }

    // assume default has been provided
    abort();
  }

}
