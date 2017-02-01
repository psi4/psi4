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

/*! \defgroup PSIO libpsio: The PSI I/O Library */

/*!
 ** \file
 ** \ingroup PSIO
 */

#include <unistd.h>
#include <cstring>
#include <cstdlib>
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/psi4-dec.h"
#include <unistd.h>

namespace psi {

std::string PSIO::getpid(void) {
  std::stringstream ss;

  if (psi::restart_id.empty()) {
    pid_t pid = ::getpid();
    ss << pid;
  }
  else
    ss << psi::restart_id;

  return ss.str();
}

std::string psio_getpid(void) {
  std::stringstream ss;

  if (psi::restart_id.empty()) {
    pid_t pid = ::getpid();
    ss << pid;
  }
  else
    ss << psi::restart_id;

  return ss.str();
}

}
