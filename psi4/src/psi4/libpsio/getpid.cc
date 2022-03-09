/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2022 The Psi4 Developers.
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

/*! \defgroup PSIO libpsio: The PSI I/O Library */

/*!
 ** \file
 ** \ingroup PSIO
 */

#ifdef _MSC_VER
#include <process.h>
#define SYSTEM_GETPID ::_getpid
#else
#include <unistd.h>
#define SYSTEM_GETPID ::getpid
#endif
#include <cstring>
#include <cstdlib>
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/psi4-dec.h"
#include <sstream>

namespace psi {

std::string PSIO::getpid() {
    std::stringstream ss;

    if (psi::restart_id.empty()) {
        ss << SYSTEM_GETPID();
    } else
        ss << psi::restart_id;

    return ss.str();
}

std::string psio_getpid() {
    std::stringstream ss;

    if (psi::restart_id.empty()) {
        ss << SYSTEM_GETPID();
    } else
        ss << psi::restart_id;

    return ss.str();
}

}  // namespace psi
