/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
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
 \ingroup PSIO
 */

#include "psi4/pragma.h"
#include <memory>
#include <cstdio>
#include <cstdlib>
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/psio.hpp"

namespace psi {

size_t PSIO::get_numvols(size_t unit) {
    std::string charnum;
    charnum = filecfg_kwd("PSI", "NVOLUME", unit);
    if (!charnum.empty()) return ((size_t)atoi(charnum.c_str()));
    charnum = filecfg_kwd("PSI", "NVOLUME", -1);
    if (!charnum.empty()) return ((size_t)atoi(charnum.c_str()));
    charnum = filecfg_kwd("DEFAULT", "NVOLUME", unit);
    if (!charnum.empty()) return ((size_t)atoi(charnum.c_str()));
    charnum = filecfg_kwd("DEFAULT", "NVOLUME", -1);
    if (!charnum.empty()) return ((size_t)atoi(charnum.c_str()));

    // assume that the default has been provided already
    abort();
}
}  // namespace psi
