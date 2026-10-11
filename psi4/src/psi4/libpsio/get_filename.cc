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
 \ingroup PSIO
 */

#include <cstdlib>
#include <string>
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/psi4-dec.h"

namespace psi {

std::string PSIO::get_filename(size_t unit, bool remove_namespace) const {
    std::string ns = "." + pid_;
    if (!(default_namespace_.empty() || remove_namespace)) ns += "." + default_namespace_;

    for (const char *kwdgrp : {"PSI", "DEFAULT"}) {
        for (int u : {static_cast<int>(unit), -1}) {
            const std::string &kval = filecfg_kwd(kwdgrp, "NAME", u);
            if (!kval.empty()) return kval + ns;
        }
    }

    // assume that the default has been provided already
    abort();
}
}  // namespace psi
